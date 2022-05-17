#ifndef CAHNHILLIARDCOMBO_CPP
#define CAHNHILLIARDCOMBO_CPP

CHC::CHC(const CH_builder &p) : CH(p), epsilon_couplings(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                epsilon_couplingsSQR(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                    inverses(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                baremat(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                oldfieldFT(Field_Wrapper<complex<double>, complex<double>>(p)),
                                oldfieldNLW(Field_Wrapper<complex<double>, complex<double>>(p)),
                                InitWeight(Field_Wrapper<complex<double>, complex<double>>(p))
{
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        IdentityWeight<complex<double>> fw;
        oldfieldFT.add_method(fw, i);
        oldfieldNLW.add_method(fw, i);
    }
}

void CHC::setup_matrices() {

    matrix<double> identitymatrix(myp.number_of_fields,myp.number_of_fields);
    for(int i  = 0  ; i < myp.number_of_fields ; i++) {
        identitymatrix(i,i) = 1.;
    }

    int iter = 0;
    for(int k1 = 0  ; k1 <= myp.N1/2 ; k1 ++) {
        for(int k2 = k1 ; k2 <= myp.N2/2 ;  k2++) {
            cout << iter << endl;
            matrix<double> temp = create_D_mat_split(k1, k2);

            baremat[iter] = temp;
        
            matrix<double> to_invert =  identitymatrix -(1-alpha)*temp-alpha*dt*temp;
            to_invert.inverse();

            inverses[iter] =  to_invert;
            
            iter++;

        
        }
    }
}

void CHC::set_interaction(double val, int i, int j) {
    epsilon_couplings(i,j) = val;
    epsilon_couplings(j,i) = val;

}

void CHC::calculate_non_linear_weight(complex<double>**input) {

    // takes as an input a real space input of fields, produces as an output the fourier transformed weight function into transformed2

    int totp =  myp.get_total();

    for(int j = 0  ; j < totp ; j++) {
        weigs.calculated_reactions[0][j] = cons1 * CUB(input[0][j]) + cons2 * SQR(input[0][j]) + cons4;
    }

    for(int i  = 1 ; i < myp.number_of_fields; i++) {
        for(int j = 0 ; j < totp ; j++ )
            weigs.calculated_reactions[i][j] = 0.0;
    }

  //  outfunc(weigs.calculated_reactions[0],"test", myp);

    transformed2.Calculate_Results(weigs.calculated_reactions);

   // outfunc(transformed2.calculated_reactions[0], "testft", myp);
}

void CHC::calculate_non_linear_weightSQR(complex<double> **input)
{ 
    

    // takes as an input a real space input of fields, produces as an output the fourier transformed weight function into transformed2

    int totp = myp.get_total();


    double a = 5.;
    double x0 = 0.5;
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        
            for (int j = 0; j < totp; j++) {
                weigs.calculated_reactions[i][j] = 0.0;

                for (int k = 0; k < myp.number_of_fields; k++)
                {
                    //weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * input[i][j] * SQR(input[k][j]);
                    if(k!=i) {
                    if(epsilon_couplingsSQR(i,k) < 0.) {
                        weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * (a * SQR(1. / cosh(a * (input[i][j] - x0))) + a * SQR(1. / cosh(a * (input[i][j] - x0))) * tanh(a * (input[k][j] - x0)));
                    }
                    // else {
                    //     weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * input[i][j] * SQR(input[k][j]);
                    // }
                    }
                    // if (input[i][j] < 0.)
                    // {
                        
                    //     cout << input[i][j] << " " << input[k][j]<< endl;
                    //     cout << epsilon_couplingsSQR(i, k) << endl;
                    //     cout << epsilon_couplingsSQR(i, k) * (a * exp(a * input[i][j]) * (1. - 1. / (1. + exp(a * input[k][j])))) / SQR(1. + exp(a * input[i][j])) << endl;
                    //     pausel();
                    // }
                }
                
            }
    }

    // add a phi^4 term to stop over crowding for the other chemicals.
    // for (int i = 1; i < myp.number_of_fields; i++)
    // {
    //     for (int j = 0; j < totp; j++)
    //         weigs.calculated_reactions[i][j] += 0.25 * CUB(input[i][j]);

    // }

    for (int j = 0; j < totp; j++)
    {
        weigs.calculated_reactions[0][j] += cons1 * CUB(input[0][j]) + cons2 * SQR(input[0][j]) + cons4;
    }

    // for (int i = 0; i < myp.number_of_fields; i++)
    // {
    //     for (int j = 0; j < totp; j++)
    //     if (input[i][j] < 0.)
    //     {
    //         for(int k  = 0; k <  myp.number_of_fields ; k++) {
    //         cout << epsilon_couplingsSQR(i, k) << " " << input[i][j] << " " << input[k][j]<< endl;
    //         // cout << epsilon_couplingsSQR(i, k) << endl;
    //         // cout << epsilon_couplingsSQR(i, k) * (a * exp(a * input[i][j]) * (1. - 1. / (1. + exp(a * input[k][j])))) / SQR(1. + exp(a * input[i][j])) << endl;
    //         }
    //         cout << "weight: " <<  weigs.calculated_reactions[i][j];

    //         int effj = j % myp.N1 ;
    //         int effi = j / myp.N1;

    //         int effin1 = effi + 1 > myp.N1 ? 0 : effi + 1;
    //         int effin2 = effi - 1 < 0 ? myp.N1 : effi - 1;

    //         int effjn1 = effj + 1 > myp.N1 ? 0 : effj + 1;
    //         int effjn2 = effj - 1 < 0 ? myp.N1 : effj - 1;

    //         int totj1 = effin1 * myp.N1 + effj;

    //         int totj2 = effin2 * myp.N1 + effj;

    //         int totj3 = effi * myp.N1 + effjn1;

    //         int totj4 = effi * myp.N1 + effjn2;

    //         cout << "n1 weight: " << weigs.calculated_reactions[i][totj1];
    //         cout << "n2 weight: " << weigs.calculated_reactions[i][totj2];
    //         cout << "n3 weight: " << weigs.calculated_reactions[i][totj3];
    //         cout << "n4 weight: " << weigs.calculated_reactions[i][totj4];
    //         pausel();
    //     }
    // }

    // for(int j = 0  ; j < totp ; j++) {
    //     if(abs(input[1][j])>2.) {
    //         cout << input[1][j] << endl;
    //         for (int k = 0; k < myp.number_of_fields; k++)
    //         {
                
    //                 cout << "k: " << epsilon_couplingsSQR(1, k) * input[1][j] * SQR(input[k][j]) << ",";
    //                 cout << input[k][j] << ",";
    //                 cout << endl;
    //         }
    //         cout << endl;
    //         cout << weigs.calculated_reactions[1][j] << endl;
    //         cout << 0.2 * CUB(input[1][j]) << endl;
    //         pausel();
    //     }
    // }
    //  outfunc(weigs.calculated_reactions[0],"test", myp);

    transformed2.Calculate_Results(weigs.calculated_reactions);

    // outfunc(transformed2.calculated_reactions[0], "testft", myp);
}

void CHC::calculate_initial_weight() {
    transformed1.Calculate_Results(fields); //calculate FT of fields
    oldfieldFT.Calculate_Results(transformed1.calculated_reactions);

    calculate_non_linear_weight(fields);
    oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);


    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for (int i = 0; i < myp.N1; i++)
    {
        for (int j = 0; j < myp.N2; j++)
        {

            double k1, k2;
            if (i <= myp.N1 / 2)
            {
                k1 = i;
            }
            else
            {
                k1 = (i - myp.N1);
            }
            if (j <= myp.N2 / 2)
            {
                k2 = j;
            }
            else
            {
                k2 = (j - myp.N2);
            }

            //take absolute values of k1 and k2
            int k1a = abs(k1);
            int k2a = abs(k2);
            if(k2a < k1a) {
                k1a = k2a;
                k2a = abs(k1);
            }
            int rel = k2a - (k1a*(1 + k1a - 2*(1 + myp.N1/2)))/2.;



            double tempor = SQR(k1) + SQR(k2);


            vector1<complex<double > > v(nof);
            for(int k = 0  ; k < nof ; k++) {
                v[k] = transformed1.calculated_reactions[k][i * myp.N2 + j];
            }

            v = baremat[rel]*v;


            for(int k = 0 ; k < nof ; k++) {
                InitWeight.calculated_reactions[k][i * myp.N2 + j] = (v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
                //oldweightFT.calculated_reactions[k][i * myp.N2 + j] = (1- alpha) * (-v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
            }
            //upd1[i * myp.N2 + j] = -2*alpha*dt*diffusion*tempor*temp1*transformed2.calculated_reactions[]
            

            
            }
    }



}

void CHC::calculate_initial_weightSQR()
{


    transformed1.Calculate_Results(fields); // calculate FT of fields
    
    //cut off large wavenumbers

    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for(int fn = 0 ; fn < nof ; fn++) {
    for (int i = 200; i < myp.N1; i++)
    {
        for (int j = 200; j < myp.N2; j++)
        {
            transformed1.calculated_reactions[fn][i * myp.N2 + j]=0.;
        }
    }

    }
    reverse_transform.Calculate_Results(transformed1.calculated_reactions);
    set_field(reverse_transform.calculated_reactions);

    oldfieldFT.Calculate_Results(transformed1.calculated_reactions);

    calculate_non_linear_weightSQR(fields);
    oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);



    for (int i = 0; i < myp.N1; i++)
    {
        for (int j = 0; j < myp.N2; j++)
        {

            double k1, k2;
            if (i <= myp.N1 / 2)
            {
                k1 = i;
            }
            else
            {
                k1 = (i - myp.N1);
            }
            if (j <= myp.N2 / 2)
            {
                k2 = j;
            }
            else
            {
                k2 = (j - myp.N2);
            }

            // take absolute values of k1 and k2
            int k1a = abs(k1);
            int k2a = abs(k2);
            if (k2a < k1a)
            {
                k1a = k2a;
                k2a = abs(k1);
            }
            int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

            double tempor = SQR(k1) + SQR(k2);

            vector1<complex<double>> v(nof);
            for (int k = 0; k < nof; k++)
            {
                v[k] = transformed1.calculated_reactions[k][i * myp.N2 + j];
            }

            v = baremat[rel] * v;

            for (int k = 0; k < nof; k++)
            {
                InitWeight.calculated_reactions[k][i * myp.N2 + j] = (v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
                // oldweightFT.calculated_reactions[k][i * myp.N2 + j] = (1- alpha) * (-v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
            }
            // upd1[i * myp.N2 + j] = -2*alpha*dt*diffusion*tempor*temp1*transformed2.calculated_reactions[]
        }
    }
}

void CHC::Update() {

    // string fieldss = "fields";
    // outfunc(fields[0], fieldss, myp);

    chems.Calculate_Results(fields); //calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);
   



    transformed3.Calculate_Results(chems.calculated_reactions);
    // string ftschem = "ftchem";
    // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //     pausel();

    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();

    transformed1.Calculate_Results(fields); // calculate FT of fields
    calculate_non_linear_weight(fields);
    
    
    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for (int i = 0; i < myp.N1; i++)
    {
        for (int j = 0; j < myp.N2; j++)
        {

            double k1, k2;
            if (i <= myp.N1 / 2)
            {
                k1 = i;
            }
            else
            {
                k1 = (i - myp.N1);
            }
            if (j <= myp.N2 / 2)
            {
                k2 = j;
            }
            else
            {
                k2 = (j - myp.N2);
            }

            //take absolute values of k1 and k2
            int k1a = abs(k1);
            int k2a = abs(k2);
            if (k2a < k1a)
            {
                k1a = k2a;
                k2a = abs(k1);
            }
            int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

            double tempor = SQR(k1) + SQR(k2);

            // vector1<complex<double> > v(nof);


            double fac =  diffusion * tempor * temp1 ;
            vector1<complex<double> > v(nof);
            vector1<complex<double> > v2(nof);
            vector1<complex<double> > v3(nof);

            for (int k = 0; k < nof; k++)
            {
                v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j] ;// oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
            }

            v2 = baremat[rel] * v2;


            

            for (int k = 0; k < nof; k++)
            {
                v[k] = -((1-alpha) + (alpha) * dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) 
                + transformed1.calculated_reactions[k][i * myp.N2 + j]//oldfieldFT.calculated_reactions[k][i * myp.N2 + j] 
                - (1-alpha)*v2[k] 
                + ((1-alpha) /* -  0.5*(alpha) * dt */) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] 
                - 0.0*2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] //zeroed because of the initial conditions we set
                + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
            }


            v3 = inverses[rel]*v;
            // if (rel == 131839)
            // {
            //     cout << inverses[rel] << endl;
            //     cout << nof << endl;
            //     for (int k = 0; k < nof; k++)
            //     {

            //         cout << "before: " << transformed1.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "1: " << -(2 * (1 - alpha) - (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) << endl;
            //         cout << "2: " << +oldfieldFT.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "3: " << -(1 - alpha) * v2[k] << endl;
            //         cout << "4: " << +(2 * (1 - alpha) + (alpha)*dt) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "5: " << +2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "6: " << +2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
            //        // cout << "7: " << 2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "after: " << v3[k] << endl;
            //         cout << endl;
            //         cout << endl;
            //     }

            //     cout << v << endl;
            //     cout << v3 << endl;
            //     pausel();
            // }
            for (int k = 0; k < nof; k++)
            {
                rules.calculated_reactions[k][i * myp.N2 + j] = v3[k];
            }
        }
    }

    oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
    oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

    reverse_transform.Calculate_Results(rules.calculated_reactions);

    GetMaximas(fields,myp);
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;

    set_field(reverse_transform.calculated_reactions);


}

void CHC::UpdateSQR()
{
    GetMaximas(fields, myp);
    GetMinimas(fields, myp);
    // string fieldss = "fields";
    // outfunc(fields[0], fieldss, myp);

    chems.Calculate_Results(fields); // calculate chemistry


    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);

    // string ftschem = "ftchem";
    // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //     pausel();

    // GetMaximas(transformed3.calculated_reactions,schem,myp);



    transformed1.Calculate_Results(fields); // calculate FT of fields
    calculate_non_linear_weightSQR(fields);

    // cout << transformed2.calculated_reactions[0][512*1024+512] << endl;

    // cout << oldfieldNLW.calculated_reactions[0][512 * 1024 + 512] << endl;

    // cout << transformed2.calculated_reactions[1][512 * 1024 + 512] << endl;

    // cout << oldfieldNLW.calculated_reactions[1][512 * 1024 + 512] << endl;

    // cout << transformed2.calculated_reactions[2][512 * 1024 + 512] << endl;

    // cout << oldfieldNLW.calculated_reactions[2][512 * 1024 + 512] << endl;

    // cout << transformed2.calculated_reactions[3][512 * 1024 + 512] << endl;

    // cout << oldfieldNLW.calculated_reactions[3][512 * 1024 + 512] << endl;

    // pausel();

    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for (int i = 0; i < myp.N1; i++)
    {
        for (int j = 0; j < myp.N2; j++)
        {

            double k1, k2;
            if (i <= myp.N1 / 2)
            {
                k1 = i;
            }
            else
            {
                k1 = (i - myp.N1);
            }
            if (j <= myp.N2 / 2)
            {
                k2 = j;
            }
            else
            {
                k2 = (j - myp.N2);
            }

            // take absolute values of k1 and k2
            int k1a = abs(k1);
            int k2a = abs(k2);
            if (k2a < k1a)
            {
                k1a = k2a;
                k2a = abs(k1);
            }
            int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

            double tempor = SQR(k1) + SQR(k2);

            // vector1<complex<double> > v(nof);

            double fac = diffusion * tempor * temp1;
            vector1<complex<double>> v(nof);
            vector1<complex<double>> v2(nof);
            vector1<complex<double>> v3(nof);

            for (int k = 0; k < nof; k++)
            {
                v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
            }

            v2 = baremat[rel] * v2;

            for (int k = 0; k < nof; k++)
            {
                v[k] = -((1 - alpha) + (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                       // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
                       - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] /* - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j]  */// zeroed because of the initial conditions we set
                       + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
            }

            v3 = inverses[rel] * v;
            // if (rel == 131839)
            // {
            //     cout << inverses[rel] << endl;
            //     cout << nof << endl;
            //     for (int k = 0; k < nof; k++)
            //     {

            //         cout << "before: " << transformed1.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "1: " << -(2 * (1 - alpha) - (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) << endl;
            //         cout << "2: " << +oldfieldFT.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "3: " << -(1 - alpha) * v2[k] << endl;
            //         cout << "4: " << +(2 * (1 - alpha) + (alpha)*dt) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "5: " << +2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "6: " << +2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
            //        // cout << "7: " << 2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
            //         cout << "after: " << v3[k] << endl;
            //         cout << endl;
            //         cout << endl;
            //     }

            //     cout << v << endl;
            //     cout << v3 << endl;
            //     pausel();
            // }
            for (int k = 0; k < nof; k++)
            {
                rules.calculated_reactions[k][i * myp.N2 + j] = v3[k];
            }
        }
    }



    reverse_transform.Calculate_Results(rules.calculated_reactions);


    bool checker = false;

    if(checker)
    for (int ifn = 0; ifn < myp.number_of_fields; ifn++)
    {
        for(int i1 = 0 ; i1 < myp.N1 ; i1++)
        for(int j1 = 0 ; j1 < myp.N2 ; j1++)
        if(reverse_transform.calculated_reactions[ifn][i1*myp.N1+j1]<0.0)
            {
                cout << ifn << endl;
                cout << fields[ifn][i1*myp.N1+j1] << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1*myp.N1+j1] << endl;

                Field_Wrapper<complex<double>, complex<double>> store1(myp);
                Field_Wrapper<complex<double>, complex<double>> store2(myp);
                Field_Wrapper<complex<double>, complex<double>> store3(myp);
                Field_Wrapper<complex<double>, complex<double>> store4(myp);
                Field_Wrapper<complex<double>, complex<double>> store5(myp);
                Field_Wrapper<complex<double>, complex<double>> store6(myp);
                for (int i = 0; i < myp.N1; i++)
                {
                    for (int j = 0; j < myp.N2; j++)
                    {

                        double k1, k2;
                        if (i <= myp.N1 / 2)
                        {
                            k1 = i;
                        }
                        else
                        {
                            k1 = (i - myp.N1);
                        }
                        if (j <= myp.N2 / 2)
                        {
                            k2 = j;
                        }
                        else
                        {
                            k2 = (j - myp.N2);
                        }

                        // take absolute values of k1 and k2
                        int k1a = abs(k1);
                        int k2a = abs(k2);
                        if (k2a < k1a)
                        {
                            k1a = k2a;
                            k2a = abs(k1);
                        }
                        int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

                        double tempor = SQR(k1) + SQR(k2);
                        double fac = diffusion * tempor * temp1;
                        vector1<complex<double>> v(nof);
                        vector1<complex<double>> v2(nof);
                        vector1<complex<double>> v3(nof);

                        for (int k = 0; k < nof; k++)
                        {
                            v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
                        }

                        v2 = baremat[rel] * v2;

                        vector1<complex<double>> vtemp1(nof);
                        vector1<complex<double>> vtemp2(nof);
                        vector1<complex<double>> vtemp3(nof);
                        vector1<complex<double>> vtemp4(nof);
                        vector1<complex<double>> vtemp5(nof);
                        vector1<complex<double>> vtemp6(nof);

                        for (int k = 0; k < nof; k++)
                        {
                            vtemp1[k] = -((1 - alpha) + (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]);
                            vtemp2[k] = + transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
                            vtemp3[k] = - (1 - alpha) * v2[k];
                            vtemp4[k]= + ((1 - alpha)) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j];
                            vtemp5[k] = + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
                            vtemp6[k] = -((1 - alpha) + (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                              // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
                                               - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] /* - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j]  */ // zeroed because of the initial conditions we set
                                               + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
                        }
                        v3 = inverses[rel] * vtemp1;

                        for (int k = 0; k < nof; k++)
                        {
                            store1.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }
                        v3 = inverses[rel] * vtemp2;

                        for (int k = 0; k < nof; k++)
                        {
                            store2.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }
                        v3 = inverses[rel] * vtemp3;

                        for (int k = 0; k < nof; k++)
                        {
                            store3.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }
                        v3 = inverses[rel] * vtemp4;

                        for (int k = 0; k < nof; k++)
                        {
                            store4.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }
                        v3 = inverses[rel] * vtemp5;

                        for (int k = 0; k < nof; k++)
                        {
                            store5.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }

                        v3 = inverses[rel] * vtemp6;

                        for (int k = 0; k < nof; k++)
                        {
                            store6.calculated_reactions[k][i * myp.N2 + j] = v3[k];
                        }
                    }
                }

                double k1, k2;
                if (i1 <= myp.N1 / 2)
                {
                    k1 = i1;
                }
                else
                {
                    k1 = (i1 - myp.N1);
                }
                if (j1 <= myp.N2 / 2)
                {
                    k2 = j1;
                }
                else
                {
                    k2 = (j1 - myp.N2);
                }

                // take absolute values of k1 and k2
                int k1a = abs(k1);
                int k2a = abs(k2);
                if (k2a < k1a)
                {
                    k1a = k2a;
                    k2a = abs(k1);
                }
                int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

                double a = 10.;

                for(int k  = 0 ; k < nof ; k++) {
                if (epsilon_couplingsSQR(ifn, k) < 0.)
                {
                    int j = i1*myp.N1+j1;
                    cout << "negative weight: ";
                    cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
                    cout << endl;
                }
                else
                {
                    int j = i1 * myp.N1 + j1;
                    cout << "positive weight: ";
                    cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
                    cout << endl;
                }
                }

                for (int k = 0; k < nof; k++)
                {
                    if (epsilon_couplingsSQR(ifn, k) < 0.)
                    {
                        int j = i1 * myp.N1 + j1 -1 ;
                        cout << "negative weight n1: ";
                        cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
                        cout << endl;
                    }
                    else
                    {
                        int j = i1 * myp.N1 + j1 -1 ;
                        cout << "positive weight n1: ";
                        cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
                        cout << endl;
                    }
                }
                for (int k = 0; k < nof; k++)
                {
                    if (epsilon_couplingsSQR(ifn, k) < 0.)
                    {
                        int j = i1 * myp.N1 + j1 + 1;
                        cout << "negative weight n2: ";
                        cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
                        cout << endl;
                    }
                    else
                    {
                        int j = i1 * myp.N1 + j1 + 1;
                        cout << "positive weight n2: ";
                        cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
                        cout << endl;
                    }
                }
                for (int k = 0; k < nof; k++)
                {
                    if (epsilon_couplingsSQR(ifn, k) < 0.)
                    {
                        int j = (i1 + 1) * myp.N1 + j1;
                        cout << "negative weight n3: ";
                        cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
                        cout << endl;
                    }
                    else
                    {
                        int j = (i1 + 1) * myp.N1 + j1 + 1;
                        cout << "positive weight n3: ";
                        cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
                        cout << endl;
                    }
                }

                for (int k = 0; k < nof; k++)
                {
                    if (epsilon_couplingsSQR(ifn, k) < 0.)
                    {
                        int j = (i1 - 1) * myp.N1 + j1 + 1;
                        cout << "negative weight n4: ";
                        cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
                        cout << endl;
                    }
                    else
                    {
                        int j = (i1 - 1) * myp.N1 + j1 + 1;
                        cout << "positive weight n4: ";
                        cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
                        cout << endl;
                    }
                }

                cout << "inv" << endl;
                cout << inverses[rel] << endl;



                reverse_transform.Calculate_Results(transformed1.calculated_reactions);
                cout << "i1: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store1.calculated_reactions);
                cout << "t1: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store2.calculated_reactions);
                cout << "t2: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store3.calculated_reactions);
                cout << "t3: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store4.calculated_reactions);
                cout << "t4: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store5.calculated_reactions);
                cout << "t5: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                reverse_transform.Calculate_Results(store6.calculated_reactions);
                cout << "t6: " << endl;
                cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

                pausel();
            }
        }

        oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
        oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

        GetMaximas(fields, myp);
        reverse_transform.GetMaximas();
        reverse_transform.GetMaximasIndex();
        reverse_transform.GetMinimas();
        reverse_transform.GetMinimasIndex();
        cout << endl;

        set_field(reverse_transform.calculated_reactions);
}
#endif /* CAHNHILLIARDCOMBO_CPP */
