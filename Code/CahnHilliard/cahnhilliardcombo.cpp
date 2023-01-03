#ifndef CAHNHILLIARDCOMBO_CPP
#define CAHNHILLIARDCOMBO_CPP

CHC::CHC(const CH_builder &p) : CH(p), phase_separators(vector1<bool>(p.number_of_fields)), epsilon_couplings(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                epsilon_couplingsSQR(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                inverses(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                baremat(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                oldfieldFT(Field_Wrapper<complex<double>, complex<double>>(p)),
                                oldfieldNLW(Field_Wrapper<complex<double>, complex<double>>(p)),
                                InitWeight(Field_Wrapper<complex<double>, complex<double>>(p)), OLD(vector<Eigen::VectorXcd>())
{
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        IdentityWeight<complex<double>> fw;
        oldfieldFT.add_method(fw, i);
        oldfieldNLW.add_method(fw, i);
    }
    M=1;

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



    for(int i  = 0 ; i < myp.number_of_fields; i++) {
        if(phase_separators[i]) {
        for (int j = 0; j < totp; j++)
        {
            weigs.calculated_reactions[i][j] = cons1 * CUB(input[i][j]) + cons2 * SQR(input[i][j]) + cons4;
        }
        }
        else{
        for(int j = 0 ; j < totp ; j++ )
            weigs.calculated_reactions[i][j] = 0.0;
        }
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
    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    //cut off form
    int cut_off = 100;
    for (int fn = 0; fn < nof; fn++)
    {
        for (int i = 0; i < myp.N1; i++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                    int k1,k2;
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
                    double tempor = SQR(k1) + SQR(k2);
                    if(tempor>cut_off)
                    transformed1.calculated_reactions[fn][i * myp.N2 + j] = 0.; // cut off the high frequency modes
            }
        }
    }
    

//    //1/f form
//     for (int fn = 0; fn < nof; fn++)
//     {
//         for (int i = 0; i < myp.N1; i++)
//         {
//                 for (int j = 0; j < myp.N2; j++)
//                 {

//                     double k1, k2;
//                     if (i <= myp.N1 / 2)
//                     {
//                     k1 = i;
//                 }
//                 else
//                 {
//                     k1 = (i - myp.N1);
//                 }
//                 if (j <= myp.N2 / 2)
//                 {
//                     k2 = j;
//                 }
//                 else
//                 {
//                     k2 = (j - myp.N2);
//                 }

//                 // take absolute values of k1 and k2
//                 int k1a = abs(k1);
//                 int k2a = abs(k2);

//                 double tempor = SQR(k1) + SQR(k2);
//                 if (tempor < 1.E-10)
//                 {
//                 }
//                 else if(tempor < 1600. && tempor >1.E-10)
//                     transformed1.calculated_reactions[fn][i * myp.N2 + j] *= 10./sqrt(tempor);

//                 else {
//                     transformed1.calculated_reactions[fn][i * myp.N2 + j]=0.;
//                 }
//                 // upd1[i * myp.N2 + j] = -2*alpha*dt*diffusion*tempor*temp1*transformed2.calculated_reactions[]
//             }
//         }
//     }


    reverse_transform.Calculate_Results(transformed1.calculated_reactions);

    //now multiply throught


    reverse_transform.rescale();

    set_field(reverse_transform.calculated_reactions);

    transformed1.Calculate_Results(fields); // calculate FT of fields

    oldfieldFT.Calculate_Results(transformed1.calculated_reactions);

    calculate_non_linear_weight(fields);
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

    transformed1.Calculate_Results(fields); // calculate FT of fields

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

template<class Q>
void CHC::UpdateNoise(Q &func, GenNoise<complex<double>> &mynoise, vector1<double> &str)
{

    // string fieldss = "fields";
    // outfunc(fields[0], fieldss, myp);

    chems.Calculate_Results(fields); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);
    // string ftschem = "ftchem";
    // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //     pausel();

    mynoise.GenFields(func, str, 1.0 / 2.0);

    transformed3.Add_Noise(mynoise);


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
                       - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] // zeroed because of the initial conditions we set
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

    oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
    oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

    reverse_transform.Calculate_Results(rules.calculated_reactions);

    GetMaximas(fields, myp);
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


Eigen::VectorXcd CHC::CalculateWeightPartial(complex<double> *input) {
        // cout << "beginning weight calculation" << endl;
        int totp = myp.get_total();
        int nof = myp.number_of_fields;
        Eigen::VectorXcd x(totp);

        for (int i = 0; i < totp; i++)
        {
            x[i] = 0.0;
        }
        int N = myp.N1;
        double dtalpha = pow(dt, alpha);
        int fn = 0;
        // transformed1.Calculate_Results(input); //find the FT of the fields
        
        complex<double> *in;
        complex<double> *out;
        fftw_plan p;

        in = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));

        out = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));

        for (int j = 0; j < totp; j++)
        {
            in[j] = cons1 * CUB(input[j]) + cons2 * SQR(input[j]) + cons4;
        }


            //  outfunc(weigs.calculated_reactions[0],"test", myp);
        p = fftw_plan_dft_2d(myp.N1, myp.N2, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

        // cout << "created plan" << endl;
        fftw_execute(p);
        for (int i1 = 0; i1 < N; i1++)
        {
            for (int j1 = 0; j1 < N; j1++)
            {

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
                // x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * cons3s * temp1 * (SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn][i1 * N + j1] - dtalpha*diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn][i1 * N + j1];

                x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (out[i1 * N + j1] / (double)myp.N1);

                // for(int fn2 = 1 ; fn2 < nof ; fn2++) {
                //     x[fn * totp + i1 * N + j1] += -epsilon_couplings(fn, fn2) * dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn2][i1 * N + j1];
                // }
            }
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);

        return x;
}


Eigen::VectorXcd CHC::CalculateWeight(complex<double> **input) {
        // cout << "beginning weight calculation" << endl;
        int totp = myp.get_total();
        int nof = myp.number_of_fields;
        Eigen::VectorXcd x(nof * totp);

        for (int i = 0; i < nof * totp; i++)
        {
            x[i] = 0.0;
        }
        int N = myp.N1;
        double dtalpha = pow(dt, alpha);
        int fn = 0;
        // transformed1.Calculate_Results(input); //find the FT of the fields
        
        complex<double> *in;
        complex<double> *out;
        fftw_plan p;

        in = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));

        out = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));

        for (int j = 0; j < totp; j++)
        {
            in[j] = cons1 * CUB(input[0][j]) + cons2 * SQR(input[0][j]) + cons4;
        }


            //  outfunc(weigs.calculated_reactions[0],"test", myp);
        p = fftw_plan_dft_2d(myp.N1, myp.N2, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

        // cout << "created plan" << endl;
        fftw_execute(p);
        
        // transformed2.Calculate_Results(weigs.calculated_reactions);


        //for(int fn = 0 ; fn < myp.number_of_fields ; fn++) {
            
            for (int i1 = 0; i1 < N; i1++)
            {
                for (int j1 = 0; j1 < N; j1++)
                {

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
                    // x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * cons3s * temp1 * (SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn][i1 * N + j1] - dtalpha*diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn][i1 * N + j1];

                    x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (out[i1 * N + j1]/(double)myp.N1);

                    // for(int fn2 = 1 ; fn2 < nof ; fn2++) {
                    //     x[fn * totp + i1 * N + j1] += -epsilon_couplings(fn, fn2) * dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * transformed1.calculated_reactions[fn2][i1 * N + j1];
                    // }

                }
            }
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);
            

            Field_Wrapper<complex<double>, complex<double>> tempf(myp);

            for(int i = 0  ; i < myp.number_of_fields ; i++) {
                FourierWeightForward2D fw;
                tempf.add_method(fw,i);
            }
            //overwrite the weight fourier transform with the field fourier transform
            tempf.Calculate_Results(input);

            
            for (int i1 = 0; i1 < N; i1++)
            {
                for (int j1 = 0; j1 < N; j1++)
                {

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
                    x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * cons3 * temp1 * (SQR(k1) + SQR(k2)) * tempf.calculated_reactions[fn][i1 * N + j1] - dtalpha * diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * tempf.calculated_reactions[fn][i1 * N + j1];

                    //x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * transformed2.calculated_reactions[fn][i1 * N + j1];

                    for (int fn2 = 1; fn2 < nof; fn2++)
                    {
                        x[fn * totp + i1 * N + j1] += -epsilon_couplings(fn, fn2) * dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * tempf.calculated_reactions[fn2][i1 * N + j1];
                    }
                }
            }

            for(int fn = 1 ; fn < nof ; fn++)
            for (int i1 = 0; i1 < N; i1++)
            {
                for (int j1 = 0; j1 < N; j1++)
                {

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
                    
                    //x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * transformed2.calculated_reactions[fn][i1 * N + j1];

                    for (int fn2 = 0; fn2 < nof; fn2++)
                    {
                        if(fn!=fn2)
                            x[fn * totp + i1 * N + j1] += -epsilon_couplings(fn, fn2) * dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * tempf.calculated_reactions[fn2][i1 * N + j1];
                        else
                            x[fn * totp + i1 * N + j1] += -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * tempf.calculated_reactions[fn2][i1 * N + j1];
                    }
                }
            }

            return x;
}

Eigen::VectorXcd CHC::CalculateRHS(complex<double> **cnp1, complex<double> **cn, complex<double> **Rn) //input are the already fourier transformed versions
{

    reverse_transform.Calculate_Results(cnp1);

    //chems.Calculate_Results(reverse_transform.calculated_reactions); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    //transformed3.Calculate_Results(chems.calculated_reactions);

    Eigen::VectorXcd wn = CalculateWeight(reverse_transform.calculated_reactions);

    //the input is the fourier transformed version of input


    double c = 1.0;
    Eigen::VectorXcd totalwn = c * wn;
    int totp = myp.get_total();
    int nof = myp.number_of_fields;



    for(int i = 0  ; i < M ; i++) {
        c = (1-(1+(1-alpha))/(double)(i+1))*c;
        //cout << c << endl;
        totalwn += c*OLD[i];

    }
    // for(int i = 0 ; i < nof ; i++) {
    //     for(int j = 0  ; j < totp ;j++) {
    //         if (abs(totalwn[i * totp + j])>1E6 ) {
    //             cout << i << " " << j << endl;
    //             pausel();
    //         }
    //     }
    // }
    //calculate
    for(int i = 0 ; i < nof ; i++) {
        for(int j = 0  ; j < totp ;j++) {
            totalwn[i*totp+j] += cn[i][j] - cnp1[i][j] + dt*Rn[i][j];
        }
    }

    return totalwn;
}

Eigen::VectorXd CHC::CalculateRHS_real(complex<double> **cnp1, complex<double> **cn, complex<double> **Rn) // input are the already fourier transformed versions
{

    //THIS FUNCTION EXCLUDES THE FINAL CONTRIBUTION OF CNP1, SUBTRACT IT AFTER A FOURIER TRANSFORM

    Eigen::VectorXcd wn = CalculateWeight(cnp1);

    // the input is the fourier transformed version of input

    int N = myp.N1 * myp.N2;
    int totN = myp.number_of_fields * N;

    double c = 1.0;
    Eigen::VectorXcd totalwn = c * wn;
    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for (int i = 0; i < M; i++)
    {
        c = (1 - (1 + (1 - alpha)) / (double)(i + 1)) * c;
        // cout << c << endl;
        totalwn += c * OLD[i];
    }

    for (int i = 0; i < nof; i++)
    {
        for (int j = 0; j < totp; j++)
        {
            totalwn[i * totp + j] += cn[i][j] + dt * Rn[i][j];
        }
    }

    Field_Wrapper<complex<double>, complex<double>> storex(myp);


    for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    {
        for (int p = 0; p < myp.N1 * myp.N2; p++)
        {
            storex.calculated_reactions[fieldno1][p] = totalwn[fieldno1 * N + p];
        }
    }



    
    reverse_transform.Calculate_Results(storex.calculated_reactions);
    Eigen::VectorXd F2(totN);
    for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    {
        for (int p = 0; p < myp.N1 * myp.N2; p++)
        {

            F2[fieldno1 * N + p] = -cnp1[fieldno1][p].real() + reverse_transform.calculated_reactions[fieldno1][p].real();
        }
    }

    return F2;
}

Eigen::VectorXcd CHC::SolveLinearProblem(Eigen::SparseMatrix<complex<double> > &unchangedJ, complex<double> **cn, complex<double> **Rn)
{
    int N = myp.N1 * myp.N2;
    int totN = myp.number_of_fields * N;
    double c = 1.0;
    Eigen::VectorXcd totalwn(totN);
    for(int i = 0 ; i < totN ; i++) {
        totalwn[i] = 0.0;
    }
    
    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    for (int i = 0; i < M; i++)
    {
        c = (1 - (1 + (1 - alpha)) / (double)(i + 1)) * c;
        // cout << c << endl;
        totalwn += c * OLD[i];
    }

    for (int i = 0; i < nof; i++)
    {
        for (int j = 0; j < totp; j++)
        {
            totalwn[i * totp + j] += cn[i][j] + dt * Rn[i][j];
        }
    }

    // Field_Wrapper<complex<double>, complex<double>> storex(myp);
    // for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    // {
    //     for (int p = 0; p < myp.N1 * myp.N2; p++)
    //     {
    //         storex.calculated_reactions[fieldno1][p] = totalwn[fieldno1 * N + p];
    //     }
    // }

    // reverse_transform.Calculate_Results(storex.calculated_reactions);
    // Eigen::VectorXd F2(totN);
    // for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    // {
    //     for (int p = 0; p < myp.N1 * myp.N2; p++)
    //     {

    //         F2[fieldno1 * N + p] = reverse_transform.calculated_reactions[fieldno1][p].real();
    //     }
    // }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<complex<double>>> solver;
    solver.compute(unchangedJ);
    Eigen::VectorXcd dxv = solver.solve(-totalwn);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;

    return dxv;
    //cg.compute(J);
}

struct cdi {
    int index;
    double Score;
};
typedef Eigen::Triplet<complex<double>> Trip;
typedef Eigen::Triplet<double> TripD;

Eigen::SparseMatrix<complex<double>> CHC::CalculateJ(complex<double> **input, bool show = false) //the input is Fourier transformed already
{
    

    int field_no = myp.number_of_fields;
    // for the first field fields
    int totN = myp.N1*myp.N2;
    //for each individual field
    // for (int i1 = 0; i1 < myp.N1; i1++) {
    //     for (int j1 = 0; j1 < myp.N2; j1++) {
        
    //     }
    // }

    vector<Trip> coeffs;

    double dtalpha = pow(dt, alpha);
    
    reverse_transform.Calculate_Results(input);


    complex<double> meanval = 0.0;
    for (int i1 = 0; i1 < myp.N1; i1++)
    {
        for (int j1 = 0; j1 < myp.N2; j1++)
        {
            int id = myp.N1 * i1 + j1;
            meanval += (3. * cons1 * SQR(reverse_transform.calculated_reactions[0][id]) + 2. * cons2 * (reverse_transform.calculated_reactions[0][id]) + cons3);
        }
    }

    //cout << meanval/(double)(totN) << endl;
    meanval /= (double)(totN);

    complex<double> *in;
    complex<double> *out;
    fftw_plan p, p2;

    in = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));

    out = (complex<double> *)fftw_malloc(myp.N1 * myp.N2 * sizeof(complex<double>));
    for (int i = 0; i < myp.N1; i++)
    {
        for (int j = 0; j < myp.N2; j++)
        {
            int id = myp.N1 * i + j;
            in[id] = (3. * cons1 * SQR(reverse_transform.calculated_reactions[0][id]) + 2. * cons2 * (reverse_transform.calculated_reactions[0][id]) + cons3);
            out[id] = 0.;
        }
    }

    p = fftw_plan_dft_2d(myp.N1, myp.N2, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

    // cout << "created plan" << endl;
    fftw_execute(p);


    // vector<int> indexes;
    vector<cdi> allt;
    for(int i = 0  ; i < totN ; i++) {
        allt.push_back({i,abs(out[i])});
        // indexes.push_back(i);
    }

    std::sort(allt.begin(), allt.end(),
              [](const auto &i, const auto &j)
              { return i.Score > j.Score; });

    // for(int i = 0 ; i < 100 ; i++)
    // cout << allt[i].index << endl;


    // pausel();
    // int iter = 0;
    // double m = abs(out[0]);
    // for(int i = 0  ; i < totN ; i++) {
    //     if(abs(out[i])>0.001*m) {
    //         iter++;
    //     }
    // }
    // cout << iter << endl;
    // pausel();
    // cout << out[0] / (double)(totN) << endl;
    // vector<complex<double> > c;
    // for(int i = 0 ; i < totN ; i++) {
    //     c.push_back(out[i]/(double)(totN));
    // }

    // int nof = 1;
    // CH_builder pg;
    // pg.number_of_fields = nof;
    // pg.N1 = myp.N1;
    // pg.N2 = myp.N2;
    // outfunc(out, "StartData", pg);
    // // cout << "done 2" << endl;
    // pausel();

    // meanval /= (double)(myp.N1*myp.N2);

    //calculate the FT of the weights
    int fmodes = 1024;
    // for each individual field
    for (int i1 = 0; i1 < myp.N1; i1++)
    {
        for (int j1 = 0; j1 < myp.N2; j1++)
        {

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
            int id = myp.N1 * i1 + j1;
            complex<double> w = -1.0 - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (meanval) - dtalpha * diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2));

            
            // complex<double> w = -1.0 - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (3. * cons1 * SQR(reverse_transform.calculated_reactions[0][id]) + 2. * cons2 * (reverse_transform.calculated_reactions[0][id]) + cons3) - dtalpha * diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2));
            // complex<double> w = -1.0 - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (3 * cons1 * SQR(input[0][id]) + 2 * cons2 * (input[0][id]) + cons3) - dtalpha * diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2));

            // if (show == true && i1 == 0 && j1 == 1)
            // {
            //     cout << "start" << endl;
            //     cout << cons1 << endl;
            //     cout << cons2 << endl;
            //     cout << reverse_transform.calculated_reactions[0][id] << endl;
            //     cout << SQR(reverse_transform.calculated_reactions[0][id]) << endl;
            //     cout << -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (3 * cons1 * SQR(reverse_transform.calculated_reactions[0][id])) << endl;
            //     cout << -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * 2 * cons2 * (reverse_transform.calculated_reactions[0][id]) << endl;
            //     cout << -dtalpha * diffusion *temp1 * (SQR(k1) + SQR(k2))  * cons3 << endl;
            //     cout << -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * (3 * cons1 * SQR(reverse_transform.calculated_reactions[0][id]) + 2 * cons2 * (reverse_transform.calculated_reactions[0][id]) + cons3) << endl;
            //     cout <<  -dtalpha * diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) << endl;
            //     cout << w << endl;
            //     cout << "done" << endl;
            // }
            coeffs.push_back(Trip(id, id, w));

            // for(int j1 = 0  ; j1 < fmodes ; j1++) {
            //     int j = allt[j1].index;
            //     //if(j!= id) {
            //         int indx =  j + id;
            //         if(indx >= totN) indx -= totN; 
            //         complex<double> w2 = - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * conj((out[j]/(double)(totN)));
            //         //if(abs(out[indx]) > 0.01*abs(out[0])) { 
            //         // coeffs.push_back(Trip(id, indx, w2));
            //         if(indx != id) coeffs.push_back(Trip(indx, id, w2));
            //         // if (indx != id)
            //         //     coeffs.push_back(Trip(id, indx, w2));
            //         // // coeffs.push_back(Trip(j, id, w2));
            //         //}                    
            //    // }
            // // }
            // for(int j = 0  ; j < totN ; j++) {
            //     if(j!=id) {
            //         int indx =  j - id;
            //         if(indx < 0) indx += totN;
            //         complex<double> w2 = - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * conj((out[indx]/(double)(totN)));
            //         //if(abs(out[indx]) > 0.01*abs(out[0])) {
            //            if( id == 3 && j ==32) {
            //             cout << indx << " " << j << endl;
            //             pausel();
            //            }
            //            coeffs.push_back(Trip(id, j, w2));
            //         //    coeffs.push_back(Trip(j,id, w2));
            //         //}
            //     }
            // }
            
            if(show)
            for (int i2 = 0; i2 < myp.N1; i2++)
            {
                for(int j2 = 0  ; j2 < myp.N2 ; j2++) {
                

                    int indx1 = i2 - i1;
                    int indx2 = j2 - j1;
                    if (indx1 < 0)
                        indx1 += myp.N1;
                    if (indx2 < 0)
                        indx2 += myp.N2;

                    int findx = indx1*myp.N2+indx2;
                    int j = i2*myp.N2 + j2;
                    if(j != id) { 
                    complex<double> w2 = -dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)) * conj((out[findx] / (double)(totN)));
                    //if(abs(out[findx]) > 0.01*abs(out[0])) {
                    coeffs.push_back(Trip(id, j, w2));
                    
                    //}
                    
                    }

                }
            } 

        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    

    

    for (int ifn = 1; ifn < myp.number_of_fields; ifn++) // cross fields
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {

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
                int id = i1 * myp.N1 + j1;

                complex<double> w = - 1.0 - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2));

                coeffs.push_back(Trip(ifn * totN + id, ifn * totN+id, w));
            }
        }
    }



    for (int ifn = 0; ifn < myp.number_of_fields; ifn++) // cross fields
    {
        for (int ifn2 = ifn+1; ifn2 < myp.number_of_fields; ifn2++) // cross fields
        {

            for (int i1 = 0; i1 < myp.N1; i1++)
            {
                for (int j1 = 0; j1 < myp.N2; j1++)
                {

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
                    int id = i1 * myp.N1 + j1;

                    complex<double> w = - dtalpha * diffusion * temp1 * epsilon_couplings(ifn,ifn2) *(SQR(k1) + SQR(k2));

                    coeffs.push_back(Trip(ifn *totN + id, ifn2 * totN + id, w));

                    coeffs.push_back(Trip(ifn2 * totN + id, ifn *totN + id, w));
                }
            }


        }
    }



    int N = myp.number_of_fields * myp.N1 * myp.N2;

    Eigen::SparseMatrix<complex<double> > A(N,N);

    A.setFromTriplets(coeffs.begin(),coeffs.end());
    

    return A;


}

void CHC::SetupFracScheme(int MM) {
    M=MM;
    int N = myp.N1 * myp.N2;
    int nof = myp.number_of_fields;
    for(int i = 1 ; i < M ; i++) {
    Eigen::VectorXcd dxv(nof*N);

    for(int j = 0 ; j < nof*N ; j++) {
        dxv[j]=0.0;
    }
    OLD.push_back(dxv);
    }

    Eigen::VectorXcd w0 = CalculateWeight(fields);

    // OLD.pop_back();
    OLD.push_back(w0);
    std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());
}


void CHC::UpdateWithNewton(bool det = true) {
    
    chems.Calculate_Results(fields); // calculate chemistry

 
    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);
    //  string ftschem = "ftchem";
    //  outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //      pausel();
    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();
    int N = myp.N1 * myp.N2;
    transformed1.Calculate_Results(fields); // calculate FT of fields

    //Eigen::VectorXcd orig(myp.number_of_fields*myp.N1*myp.N2);

    double dtalpha = pow(dt, alpha);
    Eigen::SparseMatrix<complex<double>> unchangedJ = CalculateJ(transformed1.calculated_reactions, true);
    Eigen::VectorXcd orig = SolveLinearProblem(unchangedJ, transformed1.calculated_reactions, transformed3.calculated_reactions);

    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {
                // int id = i1 + myp.N1 * j1;
                rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig(k * N + i1 * myp.N2 + j1);
            }
        }
    }

        int iter=0;
        double dx = 10.;

        bool cont = true;

        while(cont) {
            cout << iter << endl;
            //iter++;
            //Eigen::ConjugateGradient<Eigen::SparseMatrix<complex<double>>, Eigen::Lower | Eigen::Upper> cg;

            //cout << "CJ set up" << endl;
            bool usefull = false; // iter>30;
            Eigen::SparseMatrix<complex<double> > J = CalculateJ(rules.calculated_reactions,usefull);
            //cout << J.isApprox(J.transpose()) << endl;

            cout << "J calculated" << endl;

            //check for accuracy.
            //CalculateWeight(rules.calculated_reactions)

            
            Eigen::VectorXcd F = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);
            

            cout << "RHS computed" << endl;
            Eigen::VectorXcd dxv;

            // Eigen::SimplicialCholesky<Eigen::SparseMatrix<complex<double>> > chol(J);
            // Eigen::VectorXcd dxv = chol.solve(-F);
            // cg.compute(J);

            if(usefull) {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<complex<double> > > solver;
            solver.compute(J);
            dxv = solver.solve(-F);
            std::cout << "#iterations:     " << solver.iterations() << std::endl;
            std::cout << "estimated error: " << solver.error() << std::endl;
            }
            else{
            Eigen::SparseLU<Eigen::SparseMatrix<complex<double>>, Eigen::COLAMDOrdering<int>> solver;
            // fill A and b;
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(J);
            cout << "pattern analyzed" << endl;
            // Compute the numerical factorization
            solver.factorize(J);
            cout << "pattern factorized" << endl;
            // Use the factors to solve the linear system
            dxv = solver.solve(-F);
            }
            // // cout << "cg computed" << endl;
            // cout << "dxv solved" << endl;
            // Eigen::VectorXcd dxv = cg.solve(-F);
            // cout << "dxv solved" << endl;

            // if (maxdx > dx)
            // {
            //     //break;
            //     cout << "non-downward step" << endl;
            //     // for(int i = 0 ; i < 100 ; i++)
            //     // cout << rules.calculated_reactions[0][i] <<" " <<dxv[i] << endl;

            //     cout << maxdx << " " << dx << " " << it % (myp.N1*myp.N2)  << endl;

            //     if(dx < 0.01)
            //         break;


            // }
            // cout << maxdx << " " << it << endl;
            // cout << F[it] << endl;
            // cout << rules.calculated_reactions[0][it] << endl;
            // cout << dxv[it] << endl;
            // pausel();
            double lambda = 1.;
            // for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
            // {
            //     for (int i1 = 0; i1 < myp.N1; i1++)
            //     {
            //         for (int j1 = 0; j1 < myp.N2; j1++)
            //         {
            //             rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N2 + j1] + lambda*dxv[k * N + i1 * myp.N2 + j1];
            //         }
            //     }
            // }

            for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
            {
                for (int i1 = 0; i1 < myp.N1; i1++)
                {
                    for (int j1 = 0; j1 < myp.N2; j1++)
                    {
                        rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N2 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1];
                        // rules.calculated_reactions[k][i1 * myp.N2 + j1] += lambda * dxv[k * N + i1 * myp.N2 + j1];
                        // cout << rules.calculated_reactions[k][i1 * myp.N2 + j1] << " " << orig[k * N + i1 * myp.N2 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1] << endl;
                    }
                }
            }
            //pausel();




            Eigen::VectorXcd F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

            double g0 = 0.5 * F.squaredNorm();
            double g1 = 0.5 * F2.squaredNorm();
            /*             cout << endl;
                        cout << "max original"  << endl;
                        maxdx = F.cwiseAbs().maxCoeff(&it);
                        cout << maxdx << " " << it << endl;

                        cout << "max after" << endl;
                        maxdx = F2.cwiseAbs().maxCoeff(&it);
                        cout << maxdx  << " " << it << endl;
                        double g0 = 0.5 * F.squaredNorm();
                        double g1=  0.5 * F2.squaredNorm();

                        cout << g0 << " " << g1 << endl;
                        pausel(); */
            double lambda2 = 0.5;
            int iter2 = 0;
            bool hoo;
            bool failed = false;
            
            while(g1>g0) {
                cout << "g1 bigger!" << g1 << " " << g0 << endl;
                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] +(lambda2) * dxv[k * N + i1 * myp.N2 + j1];
                        }
                    }
                }
                F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

                g1 = 0.5 * F2.squaredNorm();


                iter2++;
                if(iter2>10) {
                    failed = true;
                    break;
                }
                lambda = lambda2;
                lambda2 /= 2.;
            }

            if(failed) {
                cout << "Simple Jacobian failed!" << endl;
                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1];
                        }
                    }
                }
                J = CalculateJ(rules.calculated_reactions, true);
                Eigen::BiCGSTAB<Eigen::SparseMatrix<complex<double>>> solver;
                solver.compute(J);
                dxv = solver.solve(-F);
                cout << "new descent direction solved" << endl;
                std::cout << "#iterations:     " << solver.iterations() << std::endl;
                std::cout << "estimated error: " << solver.error() << std::endl;
                lambda = 0.1;
                //lambda2 = 0.5;

                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1];
                        }
                    }
                }
                F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

                g1 = 0.5 * F2.squaredNorm();
                g0 = 0.5 * F.squaredNorm();

                // cout << "orig: " << orig[1] << endl;

                double lambda2 = 1/10.*lambda;
                while(g1 > g0) {
                    cout << "lambda2: " << lambda2 << endl;
                    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                    {
                        for (int i1 = 0; i1 < myp.N1; i1++)
                        {
                            for (int j1 = 0; j1 < myp.N2; j1++)
                            {
                                rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + lambda2 * dxv[k * N + i1 * myp.N2 + j1];
                            }
                        }
                    }
                        F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

                        g1 = 0.5 * F2.squaredNorm();

                        lambda = lambda2;
                        lambda2 /=10.;
                }

                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + lambda2 * dxv[k * N + i1 * myp.N2 + j1];
                        }
                    }
                }
                // cout << "after: " << orig[1] << endl;
                vector<cdi> allt;
                int alli = 0;
                for(double lambda3 = 0.0 ; lambda3 <= 100*lambda2 ; lambda3+=lambda2) {
                    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                    {
                        for (int i1 = 0; i1 < myp.N1; i1++)
                        {
                            for (int j1 = 0; j1 < myp.N2; j1++)
                            {
                                rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + lambda3 * dxv[k * N + i1 * myp.N2 + j1];
                            }
                        }
                    }
                    F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

                    g1 = 0.5 * F2.squaredNorm();

                    allt.push_back({alli, g1});
                    alli++;
                }

                std::sort(allt.begin(), allt.end(),
                          [](const auto &i, const auto &j)
                          { return i.Score < j.Score; });

                //cout << allt[0].index << " " << allt[0].Score << endl;
                

                lambda = allt[0].index * lambda2;

                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1];
                        }
                    }
                }

                cout << "lambda used: " << lambda << endl;
                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            orig[k * N + i1 * myp.N2 + j1] = rules.calculated_reactions[k][i1 * myp.N2 + j1];
                        }
                    }
                }

                // dx = maxdx;
                int it;
                double maxdx = (lambda * dxv).cwiseAbs().maxCoeff(&it);
                cout << maxdx << endl;
                if (maxdx < 0.05 * lambda)
                {
                    cont = false;
                }
                iter++;
                // iter2 = 0;
                // bool hoo;
                // while (g1 > g0)
                // {
                //     cout << "g1 bigger!" << g1 << " " << g0 << endl;

                //     for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                //     {
                //         for (int i1 = 0; i1 < myp.N1; i1++)
                //         {
                //             for (int j1 = 0; j1 < myp.N2; j1++)
                //             {
                //                 rules.calculated_reactions[k][i1 * myp.N2 + j1] += (lambda2 - lambda) * dxv[k * N + i1 * myp.N2 + j1];
                //             }
                //         }
                //     }
                //     F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

                //     g1 = 0.5 * F2.squaredNorm();

                //     iter2++;
                //     if (iter2 > 20)
                //     {
                //         failed = true;
                //         for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                //         {
                //             for (int i1 = 0; i1 < myp.N1; i1++)
                //             {
                //                 for (int j1 = 0; j1 < myp.N2; j1++)
                //                 {
                //                     rules.calculated_reactions[k][i1 * myp.N2 + j1] -= (lambda2 - lambda) * dxv[k * N + i1 * myp.N2 + j1];
                //                 }
                //             }
                //         }
                //         cout <<  "Failed to converge!" << endl;
                //         cout << "trying another initial guess!" << endl;
                //         // pausel();
                //         // break;
                //         this->UpdateWithNewton(false); //try another starting point
                //         //error("do this if it doesn't work");
                //     }
                //     lambda = lambda2;
                //     lambda2 /= 2.;
                // }
            }
            else{
                cout << "lambda used: " << lambda << endl;
                for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
                {
                    for (int i1 = 0; i1 < myp.N1; i1++)
                    {
                        for (int j1 = 0; j1 < myp.N2; j1++)
                        {
                            orig[k * N + i1 * myp.N2 + j1] = rules.calculated_reactions[k][i1 * myp.N2 + j1];
                        }
                    }
                }

                //dx = maxdx;
                int it;
                double maxdx = (F2).cwiseAbs().maxCoeff(&it);
                if(maxdx < 0.0001*lambda) {
                    cont = false;
                }
                iter++;
            }
        
        }
        cout << "non-linear system solved" << endl;
        //once it is complete, update the field
        Eigen::VectorXcd Fx = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);
        int ik;
        double xx =  Fx.cwiseAbs().maxCoeff(&ik);
        cout << "max diff: " << xx <<" at " << ik << endl;



        reverse_transform.Calculate_Results(rules.calculated_reactions);

        //GetMaximas(fields, myp);

        Eigen::VectorXcd w0 = CalculateWeight(reverse_transform.calculated_reactions);

        OLD.pop_back();
        OLD.push_back(w0);
        std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());

        // for(int j = 0 ; j < OLD.size() ; j++)
        // cout << OLD[j][1] << ",";
        // cout << endl;

        reverse_transform.GetMaximas();
        reverse_transform.GetMaximasIndex();
        reverse_transform.GetMinimas();
        reverse_transform.GetMinimasIndex();
        cout << endl;

        set_field(reverse_transform.calculated_reactions);

        //after the code has finished working, 
    
    
    
                    // the current field configuration is in fields
                    // we need a starting vector to solve, for which we will choose the
}

void CHC::UpdateWithNewtonCalcJ()
{
    typedef complex<double> cd;
    chems.Calculate_Results(fields); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);
    //  string ftschem = "ftchem";
    //  outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //      pausel();
    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();

    transformed1.Calculate_Results(fields); // calculate FT of fields

    double dtalpha = pow(dt, alpha);
    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {

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
                // int id = i1 + myp.N1 * j1;
                rules.calculated_reactions[k][i1 * myp.N2 + j1] = transformed1.calculated_reactions[k][i1 * myp.N2 + j1] * 1. / (1. + dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2)));
            }
        }
    }
    int iter = 0;
    double dx = 10.;
    int N = myp.N1 * myp.N2;
    while (dx > 0.01)
    {
        cout << iter << endl;
        iter++;
        // Eigen::ConjugateGradient<Eigen::SparseMatrix<complex<double>>, Eigen::Lower | Eigen::Upper> cg;

        // cout << "CJ set up" << endl;

        // Eigen::SparseMatrix<complex<double>> J = CalculateJ(rules.calculated_reactions);
        // cout << J.isApprox(J.transpose()) << endl;
        vector<Trip> coeffs;
        for (int fieldno = 0; fieldno < 1; fieldno++)
        {
            for (int jk = 0; jk < myp.N1 * myp.N2; jk++)
            {

                cd **fields1;
                cd **fields2;

                fields1 = new cd *[myp.number_of_fields];
                fields2 = new cd *[myp.number_of_fields];

                for (int i = 0; i < myp.number_of_fields; i++)
                {
                    fields1[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                    fields2[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                    for (int j = 0; j < myp.N1 * myp.N2; j++)
                    {
                        fields1[i][j] = rules.calculated_reactions[i][j];
                        fields2[i][j] = rules.calculated_reactions[i][j];
                    }
                }
                cd dc = 0.001;
                fields1[fieldno][jk] += dc;

                fields2[fieldno][jk] -= dc;

                Eigen::VectorXcd w0 = CalculateRHS(fields1, transformed1.calculated_reactions, transformed3.calculated_reactions);
                Eigen::VectorXcd w1 = CalculateRHS(fields2, transformed1.calculated_reactions, transformed3.calculated_reactions);

                for (int fieldno1 = 0; fieldno1 < 1; fieldno1++)
                {
                    for (int p = 0; p < myp.N1 * myp.N2; p++)
                    {
                        // cout << "(" << fieldno << "," << jk << "), (" << fieldno1 << "," << p << ")" << endl;
                        // cout << "J: " <<J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                        // cout << "JT: " << J.coeff(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk) << endl;
                        // cout << "WT: " << (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc) << endl;
                        // cout << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)) /J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                        // pausel();
                        cd w = (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc);
                        //myfile << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)).real() << ",";
                        if(abs(w)>1E-6)
                        // coeffs.push_back(Trip(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p, w));
                        coeffs.push_back(Trip(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk, w));

                        //J(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) = (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc);
                        }
                }

                //myfile << endl;

                for (int i = 0; i < myp.number_of_fields; i++)
                {
                    fftw_free(fields1[i]);
                    fftw_free(fields2[i]);
                }
                delete fields1;
                delete fields2;
            }
        }
        int totN = myp.N1 *myp.N2;

        for (int ifn = 1; ifn < myp.number_of_fields; ifn++) // cross fields
        {
            for (int i1 = 0; i1 < myp.N1; i1++)
            {
                for (int j1 = 0; j1 < myp.N2; j1++)
                {

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
                    int id = i1 * myp.N1 + j1;

                    complex<double> w = -1.0 - dtalpha * diffusion * temp1 * (SQR(k1) + SQR(k2));

                    coeffs.push_back(Trip(ifn * totN + id, ifn * totN + id, w));
                }
            }
        }

        for (int ifn = 0; ifn < myp.number_of_fields; ifn++) // cross fields
        {
            for (int ifn2 = ifn + 1; ifn2 < myp.number_of_fields; ifn2++) // cross fields
            {

                for (int i1 = 0; i1 < myp.N1; i1++)
                {
                    for (int j1 = 0; j1 < myp.N2; j1++)
                    {

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
                        int id = i1 * myp.N1 + j1;

                        complex<double> w = -dtalpha * diffusion * temp1 * epsilon_couplings(ifn, ifn2) * (SQR(k1) + SQR(k2));

                        coeffs.push_back(Trip(ifn * totN + id, ifn2 * totN + id, w));

                        coeffs.push_back(Trip(ifn2 * totN + id, ifn * totN + id, w));
                    }
                }
            }
        }
        Eigen::SparseMatrix<complex<double> > J(myp.number_of_fields * myp.N1 * myp.N2, myp.number_of_fields * myp.N1 * myp.N2);

        cout << "J calculated" << endl;
        J.setFromTriplets(coeffs.begin(), coeffs.end());
        // check for accuracy.
        // CalculateWeight(rules.calculated_reactions)
        

        Eigen::VectorXcd F = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

        cout << "RHS computed" << endl;

        Eigen::BiCGSTAB<Eigen::SparseMatrix<complex<double>>> solver;
        solver.compute(J);
        Eigen::VectorXcd dxv = solver.solve(-F);
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error() << std::endl;

        // cout << rules.calculated_reactions[0][0] << " " << transformed1.calculated_reactions[0][0] << " " << transformed3.calculated_reactions[0][0]  << endl;
        // Eigen::SimplicialCholesky<Eigen::SparseMatrix<complex<double>> > chol(J);
        // Eigen::VectorXcd dxv = chol.solve(-F);
        // cg.compute(J);
        // Eigen::SparseLU<Eigen::SparseMatrix<complex<double>>, Eigen::COLAMDOrdering<int>> solver;
        // // fill A and b;
        // // Compute the ordering permutation vector from the structural pattern of A
        // solver.analyzePattern(J);
        // cout << "pattern analyzed" << endl;
        // // Compute the numerical factorization
        // solver.factorize(J);
        // cout << "pattern factorized" << endl;
        // // Use the factors to solve the linear system
        // Eigen::VectorXcd dxv = solver.solve(-F);
        // // cout << "cg computed" << endl;
        // cout << "dxv solved" << endl;
        // Eigen::VectorXcd dxv = cg.solve(-F);
        // cout << "dxv solved" << endl;
        int it;
        double maxdx = dxv.cwiseAbs().maxCoeff(&it);

       /*  if (maxdx > dx)
        {
            // break;
            cout << "non-downward step" << endl;
            cout << maxdx << " " << dx << endl;

            ofstream myfile;
            myfile.open("calcJ.csv");
            for (int fieldno = 0; fieldno < 5; fieldno++)
            {
                for (int jk = 0; jk < myp.N1 * myp.N2; jk++)
                {

                    cd **fields1;
                    cd **fields2;

                    fields1 = new cd *[myp.number_of_fields];
                    fields2 = new cd *[myp.number_of_fields];

                    for (int i = 0; i < myp.number_of_fields; i++)
                    {
                        fields1[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                        fields2[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                        for (int j = 0; j < myp.N1 * myp.N2; j++)
                        {
                            fields1[i][j] = rules.calculated_reactions[i][j];
                            fields2[i][j] = rules.calculated_reactions[i][j];
                        }
                    }
                    cd dc = 0.001;
                    fields1[fieldno][jk] += dc;

                    fields2[fieldno][jk] -= dc;

                    Eigen::VectorXcd w0 = CalculateRHS(fields1, transformed1.calculated_reactions, transformed3.calculated_reactions);
                    Eigen::VectorXcd w1 = CalculateRHS(fields2, transformed1.calculated_reactions, transformed3.calculated_reactions);

                    for (int fieldno1 = 0; fieldno1 < 5; fieldno1++)
                    {
                        for (int p = 0; p < myp.N1 * myp.N2; p++)
                        {
                            // cout << "(" << fieldno << "," << jk << "), (" << fieldno1 << "," << p << ")" << endl;
                            // cout << "J: " <<J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                            // cout << "JT: " << J.coeff(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk) << endl;
                            // cout << "WT: " << (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc) << endl;
                            // cout << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)) /J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                            // pausel();
                            myfile << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)).real();
                            if (p != 5 * myp.N1 * myp.N2 - 1 )
                                myfile << ",";
                        }
                    }
                    myfile << endl;

                    for (int i = 0; i < myp.number_of_fields; i++)
                    {
                        fftw_free(fields1[i]);
                        fftw_free(fields2[i]);
                    }
                    delete fields1;
                    delete fields2;
                }
            }
            myfile.close();

            ofstream myfile2;
            myfile2.open("J.csv");
            for (int k = 0; k < 5 * myp.N1 * myp.N2; k++)
            {
                for (int k1 = 0; k1 < 5 * myp.N1 * myp.N2; k1++)
                {
                    myfile2 << J.coeff(k, k1).real();
                    if (k1 != 5*myp.N1 * myp.N2-1)
                        myfile2 << ",";
                }
                myfile2 << endl;
            }
            myfile2.close();

            pausel();
        } */
        // cout << maxdx << " " << it << endl;
        // cout << F[it] << endl;
        // cout << rules.calculated_reactions[0][it] << endl;
        // cout << dxv[it] << endl;
        // pausel();

        for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
        {
            for (int i1 = 0; i1 < myp.N1; i1++)
            {
                for (int j1 = 0; j1 < myp.N2; j1++)
                {
                    rules.calculated_reactions[k][i1 * myp.N2 + j1] += dxv[k * N + i1 * myp.N2 + j1];
                }
            }
        }

        dx = maxdx;

        cout << dx << endl;
    }

    // once it is complete, update the field
    reverse_transform.Calculate_Results(rules.calculated_reactions);

    // GetMaximas(fields, myp);

    Eigen::VectorXcd w0 = CalculateWeight(reverse_transform.calculated_reactions);

    OLD.pop_back();
    OLD.push_back(w0);
    std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());

    for (int j = 0; j < OLD.size(); j++)
        cout << OLD[j][1] << ",";
    cout << endl;

    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;

    set_field(reverse_transform.calculated_reactions);

    // after the code has finished working,

    // the current field configuration is in fields
    // we need a starting vector to solve, for which we will choose the
}


Eigen::SparseMatrix<double> CHC::PartialJCalculation(Rule_Wrapper<cd,cd,cd,cd> &r) {
    int N = myp.N1 * myp.N2;
    int totN = N * myp.number_of_fields;
    double dtalpha = pow(dt, alpha);
    vector<TripD> coeffs;
    coeffs.reserve(N);

    {
        for (int index = 0; index < N; index++)
        {
            int fieldno = floor(index / N);
            int jk = index % N;

            cd *fields1;
            cd *fields2;
            cd *fields3;
            cd *fields4;

            fields1 = new cd;
            fields2 = new cd;
            fields3 = new cd;
            fields4 = new cd;

            fields1= (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
            fields2= (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
            fields3 = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
            fields4 = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));

            for (int j = 0; j < myp.N1 * myp.N2; j++)
            {
                fields1[j] = r.calculated_reactions[0][j];
                fields2[j] = r.calculated_reactions[0][j];
            }
            

            double dc = 0.001;
            fields1[jk] += dc;

            fields2[jk] -= dc;
            // Eigen::VectorXcd w0 = CalculateRHS_real(fields1, transformed1.calculated_reactions, transformed3.calculated_reactions);
            // Eigen::VectorXcd w1 = CalculateRHS_real(fields2, transformed1.calculated_reactions, transformed3.calculated_reactions);

            Eigen::VectorXcd w0 = CalculateWeightPartial(fields1);

            Eigen::VectorXcd w1 = CalculateWeightPartial(fields2);


            for (int fieldno1 = 0; fieldno1 < 1; fieldno1++)
            {
                for (int p = 0; p < myp.N1 * myp.N2; p++)
                {
                    // cout << "(" << fieldno << "," << jk << "), (" << fieldno1 << "," << p << ")" << endl;
                    // cout << "J: " <<J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                    // cout << "JT: " << J.coeff(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk) << endl;
                    // cout << "WT: " << (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc) << endl;
                    // cout << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)) /J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                    // pausel();
                    cd w = (w0[ p] - w1[ p]) / (2. * dc);
                    // myfile << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)).real() << ",";
                    // if(fieldno == fieldno1 && jk == p) w -= 1.0;
                    fields3[p] = w;
                    // if (abs(w) > 1E-6)
                    //     // coeffs.push_back(Trip(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p, w));
                    //     coeffs_private.push_back(Trip(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk, w));

                    // J(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) = (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc);
                }
            }

            //reverse_transform.Calculate_Results(store.calculated_reactions);

            fftw_plan p2;

            p2 = fftw_plan_dft_2d(myp.N1, myp.N2, reinterpret_cast<fftw_complex *>(fields3), reinterpret_cast<fftw_complex *>(fields4), FFTW_BACKWARD, FFTW_ESTIMATE);

            fftw_execute(p2);
            double corr = 1. / (double)(myp.N1);
            int end = myp.get_total();
            for (int i = 0; i < end; i++)
            {
                fields4[i] *= corr;
            }

            fftw_destroy_plan(p2);

   
                for (int p = 0; p < myp.N1 * myp.N2; p++)
                {
                    if (abs(fields4[p]) > 1E-3)
                    {
                        double w = fields4[p].real();
                        coeffs.push_back(TripD(jk, p, w));
                    }
                }
            
            // myfile << endl;

            fftw_free(fields1);
            fftw_free(fields2);
            fftw_free(fields3);
            fftw_free(fields4);
        }
    }

    Eigen::SparseMatrix<double> J(myp.N1 * myp.N2, myp.N1 * myp.N2);

    J.setFromTriplets(coeffs.begin(), coeffs.end());

    return J;
}

Eigen::SparseMatrix<double> CHC::CalculateInitialJ()
{

    typedef complex<double> cd;
    chems.Calculate_Results(fields); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);
    //  string ftschem = "ftchem";
    //  outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //      pausel();
    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();

    transformed1.Calculate_Results(fields); // calculate FT of fields
    int N = myp.N1*myp.N2;

    double dtalpha = pow(dt, alpha);
    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {

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
                // int id = i1 + myp.N1 * j1;
                rules.calculated_reactions[k][i1 * myp.N2 + j1] = fields[k][i1 * myp.N2 + j1];
            }
        }
    }

    

//calculate J
    cout << "beginning J calculation" << endl;
    vector<TripD> coeffs;
    int totN = myp.number_of_fields * myp.N1 * myp.N2;
    
	{
	    vector<TripD> coeffs_private;
        coeffs_private.reserve(N);
        #pragma omp for nowait schedule(static)
        for(int index = 0  ; index < totN ; index++)
        {
                int fieldno = floor(index/N);
                int jk = index % N;

                
                cout << fieldno << " " << jk << endl;
                cd **fields1;
                cd **fields2;

                fields1 = new cd *[myp.number_of_fields];
                fields2 = new cd *[myp.number_of_fields];

                for (int i = 0; i < myp.number_of_fields; i++)
                {
                    fields1[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                    fields2[i] = (cd *)fftw_malloc(myp.N1 * myp.N2 * sizeof(cd));
                    for (int j = 0; j < myp.N1 * myp.N2; j++)
                    {
                        fields1[i][j] = rules.calculated_reactions[i][j];
                        fields2[i][j] = rules.calculated_reactions[i][j];
                    }
                }
                
                double dc = 0.001;
                fields1[fieldno][jk] += dc;

                fields2[fieldno][jk] -= dc;
                // Eigen::VectorXcd w0 = CalculateRHS_real(fields1, transformed1.calculated_reactions, transformed3.calculated_reactions);
                // Eigen::VectorXcd w1 = CalculateRHS_real(fields2, transformed1.calculated_reactions, transformed3.calculated_reactions);
                
                Eigen::VectorXcd w0 = CalculateWeight(fields1);
                
                Eigen::VectorXcd w1 = CalculateWeight(fields2);

                Field_Wrapper<complex<double>, complex<double> > store(myp);

                for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
                {
                    for (int p = 0; p < myp.N1 * myp.N2; p++)
                    {
                        // cout << "(" << fieldno << "," << jk << "), (" << fieldno1 << "," << p << ")" << endl;
                        // cout << "J: " <<J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                        // cout << "JT: " << J.coeff(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk) << endl;
                        // cout << "WT: " << (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc) << endl;
                        // cout << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)) /J.coeff(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) << endl;
                        // pausel();
                        cd w = (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc);
                        // myfile << ((w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc)).real() << ",";
                        //if(fieldno == fieldno1 && jk == p) w -= 1.0;
                        store.calculated_reactions[fieldno1][p] = w;
                        // if (abs(w) > 1E-6)
                        //     // coeffs.push_back(Trip(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p, w));
                        //     coeffs_private.push_back(Trip(myp.N1 * myp.N2 * fieldno1 + p, fieldno * myp.N1 * myp.N2 + jk, w));

                        // J(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p) = (w0[myp.N1 * myp.N2 * fieldno1 + p] - w1[myp.N1 * myp.N2 * fieldno1 + p]) / (2. * dc);
                    }
                }

                reverse_transform.Calculate_Results(store.calculated_reactions);

                for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
                {
                    for (int p = 0; p < myp.N1 * myp.N2; p++)
                    {
                        if (abs(reverse_transform.calculated_reactions[fieldno1][p]) > 1E-3) {
                         double w = reverse_transform.calculated_reactions[fieldno1][p].real();
                         if(fieldno1 == fieldno && p == jk) w-=1.0; //account for the identity matrix
                            coeffs.push_back(TripD(fieldno * myp.N1 * myp.N2 + jk, myp.N1 * myp.N2 * fieldno1 + p, w));
                        }
                    }
                }
                // myfile << endl;

                for (int i = 0; i < myp.number_of_fields; i++)
                {
                    fftw_free(fields1[i]);
                    fftw_free(fields2[i]);
                }

                delete[] fields1;
                delete[] fields2;

            
        }

    }

    Eigen::SparseMatrix<double> J(myp.number_of_fields * myp.N1 * myp.N2, myp.number_of_fields * myp.N1 * myp.N2);

    J.setFromTriplets(coeffs.begin(), coeffs.end());
    Eigen::SparseMatrix<double> modJ = this->PartialJCalculation(rules);

    Eigen::SparseMatrix<double> unchangedJ(J);

    for (int k = 0; k < modJ.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(modJ, k); it; ++it)
        {
            double val = it.value();
            int row = it.row();   // row index
            int col = it.col();   // col index (here it is equal to k)
            
            //double val2 = J.coeff(row,col);
            unchangedJ.coeffRef(row,col) -= val;

            
        }


    return unchangedJ;

    /*
    cout << J.coeff(0,0) << " " << modJ.coeff(0,0) << " " << unchangedJ.coeff(0,0) << " " << endl;
    pausel();

    Eigen::VectorXcd F = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

    
    
    Field_Wrapper<complex<double>, complex<double>> store(myp);
    for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    {
        for (int p = 0; p < myp.N1 * myp.N2; p++)
        {
            store.calculated_reactions[fieldno1][p] = F[fieldno1*N+p];
        }
    }

    reverse_transform.Calculate_Results(store.calculated_reactions);



    Eigen::VectorXd F2(totN);
    for (int fieldno1 = 0; fieldno1 < myp.number_of_fields; fieldno1++)
    {
        for (int p = 0; p < myp.N1 * myp.N2; p++)
        {

            F2[fieldno1 * N + p] = - rules.calculated_reactions[fieldno1][p].real() + reverse_transform.calculated_reactions[fieldno1][p].real();
        }
    }
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(J);
    Eigen::VectorXcd dxv = solver.solve(-F2);
    */

    // cout << "RHS computed" << endl;
    // Eigen::VectorXcd dxv;

    // Eigen::BiCGSTAB<Eigen::SparseMatrix<complex<double>>> solver;
    // solver.compute(J);
    // dxv = solver.solve(-F);
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error() << std::endl;
    // int N = myp.N1 * myp.N2;

    // while(true) {
    //     double lambda = 1.;
    //     double x1norm = 0.0;
    //     for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    //     {
    //         for (int i1 = 0; i1 < myp.N1; i1++)
    //         {
    //             for (int j1 = 0; j1 < myp.N2; j1++)
    //             {
    //                 x1norm += SQR(abs(fields[k][i1 * myp.N2 + j1]));
    //                 rules.calculated_reactions[k][i1 * myp.N2 + j1] += lambda * dxv[k * N + i1 * myp.N2 + j1];
    //             }
    //         }
    //     }

    //     Eigen::VectorXcd F2 = CalculateRHS(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

    //     Eigen::VectorXcd dF = (F2 - F)/x1norm;

    //     Eigen::VectorXcd d2 = (J*dxv)/x1norm;

    //     dF -= d2;

    //     Eigen::SparseMatrix<complex<double>> Ju;

    // }
}



#endif /* CAHNHILLIARDCOMBO_CPP */

/*
void CHC::UpdateWithNewtonGivenJ(Eigen::SparseMatrix<double> &unchangedJ)
{
    int N = myp.N1 * myp.N2;
    int totN = myp.number_of_fields * N;


    chems.Calculate_Results(fields); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    transformed3.Calculate_Results(chems.calculated_reactions);
    //  string ftschem = "ftchem";
    //  outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //      pausel();
    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();

    transformed1.Calculate_Results(fields); // calculate FT of fields

    Eigen::VectorXd F(totN);

    double dtalpha = pow(dt, alpha);
    Eigen::VectorXd orig = SolveLinearProblem(unchangedJ,transformed1.calculated_reactions,transformed3.calculated_reactions);


    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                // int id = i1 + myp.N1 * j1;
                rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig(k * N + i1 * myp.N2 + j1);
            }
        }
    }
    F = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);
    int ik;
    double xx = F.cwiseAbs().maxCoeff(&ik);
    cout << "max diff: " << xx << " at " << ik << endl;
    

  

    int iter = 0;
    double dx = 10.;
    bool cont = true;

    while(cont) {
    Eigen::SparseMatrix<double> JforThis(unchangedJ);
    Eigen::SparseMatrix<double> modJ = this->PartialJCalculation(rules);

    for (int k = 0; k < modJ.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(modJ, k); it; ++it)
        {
            double val = it.value();
            int row = it.row(); // row index
            int col = it.col(); // col index (here it is equal to k)

            // double val2 = J.coeff(row,col);
            JforThis.coeffRef(row, col) += val;
        }

    Eigen::VectorXd F = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);


    int ik;
    double xx = F.cwiseAbs().maxCoeff(&ik);
    cout << "max diff: " << xx << " at " << ik << endl;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(JforThis);
    Eigen::VectorXd dxv = solver.solve(-F);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;
    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    // // fill A and b;
    // // Compute the ordering permutation vector from the structural pattern of A
    // solver.analyzePattern(JforThis);
    // cout << "pattern analyzed" << endl;
    // // Compute the numerical factorization
    // solver.factorize(JforThis);
    // cout << "pattern factorized" << endl;
    // // Use the factors to solve the linear system
    // Eigen::VectorXd dxv = solver.solve(-F2);

    double lambda = 1.0;

    lambda *= (1/(1+floor(iter/100)));

    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {
                rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N2 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1];
                //rules.calculated_reactions[k][i1 * myp.N2 + j1] += lambda * dxv[k * N + i1 * myp.N2 + j1];
                // cout << rules.calculated_reactions[k][i1 * myp.N2 + j1] << " " << orig[k * N + i1 * myp.N2 + j1] + lambda * dxv[k * N + i1 * myp.N2 + j1] << endl;
            }
        }
    }
    Eigen::VectorXd F2 = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

    double g0 = 0.5 * F.squaredNorm();
    double g1 = 0.5 * F2.squaredNorm();
    cout << "gs: " << g0 << " " << g1 << endl;
    // double lambda2 = 0.5;
    // while (g1 > g0)
    // {
    //     cout << "g1 bigger!" << g1 << " " << g0 << endl;
    //     for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    //     {
    //         for (int i1 = 0; i1 < myp.N1; i1++)
    //         {
    //             for (int j1 = 0; j1 < myp.N2; j1++)
    //             {
    //                 rules.calculated_reactions[k][i1 * myp.N2 + j1] = orig[k * N + i1 * myp.N1 + j1] + (lambda2)*dxv[k * N + i1 * myp.N2 + j1];
    //             }
    //         }
    //     }
    //     F2 = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);

    //     g1 = 0.5 * F2.squaredNorm();

    //     lambda = lambda2;
    //     lambda2 /= 2.;
    // }

    cout << "lambda used: " << lambda << endl;
    for (int k = 0; k < myp.number_of_fields; k++) // set up our initial guess at the solution
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j1 = 0; j1 < myp.N2; j1++)
            {
                orig[k * N + i1 * myp.N2 + j1] = rules.calculated_reactions[k][i1 * myp.N2 + j1].real();
            }
        }
    }

    int it;
    double maxdx = (F2).cwiseAbs().maxCoeff(&it);
    cout << "F2 maximum: " << maxdx << endl;
    if (maxdx < 0.005 * lambda)
    {
        cont = false;
    }
    iter++;

    }
    cout << "non-linear system solved" << endl;
    // once it is complete, update the field
    // Eigen::VectorXcd Fx = CalculateRHS_real(rules.calculated_reactions, transformed1.calculated_reactions, transformed3.calculated_reactions);
    // int ik;
    // double xx = F2.cwiseAbs().maxCoeff(&ik);
    // cout << "max diff: " << xx << " at " << ik << endl;


    Eigen::VectorXcd w0 = CalculateWeight(rules.calculated_reactions);

    OLD.pop_back();
    OLD.push_back(w0);
    std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());

    // for(int j = 0 ; j < OLD.size() ; j++)
    // cout << OLD[j][1] << ",";
    // cout << endl;

    rules.GetMaximas();
    rules.GetMaximasIndex();
    rules.GetMinimas();
    rules.GetMinimasIndex();
    rules.GetTotal();
    cout << endl;

    set_field(rules.calculated_reactions);
}
*/