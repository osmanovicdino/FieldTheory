#ifndef CAHNHILLIARDCOMBO_CPP
#define CAHNHILLIARDCOMBO_CPP

CHC::CHC(const CH_builder &p) : CH(p), epsilon_couplings(matrix<double>(p.number_of_fields, p.number_of_fields)),
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
                - 0.0*2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j]
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
#endif /* CAHNHILLIARDCOMBO_CPP */
