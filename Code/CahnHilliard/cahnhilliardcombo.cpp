#ifndef CAHNHILLIARDCOMBO_CPP
#define CAHNHILLIARDCOMBO_CPP

template<class T>
CHC<T>::CHC(const CH_builder &p) : CH<T>(p), phase_separators(vector1<bool>(p.number_of_fields)),
                                diffusion(vector1<double>(p.number_of_fields)),
                                epsilon(vector1<double>(p.number_of_fields)),
                                c0(vector1<double>(p.number_of_fields)),
                                c1(vector1<double>(p.number_of_fields)),
                                cons1(vector1<double>(p.number_of_fields)),
                                cons2(vector1<double>(p.number_of_fields)),
                                cons3(vector1<double>(p.number_of_fields)),
                                cons4(vector1<double>(p.number_of_fields)),
                                cons3s(vector1<double>(p.number_of_fields)),
                                epsilon_couplings(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                epsilon_couplingsSQR(matrix<double>(p.number_of_fields, p.number_of_fields)),
                                inverses(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                baremat(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                oldfieldFT(Field_Wrapper<T, T>(p)),
                                oldfieldNLW(Field_Wrapper<T, T>(p)),
                                InitWeight(Field_Wrapper<T, T>(p))
{
    for (int i = 0; i < this->myp.number_of_fields; i++)
    {
        IdentityWeight<T> fw;
        oldfieldFT.add_method(fw, i);
        oldfieldNLW.add_method(fw, i);
    }
    //M=1;

}

template<class T>
void CHC<T>::setup_matrices() {

    matrix<double> identitymatrix(this->myp.number_of_fields, this->myp.number_of_fields);
    for(int i  = 0  ; i < this->myp.number_of_fields ; i++) {
        identitymatrix(i,i) = 1.;
    }

    if constexpr (std::is_same_v<T, complex<double> >) {
    int iter = 0;
    for(int k1 = 0  ; k1 <= this->myp.N1/2 ; k1 ++) {
        for (int k2 = k1; k2 <= this->myp.N2 / 2; k2++)
        {
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
    else if constexpr (std::is_same_v<T, double>)
    {
    baremat.resize(this->myp.N1*this->myp.N2); 
    inverses.resize(this->myp.N1*this->myp.N2); 
    
    int iter = 0;
    
    for (int k1 = 0; k1 < this->myp.N1 ; k1++)
    {
        for (int k2 = 0; k2 < this->myp.N2 ; k2++)
        {

            matrix<double> temp = create_D_mat_split(k1, k2);

            baremat[iter] = temp;

            matrix<double> to_invert = identitymatrix - (1 - alpha) * temp - alpha * dt * temp;
            to_invert.inverse();

            inverses[iter] = to_invert;

            iter++;
        }
        }
    }
    else{

    }


}


template <class T>
void CHC<T>::set_interaction(double val, int i, int j)
{
    epsilon_couplings(i,j) = val;
    epsilon_couplings(j,i) = val;
}

template <class T>
void CHC<T>::calculate_non_linear_weight(T **input)
{

    // takes as an input a real space input of fields, produces as an output the fourier transformed weight function into transformed2

    int totp =  this->myp.get_total();



    for(int i  = 0 ; i < this->myp.number_of_fields; i++) {
        if(phase_separators[i]) {
        for (int j = 0; j < totp; j++)
        {
            this->weigs.calculated_reactions[i][j] = cons1[i] * CUB(input[i][j]) + cons2[i] * SQR(input[i][j]) + cons4[i];
        }
        }
        else{
        for(int j = 0 ; j < totp ; j++ )
            this->weigs.calculated_reactions[i][j] = 0.0;
        }
    }

  //  outfunc(weigs.calculated_reactions[0],"test", myp);

    this->transformed2.Calculate_Results(this->weigs.calculated_reactions);

    // outfunc(transformed2.calculated_reactions[0], "testft", myp);
}

// template <class T>
// void CHC<T>::calculate_non_linear_weightSQR(T **input)
// { 
    

//     // takes as an input a real space input of fields, produces as an output the fourier transformed weight function into transformed2

//     int totp = myp.get_total();


//     double a = 5.;
//     double x0 = 0.5;
//     for (int i = 0; i < myp.number_of_fields; i++)
//     {
        
//             for (int j = 0; j < totp; j++) {
//                 weigs.calculated_reactions[i][j] = 0.0;

//                 for (int k = 0; k < myp.number_of_fields; k++)
//                 {
//                     //weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * input[i][j] * SQR(input[k][j]);
//                     if(k!=i) {
//                     if(epsilon_couplingsSQR(i,k) < 0.) {
//                         weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * (a * SQR(1. / cosh(a * (input[i][j] - x0))) + a * SQR(1. / cosh(a * (input[i][j] - x0))) * tanh(a * (input[k][j] - x0)));
//                     }
//                     // else {
//                     //     weigs.calculated_reactions[i][j] += epsilon_couplingsSQR(i, k) * input[i][j] * SQR(input[k][j]);
//                     // }
//                     }
//                     // if (input[i][j] < 0.)
//                     // {
                        
//                     //     cout << input[i][j] << " " << input[k][j]<< endl;
//                     //     cout << epsilon_couplingsSQR(i, k) << endl;
//                     //     cout << epsilon_couplingsSQR(i, k) * (a * exp(a * input[i][j]) * (1. - 1. / (1. + exp(a * input[k][j])))) / SQR(1. + exp(a * input[i][j])) << endl;
//                     //     pausel();
//                     // }
//                 }
                
//             }
//     }

//     // add a phi^4 term to stop over crowding for the other chemicals.
//     // for (int i = 1; i < myp.number_of_fields; i++)
//     // {
//     //     for (int j = 0; j < totp; j++)
//     //         weigs.calculated_reactions[i][j] += 0.25 * CUB(input[i][j]);

//     // }
//     for (int i = 0; i < myp.number_of_fields; i++)
//     {
//         if(phase_separators[i])
//         for (int j = 0; j < totp; j++)
//         {
//             weigs.calculated_reactions[i][j] += cons1[i] * CUB(input[i][j]) + cons2[i] * SQR(input[i][j]) + cons4[i];
//         }
//     }

//     // for (int i = 0; i < myp.number_of_fields; i++)
//     // {
//     //     for (int j = 0; j < totp; j++)
//     //     if (input[i][j] < 0.)
//     //     {
//     //         for(int k  = 0; k <  myp.number_of_fields ; k++) {
//     //         cout << epsilon_couplingsSQR(i, k) << " " << input[i][j] << " " << input[k][j]<< endl;
//     //         // cout << epsilon_couplingsSQR(i, k) << endl;
//     //         // cout << epsilon_couplingsSQR(i, k) * (a * exp(a * input[i][j]) * (1. - 1. / (1. + exp(a * input[k][j])))) / SQR(1. + exp(a * input[i][j])) << endl;
//     //         }
//     //         cout << "weight: " <<  weigs.calculated_reactions[i][j];

//     //         int effj = j % myp.N1 ;
//     //         int effi = j / myp.N1;

//     //         int effin1 = effi + 1 > myp.N1 ? 0 : effi + 1;
//     //         int effin2 = effi - 1 < 0 ? myp.N1 : effi - 1;

//     //         int effjn1 = effj + 1 > myp.N1 ? 0 : effj + 1;
//     //         int effjn2 = effj - 1 < 0 ? myp.N1 : effj - 1;

//     //         int totj1 = effin1 * myp.N1 + effj;

//     //         int totj2 = effin2 * myp.N1 + effj;

//     //         int totj3 = effi * myp.N1 + effjn1;

//     //         int totj4 = effi * myp.N1 + effjn2;

//     //         cout << "n1 weight: " << weigs.calculated_reactions[i][totj1];
//     //         cout << "n2 weight: " << weigs.calculated_reactions[i][totj2];
//     //         cout << "n3 weight: " << weigs.calculated_reactions[i][totj3];
//     //         cout << "n4 weight: " << weigs.calculated_reactions[i][totj4];
//     //         pausel();
//     //     }
//     // }

//     // for(int j = 0  ; j < totp ; j++) {
//     //     if(abs(input[1][j])>2.) {
//     //         cout << input[1][j] << endl;
//     //         for (int k = 0; k < myp.number_of_fields; k++)
//     //         {
                
//     //                 cout << "k: " << epsilon_couplingsSQR(1, k) * input[1][j] * SQR(input[k][j]) << ",";
//     //                 cout << input[k][j] << ",";
//     //                 cout << endl;
//     //         }
//     //         cout << endl;
//     //         cout << weigs.calculated_reactions[1][j] << endl;
//     //         cout << 0.2 * CUB(input[1][j]) << endl;
//     //         pausel();
//     //     }
//     // }
//     //  outfunc(weigs.calculated_reactions[0],"test", myp);

//     transformed2.Calculate_Results(weigs.calculated_reactions);

//     // outfunc(transformed2.calculated_reactions[0], "testft", myp);
// }


// template<class T>
// double CHC<T>::getk(int i, int j, int N1, int N2, int &rel) {

//     if constexpr (std::is_same_v<T, complex<double>>)
//     {
//         int k1, k2;
//         if (i <= N1 / 2)
//         {
//         k1 = i;
//         }
//         else
//         {
//         k1 = (i - N1);
//         }
//         if (j <= N2 / 2)
//         {
//         k2 = j;
//         }
//         else
//         {
//         k2 = (j - N2);
//         }
//         double tempor = SQR(k1) + SQR(k2);

//         int k1a = abs(k1);
//         int k2a = abs(k2);
//         if (k2a < k1a)
//         {
//         k1a = k2a;
//         k2a = abs(k1);
//         }
//         rel = k2a - (k1a * (1 + k1a - 2 * (1 + N1 / 2))) / 2.;


//         return tempor;
//     }
//     else if constexpr(std::is_same_v<T,double>) {
//         int k1 = i;
//         int k2 = j;

//         double tempor = SQR(k1) + SQR(k2);
//         rel=k1*N2+k2;
//         return tempor;
//     }
//     else{
//         error("something weird");

//     }
// }

template<>
double CHC<complex<double> >::getk(int i, int j, int N1, int N2, int &rel)
{

    // if constexpr (std::is_same_v<T, complex<double>>)
    // {
        int k1, k2;
        if (i <= N1 / 2)
        {
        k1 = i;
        }
        else
        {
        k1 = (i - N1);
        }
        if (j <= N2 / 2)
        {
        k2 = j;
        }
        else
        {
        k2 = (j - N2);
        }
        double tempor = SQR(k1) + SQR(k2);

        int k1a = abs(k1);
        int k2a = abs(k2);
        if (k2a < k1a)
        {
        k1a = k2a;
        k2a = abs(k1);
        }
        rel = k2a - (k1a * (1 + k1a - 2 * (1 + N1 / 2))) / 2.;

        return tempor;
    }

    template <>
    double CHC<double>::getk(int i, int j, int N1, int N2, int &rel)
    {

        // if constexpr (std::is_same_v<T, complex<double>>)
        // {
            int k1 = i;
            int k2 = j;

            double tempor = SQR(k1) + SQR(k2);
            rel = k1 * N2 + k2;
            return tempor;
    }
    // else if constexpr (std::is_same_v<T, double>)
    // {
    //     int k1 = i;
    //     int k2 = j;

    //     double tempor = SQR(k1) + SQR(k2);
    //     rel = k1 * N2 + k2;
    //     return tempor;
    // }
    // else
    // {
    //     error("something weird");
    // }


template<class T>
void CHC<T>::calculate_initial_weight(int cut_offf) {
    this->transformed1.Calculate_Results(this->fields); //calculate FT of fields
    int totp = this->myp.get_total();
    int nof = this->myp.number_of_fields;

    //cut off form
    int cut_off = cut_offf;
    for (int fn = 0; fn < nof; fn++)
    {
        for (int i = 0; i < this->myp.N1; i++)
        {
            for (int j = 0; j < this->myp.N2; j++)
            {
            int rel;
            double tempor = getk(i, j, this->myp.N1, this->myp.N2, rel);
            if (tempor > cut_off)
            {
               // cout << tempor << endl;
                this->transformed1.calculated_reactions[fn][i * this->myp.N2 + j] = 0.; // cut off the high frequency modes
                    
                    }
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
    this->reverse_transform.Calculate_Results(this->transformed1.calculated_reactions);

    this->set_field(this->reverse_transform.calculated_reactions);


    // outfunc(fields[0], "fi1", myp);
    // outfunc(fields[1], "fi2", myp);

    // pausel();

    // reverse_transform.rescale();


//    set_field(reverse_transform.calculated_reactions);

    this->transformed1.Calculate_Results(this->fields); // calculate FT of fields

    oldfieldFT.Calculate_Results(this->transformed1.calculated_reactions);

    calculate_non_linear_weight(this->fields);
    oldfieldNLW.Calculate_Results(this->transformed2.calculated_reactions);

    for (int i = 0; i < this->myp.N1; i++)
    {
        for (int j = 0; j < this->myp.N2; j++)
        {

            // double k1, k2;
            // if (i <= this->myp.N1 / 2)
            // {
            //     k1 = i;
            // }
            // else
            // {
            //     k1 = (i - this->myp.N1);
            // }
            // if (j <= this->myp.N2 / 2)
            // {
            //     k2 = j;
            // }
            // else
            // {
            //     k2 = (j - this->myp.N2);
            // }

            // //take absolute values of k1 and k2
            // int k1a = abs(k1);
            // int k2a = abs(k2);
            // if(k2a < k1a) {
            //     k1a = k2a;
            //     k2a = abs(k1);
            // }
            // int rel = k2a - (k1a * (1 + k1a - 2 * (1 + this->myp.N1 / 2))) / 2.;

            // double tempor = SQR(k1) + SQR(k2);
            int rel;
            double tempor = getk(i, j, this->myp.N1, this->myp.N2, rel);

            vector1<T> v(nof);
            for(int k = 0  ; k < nof ; k++) {
                v[k] = this->transformed1.calculated_reactions[k][i * this->myp.N2 + j];
            }

            v = baremat[rel]*v;


            for(int k = 0 ; k < nof ; k++) {
                InitWeight.calculated_reactions[k][i * this->myp.N2 + j] = (v[k] - diffusion[k] * tempor * temp1 * this->transformed2.calculated_reactions[k][i * this->myp.N2 + j]);
                //oldweightFT.calculated_reactions[k][i * myp.N2 + j] = (1- alpha) * (-v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
            }
            //upd1[i * myp.N2 + j] = -2*alpha*dt*diffusion*tempor*temp1*transformed2.calculated_reactions[]
            

            
            }
    }

    this->reverse_transform.Calculate_Results(InitWeight.calculated_reactions);

    // outfunc(reverse_transform.calculated_reactions[0], "fi1", myp);
    // outfunc(reverse_transform.calculated_reactions[1], "fi2", myp);
    cout << "done" << endl;

}

// template <class T>
// void CHC<T>::calculate_initial_weightSQR()
// {

//     this->transformed1.Calculate_Results(this->fields); // calculate FT of fields

//     //cut off large wavenumbers

//     int totp = this->myp.get_total();
//     int nof = this->myp.number_of_fields;

//     for(int fn = 0 ; fn < nof ; fn++) {
//     for (int i = 200; i < myp.N1; i++)
//     {
//         for (int j = 200; j < myp.N2; j++)
//         {
//             transformed1.calculated_reactions[fn][i * myp.N2 + j]=0.;
//         }
//     }

//     }
//     reverse_transform.Calculate_Results(transformed1.calculated_reactions);
//     set_field(reverse_transform.calculated_reactions);

//     transformed1.Calculate_Results(fields); // calculate FT of fields

//     oldfieldFT.Calculate_Results(transformed1.calculated_reactions);

//     calculate_non_linear_weightSQR(fields);
//     oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);



//     for (int i = 0; i < myp.N1; i++)
//     {
//         for (int j = 0; j < myp.N2; j++)
//         {

//             double k1, k2;
//             if (i <= myp.N1 / 2)
//             {
//                 k1 = i;
//             }
//             else
//             {
//                 k1 = (i - myp.N1);
//             }
//             if (j <= myp.N2 / 2)
//             {
//                 k2 = j;
//             }
//             else
//             {
//                 k2 = (j - myp.N2);
//             }

//             // take absolute values of k1 and k2
//             int k1a = abs(k1);
//             int k2a = abs(k2);
//             if (k2a < k1a)
//             {
//                 k1a = k2a;
//                 k2a = abs(k1);
//             }
//             int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

//             double tempor = SQR(k1) + SQR(k2);

//             vector1<T> v(nof);
//             for (int k = 0; k < nof; k++)
//             {
//                 v[k] = transformed1.calculated_reactions[k][i * myp.N2 + j];
//             }

//             v = baremat[rel] * v;

//             for (int k = 0; k < nof; k++)
//             {
//                 InitWeight.calculated_reactions[k][i * myp.N2 + j] = (v[k] - diffusion[k] * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
//                 // oldweightFT.calculated_reactions[k][i * myp.N2 + j] = (1- alpha) * (-v[k] - diffusion * tempor * temp1 * transformed2.calculated_reactions[k][i * myp.N2 + j]);
//             }
//             // upd1[i * myp.N2 + j] = -2*alpha*dt*diffusion*tempor*temp1*transformed2.calculated_reactions[]
//         }
//     }
// }

template<class T>
void CHC<T>::Update() {

    // string fieldss = "fields";
    // outfunc(fields[0], fieldss, myp);

    this->chems.Calculate_Results(this->fields); // calculate chemistry

    this->transformed3.Calculate_Results(this->chems.calculated_reactions);

    // string ftschem = "ftchem";
    // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
    //     pausel();
    // GetMaximas(transformed3.calculated_reactions,schem,myp);
    // pausel();

    this->transformed1.Calculate_Results(this->fields); // calculate FT of fields
    calculate_non_linear_weight(this->fields);

    int totp = this->myp.get_total();
    int nof = this->myp.number_of_fields;

    for (int i = 0; i < this->myp.N1; i++)
    {
            for (int j = 0; j < this->myp.N2; j++)
            {

            // double k1, k2;
            // if (i <= this->myp.N1 / 2)
            // {
            //     k1 = i;
            // }
            // else
            // {
            //     k1 = (i - this->myp.N1);
            // }
            // if (j <= this->myp.N2 / 2)
            // {
            //     k2 = j;
            // }
            // else
            // {
            //     k2 = (j - this->myp.N2);
            // }

            // //take absolute values of k1 and k2
            // int k1a = abs(k1);
            // int k2a = abs(k2);
            // if (k2a < k1a)
            // {
            //     k1a = k2a;
            //     k2a = abs(k1);
            // }
            // int rel = k2a - (k1a * (1 + k1a - 2 * (1 + this->myp.N1 / 2))) / 2.;

            // double tempor = SQR(k1) + SQR(k2);
            int rel;
            double tempor = getk(i, j, this->myp.N1, this->myp.N2, rel);

            // vector1<complex<double> > v(nof);


            double fac =  tempor * temp1 ;
            vector1<T > v(nof);
            vector1<T > v2(nof);
            vector1<T > v3(nof);

            for (int k = 0; k < nof; k++)
            {
                v2[k] = this->transformed1.calculated_reactions[k][i * this->myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
            }

            v2 = baremat[rel] * v2;


            

            for (int k = 0; k < nof; k++)
            {
                v[k] = -((1-alpha) + (alpha) * dt) * diffusion[k]*fac * (this->transformed2.calculated_reactions[k][i * this->myp.N2 + j]) 
                + this->transformed1.calculated_reactions[k][i * this->myp.N2 + j]//oldfieldFT.calculated_reactions[k][i * myp.N2 + j] 
                - (1-alpha)*v2[k] 
                + ((1-alpha) /* -  0.5*(alpha) * dt */) *diffusion[k] * fac * oldfieldNLW.calculated_reactions[k][i * this->myp.N2 + j] 
                - 0.0*2 * alpha * dt * InitWeight.calculated_reactions[k][i * this->myp.N2 + j] //zeroed because of the initial conditions we set
                + dt * this->transformed3.calculated_reactions[k][i * this->myp.N2 + j];
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
                this->rules.calculated_reactions[k][i * this->myp.N2 + j] = v3[k];
            }
        }
    }

    oldfieldFT.Calculate_Results(this->transformed1.calculated_reactions);
    oldfieldNLW.Calculate_Results(this->transformed2.calculated_reactions);

    this->reverse_transform.Calculate_Results(this->rules.calculated_reactions);

    GetMaximas(this->fields, this->myp);
    this->reverse_transform.GetMaximas();
    this->reverse_transform.GetMaximasIndex();
    this->reverse_transform.GetMinimas();
    this->reverse_transform.GetMinimasIndex();
    cout << endl;

    this->set_field(this->reverse_transform.calculated_reactions);
}

// template<class Q, class T>
// void CHC<T>::UpdateNoise(Q &func, GenNoise<T> &mynoise, vector1<double> &str)
// {

//     // string fieldss = "fields";
//     // outfunc(fields[0], fieldss, myp);

//     chems.Calculate_Results(fields); // calculate chemistry

//     // string schem = "chem";
//     // outfunc(chems.calculated_reactions[0],schem,myp);

//     transformed3.Calculate_Results(chems.calculated_reactions);
//     // string ftschem = "ftchem";
//     // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
//     //     pausel();

//     mynoise.GenFields(func, str, 1.0 / 2.0);

//     transformed3.Add_Noise(mynoise);


//     // GetMaximas(transformed3.calculated_reactions,schem,myp);
//     // pausel();

//     transformed1.Calculate_Results(fields); // calculate FT of fields
//     calculate_non_linear_weight(fields);

//     int totp = myp.get_total();
//     int nof = myp.number_of_fields;

//     for (int i = 0; i < myp.N1; i++)
//     {
//         for (int j = 0; j < myp.N2; j++)
//         {

//             double k1, k2;
//             if (i <= myp.N1 / 2)
//             {
//                 k1 = i;
//             }
//             else
//             {
//                 k1 = (i - myp.N1);
//             }
//             if (j <= myp.N2 / 2)
//             {
//                 k2 = j;
//             }
//             else
//             {
//                 k2 = (j - myp.N2);
//             }

//             // take absolute values of k1 and k2
//             int k1a = abs(k1);
//             int k2a = abs(k2);
//             if (k2a < k1a)
//             {
//                 k1a = k2a;
//                 k2a = abs(k1);
//             }
//             int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

//             double tempor = SQR(k1) + SQR(k2);

//             // vector1<complex<double> > v(nof);

//             double fac = tempor * temp1;
//             vector1<T> v(nof);
//             vector1<T> v2(nof);
//             vector1<T> v3(nof);

//             for (int k = 0; k < nof; k++)
//             {
//                 v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
//             }

//             v2 = baremat[rel] * v2;

//             for (int k = 0; k < nof; k++)
//             {
//                 v[k] = -((1 - alpha) + (alpha)*dt) * diffusion[k] * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                       // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
//                        - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * diffusion[k] * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] // zeroed because of the initial conditions we set
//                        + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
//             }

//             v3 = inverses[rel] * v;
//             // if (rel == 131839)
//             // {
//             //     cout << inverses[rel] << endl;
//             //     cout << nof << endl;
//             //     for (int k = 0; k < nof; k++)
//             //     {

//             //         cout << "before: " << transformed1.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "1: " << -(2 * (1 - alpha) - (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) << endl;
//             //         cout << "2: " << +oldfieldFT.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "3: " << -(1 - alpha) * v2[k] << endl;
//             //         cout << "4: " << +(2 * (1 - alpha) + (alpha)*dt) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "5: " << +2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "6: " << +2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //        // cout << "7: " << 2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "after: " << v3[k] << endl;
//             //         cout << endl;
//             //         cout << endl;
//             //     }

//             //     cout << v << endl;
//             //     cout << v3 << endl;
//             //     pausel();
//             // }
//             for (int k = 0; k < nof; k++)
//             {
//                 rules.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//             }
//         }
//     }

//     oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
//     oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

//     reverse_transform.Calculate_Results(rules.calculated_reactions);

//     GetMaximas(fields, myp);
//     reverse_transform.GetMaximas();
//     reverse_transform.GetMaximasIndex();
//     reverse_transform.GetMinimas();
//     reverse_transform.GetMinimasIndex();
//     cout << endl;

//     set_field(reverse_transform.calculated_reactions);
// }

// template<class T>
// void CHC<T>::UpdateSQR()
// {
//     GetMaximas(fields, myp);
//     GetMinimas(fields, myp);
//     // string fieldss = "fields";
//     // outfunc(fields[0], fieldss, myp);

//     chems.Calculate_Results(fields); // calculate chemistry


//     // string schem = "chem";
//     // outfunc(chems.calculated_reactions[0],schem,myp);

//     transformed3.Calculate_Results(chems.calculated_reactions);

//     // string ftschem = "ftchem";
//     // outfunc(transformed3.calculated_reactions[0], ftschem, myp);
//     //     pausel();

//     // GetMaximas(transformed3.calculated_reactions,schem,myp);



//     transformed1.Calculate_Results(fields); // calculate FT of fields
//     calculate_non_linear_weightSQR(fields);

//     // cout << transformed2.calculated_reactions[0][512*1024+512] << endl;

//     // cout << oldfieldNLW.calculated_reactions[0][512 * 1024 + 512] << endl;

//     // cout << transformed2.calculated_reactions[1][512 * 1024 + 512] << endl;

//     // cout << oldfieldNLW.calculated_reactions[1][512 * 1024 + 512] << endl;

//     // cout << transformed2.calculated_reactions[2][512 * 1024 + 512] << endl;

//     // cout << oldfieldNLW.calculated_reactions[2][512 * 1024 + 512] << endl;

//     // cout << transformed2.calculated_reactions[3][512 * 1024 + 512] << endl;

//     // cout << oldfieldNLW.calculated_reactions[3][512 * 1024 + 512] << endl;

//     // pausel();

//     int totp = myp.get_total();
//     int nof = myp.number_of_fields;

//     for (int i = 0; i < myp.N1; i++)
//     {
//         for (int j = 0; j < myp.N2; j++)
//         {

//             double k1, k2;
//             if (i <= myp.N1 / 2)
//             {
//                 k1 = i;
//             }
//             else
//             {
//                 k1 = (i - myp.N1);
//             }
//             if (j <= myp.N2 / 2)
//             {
//                 k2 = j;
//             }
//             else
//             {
//                 k2 = (j - myp.N2);
//             }

//             // take absolute values of k1 and k2
//             int k1a = abs(k1);
//             int k2a = abs(k2);
//             if (k2a < k1a)
//             {
//                 k1a = k2a;
//                 k2a = abs(k1);
//             }
//             int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

//             double tempor = SQR(k1) + SQR(k2);

//             // vector1<complex<double> > v(nof);

//             double fac = tempor * temp1;
//             vector1<T> v(nof);
//             vector1<T> v2(nof);
//             vector1<T> v3(nof);

//             for (int k = 0; k < nof; k++)
//             {
//                 v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
//             }

//             v2 = baremat[rel] * v2;

//             for (int k = 0; k < nof; k++)
//             {
//                 v[k] = -((1 - alpha) + (alpha)*dt) * diffusion[k] * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                              // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
//                        - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * diffusion[k] * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] /* - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j]  */ // zeroed because of the initial conditions we set
//                        + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
//             }

//             v3 = inverses[rel] * v;
//             // if (rel == 131839)
//             // {
//             //     cout << inverses[rel] << endl;
//             //     cout << nof << endl;
//             //     for (int k = 0; k < nof; k++)
//             //     {

//             //         cout << "before: " << transformed1.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "1: " << -(2 * (1 - alpha) - (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) << endl;
//             //         cout << "2: " << +oldfieldFT.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "3: " << -(1 - alpha) * v2[k] << endl;
//             //         cout << "4: " << +(2 * (1 - alpha) + (alpha)*dt) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "5: " << +2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "6: " << +2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //        // cout << "7: " << 2 * dt * transformed3.calculated_reactions[k][i * myp.N2 + j] << endl;
//             //         cout << "after: " << v3[k] << endl;
//             //         cout << endl;
//             //         cout << endl;
//             //     }

//             //     cout << v << endl;
//             //     cout << v3 << endl;
//             //     pausel();
//             // }
//             for (int k = 0; k < nof; k++)
//             {
//                 rules.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//             }
//         }
//     }



//     reverse_transform.Calculate_Results(rules.calculated_reactions);


//     // bool checker = false;

//     // if(checker)
//     // for (int ifn = 0; ifn < myp.number_of_fields; ifn++)
//     // {
//     //     for(int i1 = 0 ; i1 < myp.N1 ; i1++)
//     //     for(int j1 = 0 ; j1 < myp.N2 ; j1++)
//     //     if(reverse_transform.calculated_reactions[ifn][i1*myp.N1+j1]<0.0)
//     //         {
//     //             cout << ifn << endl;
//     //             cout << fields[ifn][i1*myp.N1+j1] << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1*myp.N1+j1] << endl;

//     //             Field_Wrapper<complex<double>, complex<double>> store1(myp);
//     //             Field_Wrapper<complex<double>, complex<double>> store2(myp);
//     //             Field_Wrapper<complex<double>, complex<double>> store3(myp);
//     //             Field_Wrapper<complex<double>, complex<double>> store4(myp);
//     //             Field_Wrapper<complex<double>, complex<double>> store5(myp);
//     //             Field_Wrapper<complex<double>, complex<double>> store6(myp);
//     //             for (int i = 0; i < myp.N1; i++)
//     //             {
//     //                 for (int j = 0; j < myp.N2; j++)
//     //                 {

//     //                     double k1, k2;
//     //                     if (i <= myp.N1 / 2)
//     //                     {
//     //                         k1 = i;
//     //                     }
//     //                     else
//     //                     {
//     //                         k1 = (i - myp.N1);
//     //                     }
//     //                     if (j <= myp.N2 / 2)
//     //                     {
//     //                         k2 = j;
//     //                     }
//     //                     else
//     //                     {
//     //                         k2 = (j - myp.N2);
//     //                     }

//     //                     // take absolute values of k1 and k2
//     //                     int k1a = abs(k1);
//     //                     int k2a = abs(k2);
//     //                     if (k2a < k1a)
//     //                     {
//     //                         k1a = k2a;
//     //                         k2a = abs(k1);
//     //                     }
//     //                     int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

//     //                     double tempor = SQR(k1) + SQR(k2);
//     //                     vector1<double> fac =  (tempor * temp1) * diffusion;
//     //                     vector1<complex<double>> v(nof);
//     //                     vector1<complex<double>> v2(nof);
//     //                     vector1<complex<double>> v3(nof);

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
//     //                     }

//     //                     v2 = baremat[rel] * v2;

//     //                     vector1<complex<double>> vtemp1(nof);
//     //                     vector1<complex<double>> vtemp2(nof);
//     //                     vector1<complex<double>> vtemp3(nof);
//     //                     vector1<complex<double>> vtemp4(nof);
//     //                     vector1<complex<double>> vtemp5(nof);
//     //                     vector1<complex<double>> vtemp6(nof);

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         vtemp1[k] = -((1 - alpha) + (alpha)*dt) * fac[k] * (transformed2.calculated_reactions[k][i * myp.N2 + j]);
//     //                         vtemp2[k] = + transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
//     //                         vtemp3[k] = - (1 - alpha) * v2[k];
//     //                         vtemp4[k] = +((1 - alpha)) * fac[k] * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j];
//     //                         vtemp5[k] = + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
//     //                         vtemp6[k] = -((1 - alpha) + (alpha)*dt) * fac[k] * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                              // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
//     //                                     - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * fac[k] * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] /* - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j]  */ // zeroed because of the initial conditions we set
//     //                                     + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];
//     //                     }
//     //                     v3 = inverses[rel] * vtemp1;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store1.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }
//     //                     v3 = inverses[rel] * vtemp2;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store2.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }
//     //                     v3 = inverses[rel] * vtemp3;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store3.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }
//     //                     v3 = inverses[rel] * vtemp4;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store4.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }
//     //                     v3 = inverses[rel] * vtemp5;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store5.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }

//     //                     v3 = inverses[rel] * vtemp6;

//     //                     for (int k = 0; k < nof; k++)
//     //                     {
//     //                         store6.calculated_reactions[k][i * myp.N2 + j] = v3[k];
//     //                     }
//     //                 }
//     //             }

//     //             double k1, k2;
//     //             if (i1 <= myp.N1 / 2)
//     //             {
//     //                 k1 = i1;
//     //             }
//     //             else
//     //             {
//     //                 k1 = (i1 - myp.N1);
//     //             }
//     //             if (j1 <= myp.N2 / 2)
//     //             {
//     //                 k2 = j1;
//     //             }
//     //             else
//     //             {
//     //                 k2 = (j1 - myp.N2);
//     //             }

//     //             // take absolute values of k1 and k2
//     //             int k1a = abs(k1);
//     //             int k2a = abs(k2);
//     //             if (k2a < k1a)
//     //             {
//     //                 k1a = k2a;
//     //                 k2a = abs(k1);
//     //             }
//     //             int rel = k2a - (k1a * (1 + k1a - 2 * (1 + myp.N1 / 2))) / 2.;

//     //             double a = 10.;

//     //             for(int k  = 0 ; k < nof ; k++) {
//     //             if (epsilon_couplingsSQR(ifn, k) < 0.)
//     //             {
//     //                 int j = i1*myp.N1+j1;
//     //                 cout << "negative weight: ";
//     //                 cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
//     //                 cout << endl;
//     //             }
//     //             else
//     //             {
//     //                 int j = i1 * myp.N1 + j1;
//     //                 cout << "positive weight: ";
//     //                 cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
//     //                 cout << endl;
//     //             }
//     //             }

//     //             for (int k = 0; k < nof; k++)
//     //             {
//     //                 if (epsilon_couplingsSQR(ifn, k) < 0.)
//     //                 {
//     //                     int j = i1 * myp.N1 + j1 -1 ;
//     //                     cout << "negative weight n1: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
//     //                     cout << endl;
//     //                 }
//     //                 else
//     //                 {
//     //                     int j = i1 * myp.N1 + j1 -1 ;
//     //                     cout << "positive weight n1: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
//     //                     cout << endl;
//     //                 }
//     //             }
//     //             for (int k = 0; k < nof; k++)
//     //             {
//     //                 if (epsilon_couplingsSQR(ifn, k) < 0.)
//     //                 {
//     //                     int j = i1 * myp.N1 + j1 + 1;
//     //                     cout << "negative weight n2: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
//     //                     cout << endl;
//     //                 }
//     //                 else
//     //                 {
//     //                     int j = i1 * myp.N1 + j1 + 1;
//     //                     cout << "positive weight n2: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
//     //                     cout << endl;
//     //                 }
//     //             }
//     //             for (int k = 0; k < nof; k++)
//     //             {
//     //                 if (epsilon_couplingsSQR(ifn, k) < 0.)
//     //                 {
//     //                     int j = (i1 + 1) * myp.N1 + j1;
//     //                     cout << "negative weight n3: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
//     //                     cout << endl;
//     //                 }
//     //                 else
//     //                 {
//     //                     int j = (i1 + 1) * myp.N1 + j1 + 1;
//     //                     cout << "positive weight n3: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
//     //                     cout << endl;
//     //                 }
//     //             }

//     //             for (int k = 0; k < nof; k++)
//     //             {
//     //                 if (epsilon_couplingsSQR(ifn, k) < 0.)
//     //                 {
//     //                     int j = (i1 - 1) * myp.N1 + j1 + 1;
//     //                     cout << "negative weight n4: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * (fields[ifn][j] * SQR(fields[k][j]) * (a * fields[ifn][j] * SQR(1. / cosh(a * fields[ifn][j])) * (1. + tanh(a * fields[k][j])) + 2. * (tanh(a * fields[k][j]) + tanh(a * fields[ifn][j]) * (1. + tanh(a * fields[k][j])))));
//     //                     cout << endl;
//     //                 }
//     //                 else
//     //                 {
//     //                     int j = (i1 - 1) * myp.N1 + j1 + 1;
//     //                     cout << "positive weight n4: ";
//     //                     cout << epsilon_couplingsSQR(ifn, k) * fields[ifn][j] * SQR(fields[k][j]);
//     //                     cout << endl;
//     //                 }
//     //             }

//     //             cout << "inv" << endl;
//     //             cout << inverses[rel] << endl;



//     //             reverse_transform.Calculate_Results(transformed1.calculated_reactions);
//     //             cout << "i1: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store1.calculated_reactions);
//     //             cout << "t1: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store2.calculated_reactions);
//     //             cout << "t2: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store3.calculated_reactions);
//     //             cout << "t3: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store4.calculated_reactions);
//     //             cout << "t4: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store5.calculated_reactions);
//     //             cout << "t5: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             reverse_transform.Calculate_Results(store6.calculated_reactions);
//     //             cout << "t6: " << endl;
//     //             cout << reverse_transform.calculated_reactions[ifn][i1 * myp.N1 + j1] << endl;

//     //             pausel();
//     //         }
//     //     }

//         oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
//         oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

//         GetMaximas(fields, myp);
//         reverse_transform.GetMaximas();
//         reverse_transform.GetMaximasIndex();
//         reverse_transform.GetMinimas();
//         reverse_transform.GetMinimasIndex();
//         cout << endl;

//         set_field(reverse_transform.calculated_reactions);
// }



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