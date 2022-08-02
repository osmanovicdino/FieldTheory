#ifndef CAHNHILLIARDDOUBLE_CPP
#define CAHNHILLIARDDOUBLE_CPP

CHD::CHD(const CH_builder &p) : CH(p), inverses(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4))),
                                baremat(vector<matrix<double>>((1 + p.N1 / 2) * (1 + p.N2 / 4)))
{

}

void CHD::setup_matrices()
{

    //matrix<double> identitymatrix(myp.number_of_fields, myp.number_of_fields);

    int iter = 0;
    for (int k1 = 0; k1 <= myp.N1 / 2; k1++)
    {
        for (int k2 = k1; k2 <= myp.N2 / 2; k2++)
        {
            
            matrix<double> temp = create_D_mat_split(k1, k2);


            baremat[iter] = temp;

            matrix<double> to_invert = temp;
            to_invert.inverse();

            inverses[iter] = to_invert;

            iter++;
        }
    }
}

void CHD::set_interaction_and_diffusion(double x12, double x13, double x23, double D1, double D2) {

    diffusion1 = D1;
    diffusion2 = D2;

    chi_12 = x12;
    chi_13 = x13;
    chi_23 = x23;

    /*
    string scr = "/home/dino/Documents/Chemistry/TwoFluid/scripts/fm.wls";
    stringstream p1,p2,p3;

    p1 << chi_12;
    p2 << chi_13;
    p3 << chi_23;

    string param1 = p1.str();
    string param2 = p2.str();
    string param3 = p3.str();

    string gap = " ";

    string total_command = scr + gap + param1 + gap + param2 + gap + param3;


    cout << total_command << endl;
    
    system(total_command.c_str());

    double T;
    bool err1;
    matrix<double> res = importcsv("/home/dino/Documents/Chemistry/TwoFluid/scripts/output.csv",T,err1);

    x0 = res(0,0);
    y0 = res(1,0);
    */

    x0 = 1./3.;
    y0 = 1./3.;

    cout << x0 << " " << y0 << endl;

    ml1 = diffusion1*(2 * (1 / (2. * x0) - 1 / (2. * (-1 + x0 + y0)) - chi_13));
    ml2 = diffusion1*(-(1 / (-1 + x0 + y0)) + chi_12 - chi_13 - chi_23);
    ml3 = (diffusion2/diffusion1)*ml2;
    ml4 = diffusion2*(2 * (1 / (2. * y0) - 1 / (2. * (-1 + x0 + y0)) - chi_23));

    sl2 = diffusion1*0.5*(chi_12 - chi_13 - chi_23);
    sl1 = -diffusion1*chi_13;
    sl4 = -diffusion2*chi_23;
    sl3 = diffusion2 * 0.5 * (chi_12 - chi_13 - chi_23);

    // cout << ml1 << "," << ml2 << endl;
    // cout << ml3 << "," << ml4 << endl;
    // cout << endl;
    // cout << sl1 << "," << sl2 << endl;
    // cout << sl3 << "," << sl4 << endl;
}

// void CHD::set_diffusion(double D1, double D2) {
//     diffusion1 = D1;
//     diffusion2 = D2;
// }

void CHD::calculate_non_linear_weight(complex<double> **input)
{
    int totp = myp.get_total();

    for (int j = 0; j < totp; j++)
    {
        weigs.calculated_reactions[0][j] = 3. * (-0.16666666666666666 * 1 / SQR(x0) + 1 / (6. * SQR(-1 + x0 + y0))) * SQR(input[0][j]) +
                                           4. * (1 / (12. * CUB(x0)) - 1 / (12. * CUB(-1 + x0 + y0))) * CUB(input[0][j]) +
                                           (input[0][j] / SQR(-1 + x0 + y0) - SQR(input[0][j]) / CUB(-1 + x0 + y0)) * input[1][j] +
                                           (1 / (2. * SQR(-1 + x0 + y0)) - input[0][j] / CUB(-1 + x0 + y0)) * SQR(input[1][j]) - CUB(input[1][j]) / (3. * CUB(-1 + x0 + y0));
        weigs.calculated_reactions[1][j] = SQR(input[0][j]) / (2. * SQR(-1 + x0 + y0)) - CUB(input[0][j]) / (3. * CUB(-1 + x0 + y0)) +
                                           2. * (input[0][j] / (2. * SQR(-1 + x0 + y0)) - SQR(input[0][j]) / (2. * CUB(-1 + x0 + y0))) * input[1][j] +
                                           3. * (-0.16666666666666666 * 1 / SQR(y0) + 1 / (6. * SQR(-1 + x0 + y0)) - input[0][j] / (3. * CUB(-1 + x0 + y0))) * SQR(input[1][j]) +
                                           4. * (1 / (12. * CUB(y0)) - 1 / (12. * CUB(-1 + x0 + y0))) * CUB(input[1][j]);
    }

    transformed2.Calculate_Results(weigs.calculated_reactions);
}

void CHD::Update()
{

    // string fieldss = "fields";
    // outfunc(fields[0], fieldss, myp);

    //chems.Calculate_Results(fields); // calculate chemistry

    // string schem = "chem";
    // outfunc(chems.calculated_reactions[0],schem,myp);

    //transformed3.Calculate_Results(chems.calculated_reactions);
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

            double fac1 = diffusion1 * tempor * temp1;
            double fac2 = diffusion2 * tempor * temp1;
            vector1<complex<double>> v(nof);
            vector1<complex<double>> v2(nof);
            vector1<complex<double>> v3(nof);

            for (int k = 0; k < nof; k++)
            {
                v2[k] = transformed1.calculated_reactions[k][i * myp.N2 + j]; // oldfieldFT.calculated_reactions[k][i * myp.N2 + j];
            }

            //v2 = baremat[rel] * v2;

            

            v[0] = v2[0] - dt * fac1 * transformed2.calculated_reactions[0][i * myp.N2 + j];
            v[1] = v2[1] - dt * fac2 * transformed2.calculated_reactions[1][i * myp.N2 + j];
            // v[k] = -((1 - alpha) + (alpha)*dt) * fac * (transformed2.calculated_reactions[k][i * myp.N2 + j]) + transformed1.calculated_reactions[k][i * myp.N2 + j]                                                       // oldfieldFT.calculated_reactions[k][i * myp.N2 + j]
            //        - (1 - alpha) * v2[k] + ((1 - alpha) /* -  0.5*(alpha) * dt */) * fac * oldfieldNLW.calculated_reactions[k][i * myp.N2 + j] - 0.0 * 2 * alpha * dt * InitWeight.calculated_reactions[k][i * myp.N2 + j] // zeroed because of the initial conditions we set
            //        + dt * transformed3.calculated_reactions[k][i * myp.N2 + j];

            

            v3 = inverses[rel] * v;
            if (rel == 131839)
            {
                cout << "bare mat: " << endl;
                cout << baremat[rel] << endl;
                cout << "inverse mat: " << endl;
                cout << inverses[rel] << endl;
                
                for (int k = 0; k < nof; k++)
                {

                    cout << transformed2.calculated_reactions[k][i * myp.N2 + j] << endl;
                }
                cout << "old: " << endl;
                cout << v2 << endl;
                cout << "new: " << endl;
                cout << v << endl;
                pausel();
            }
            for (int k = 0; k < nof; k++)
            {
                rules.calculated_reactions[k][i * myp.N2 + j] = v3[k];
            }
        }
    }

    // string gilename1 = "gello1";
    // string gilename2 = "gello2";

    // string filename1 = "hello1";
    // string filename2 = "hello2";

    // outfunc(transformed1.calculated_reactions[0], gilename1, myp);

    // outfunc(transformed1.calculated_reactions[1], gilename2, myp);
    
    // outfunc(rules.calculated_reactions[0], filename1, myp);

    // outfunc(rules.calculated_reactions[1], filename2, myp);
    // pausel();
    // oldfieldFT.Calculate_Results(transformed1.calculated_reactions);
    // oldfieldNLW.Calculate_Results(transformed2.calculated_reactions);

    reverse_transform.Calculate_Results(rules.calculated_reactions);

    GetMaximas(fields, myp);
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;

    set_field(reverse_transform.calculated_reactions);
}

#endif