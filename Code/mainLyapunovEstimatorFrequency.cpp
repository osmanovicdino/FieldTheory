#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <random>
//#include <mutex>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <chrono>
#include <type_traits>
//#include <stdafx>
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "DataStructures/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "CahnHilliard/cahnhilliard.h"

#include "fftw3.h"

using namespace std;

int main(int argc, char **argv)
{
    srand(time(NULL));
    string importstring;
    string importstring1;
    string importstring2;
    double val;
    int iter;
    int freq1,freq2;
    if (argc == 8)
    {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();
        stringstream ss2,ss3;
        ss2 << argv[2];
        ss3 << argv[3];
        val = atof(argv[4]);
        iter = atof(argv[5]);
        freq1 = atof(argv[6]);
        freq2 = atof(argv[7]);
        importstring1 = ss2.str();
        importstring2 = ss3.str();
    }
    else
    {
        error("no");
    }

    double T;
    bool err1;
    matrix<double> mat1 = importcsv(importstring, T, err1);


    int n = mat1(0, 0);

    vector1<bool> ps(n);

    for(int i = 0 ; i < n ; i++) {
        ps[i] = (bool)mat1(1,i);
    }


    vector1<double> simparams(4);
    for (int i = 0; i < 5; i++)
        simparams[i] = mat1(2, i);


    matrix<double> phasesepsparams(n,5);

    for(int j = 0 ; j < n ; j++ )
        for (int i = 0; i < 5; i++)
            phasesepsparams(j,i) = mat1(3+j, i);


    vector1<double> epsi((n) * (n - 1) / 2);

    for (int i = 0; i < (n) * (n - 1) / 2; i++)
    {
        epsi[i] = mat1(3+n, i);
    }

    int k = 0;
    matrix<double> epsa(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            epsa(i, j) = epsi[k];
            epsa(j, i) = epsi[k];
            k++;
        }
    }

    vector1<double> init(n);

    for (int i = 0; i < n; i++)
        init[i] = mat1(4+n, i);


    // if (argc == 5)
    // {
    //     c0 = atof(argv[1]);
    //     c1 = atof(argv[2]);
    //     eps = atof(argv[3]);
    //     dens = atof(argv[4]);
    // }
    // else {
    //     c0 = 0.2;
    //     c1 = 0.8;
    //     eps = 0.2;
    //     dens = 0.5;
    // }

    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc, myc> FWCC;
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    CH_builder p;
    int nof = n;
    p.number_of_fields = nof;
    p.N1 = 256;
    p.N2 = 256;

    CHC<double> a(p);
    CHC<double> b(p);

    for(int i = 0 ; i < nof ; i++) {
        for(int j = i+1  ; j < nof ; j++) {
            //cout << i << " " << j << endl;
             a.set_interaction(epsa(i,j), i, j);
             b.set_interaction(epsa(i, j), i, j);
        }
    }

    //ps[1] = true;


    a.set_phase_separators(ps);
    b.set_phase_separators(ps);

    // double nuc = 2.77778;
    for(int i = 0 ; i < n ; i++) {
        
        if(ps[i]) {
        a.set_diffusion(phasesepsparams(i,0), i);
        a.set_epsilon(phasesepsparams(i,1), i);
        a.set_c0_c1(phasesepsparams(i,2), phasesepsparams(i,3), i, phasesepsparams(i,4));

        b.set_diffusion(phasesepsparams(i, 0), i);
        b.set_epsilon(phasesepsparams(i, 1), i);
        b.set_c0_c1(phasesepsparams(i, 2), phasesepsparams(i, 3), i, phasesepsparams(i, 4));
        }
        else{
        a.set_diffusion(phasesepsparams(i, 0), i);
        b.set_diffusion(phasesepsparams(i, 0), i);
        }
    }




    // nuc = 1*0.563958;
    // a.set_diffusion(phasesepsparams[0], 1);
    // a.set_epsilon(0.664,1);
    // a.set_c0_c1(0.395,3.,1,nuc);

    // a.set_diffusion(phasesepsparams[0], 1);
    double L = simparams[0];
    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);
    a.set_alpha(simparams[1]);
    a.set_dt(simparams[2]);

    b.set_temp1(temp1);
    b.set_alpha(simparams[1]);
    b.set_dt(simparams[2]);

    double rate_multiplier = simparams[3];

    

    Field_Wrapper<double,double> my_chemsitry(p);
    NoWeight<double,double> nw;
    for (int j = 5+n; j < mat1.getnrows(); j++)
    {
        int no_chem = mat1(j, 0);


        if (no_chem == 0)
        {
            my_chemsitry.add_method(nw, j - (5+n));
        }
        else
        {
            MultipleReactions<double> c6(no_chem);

            double tot = 0.0;

            for (int i = 0; i < no_chem; i++)
            {
                // cout << i << endl;
                vector1<int> jpow(nof);
                for (int k = i * (nof + 1) + 2; k < i * (nof + 1) + 2 + nof; k++)
                {

                    jpow[k - (i * (nof + 1) + 2)] = (int)mat1(j, k);
                }
                // cout << jpow << endl;

                GenericChemistry<double> c6_0(rate_multiplier*mat1(j, i * (nof + 1) + 1), jpow);
                c6.add_chemical_reaction(c6_0, i);
                double tot1 = 1.0;
                for (int k = 0; k < nof; k++)
                {
                    tot1 *= Power(init[k], jpow[k]);
                    
                }

                // cout << mat1(j, i * (nof + 1) + 1) <<","<< tot1 << "," <<mat1(j, i * (nof + 1) + 1) * tot1 << endl;

                tot += mat1(j, i * (nof + 1) + 1) * tot1;
                
                
            }
            // cout << tot << endl;
            

            my_chemsitry.add_method(c6, j - (5+n));
        }
    }
    // pausel();
    a.set_chems(my_chemsitry);
    b.set_chems(my_chemsitry);

    a.setup_matrices();
    b.setup_matrices();

    // vector<matrix<double>> v;

    // for (int j = 0; j < nof; j++)
    // {
    //     matrix<double> field1(p.N1, p.N2);//1.E-012double)RAND_MAX) - 1.);

    //         for (int j = 0; j < p.N2; j++)
    //         {
    //            //double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);

    //             v[0](i, j) = mati1(i,j);
    //             v[1](i, j) = mati2(i,j);
    //         }
    //     }
    
    /* 

    {
        int centroidx = 128;
        int centroidy = 128;
        double x1 = init[0];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                if(SQR(i-centroidx)+SQR(j-centroidy) < 2500) {
                v[0](i, j) = x1 + gt * x1 * r1;
                }
                else{
                v[0](i, j) = 0.01;
                }
            }
        }
    }
    for (int lk = 1; lk < 5; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1 / 2; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = x1 + gt * x1 * r1;
            }
        }

        for (int i = p.N1 / 2; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = 0.01;
            }
        }
    }

    {
        int centroidx = 1024-128;
        int centroidy = 128;
        double x1 = init[5];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                if (SQR(i - centroidx) + SQR(j - centroidy) < 2500)
                {
                v[5](i, j) = x1 + gt * x1 * r1;
                }
                else
                {
                v[5](i, j) = 0.01;
                }
            }
        }
    }
    for (int lk = 6; lk < nof; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1 / 2; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = 0.01;
            }
        }

        for (int i = p.N1 / 2; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = x1 + gt * x1 * r1;
            }
        }
    } */

    // for (int lk = 0; lk < nof; lk++)
    // {  
    //     if(lk!=0 && lk!=5) {
    //     double x1 = init[lk];
    //     for (int i = 0; i < p.N1; i++)
    //     {
    //         for (int j = 0; j < p.N2; j++)
    //         {
    //             double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
    //             v[lk](i, j) = x1 + gt * x1 * r1;
    //         }
    //     }
    //     }
        
    // }

    // string importstring1 = "/home/dino/External/Waves/chemistry19/field0res_i=3078_real.csv";
    matrix<double> mat3 = importcsv(importstring1, T, err1);

    // string importstring2="/home/dino/External/Waves/chemistry19/field1res_i=3078_real.csv";
    matrix<double> mat4 = importcsv(importstring2, T, err1);

    // for (int lk = 0; lk < nof; lk++)
    // {
    //     a.set_field(v[lk], lk);
    // }
    a.set_field(mat3, 0);
    a.set_field(mat4, 1);

    matrix<double> perturb1(256,256);
    matrix<double> perturb2(256,256);

    for (int i = 0; i < p.N1; i++)
    {
        

        for (int j = 0; j < p.N2; j++)
        {
            // double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            double r2 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            perturb1(i, j) = r1;
            perturb2(i, j) = r2;
        }
    }

    b.set_field(perturb1,0);
    b.set_field(perturb2,1);
    b.transformed1.Calculate_Results(b.fields);


    int cut_off1 = freq1;
    int cut_off2 = freq2;

    for (int fn = 0; fn < nof; fn++)
    {
        for (int i = 0; i < b.myp.N1; i++)
        {
            for (int j = 0; j < b.myp.N2; j++)
            {
                int rel;
                double tempor = b.getk(i, j, b.myp.N1, b.myp.N2, rel);
                if (tempor < cut_off1 || tempor > cut_off2)
                {
                    // cout << tempor << endl;
                    b.transformed1.calculated_reactions[fn][i * b.myp.N2 + j] = 0.; // cut off the high frequency modes
                }
            }
        }
    }

    
    b.reverse_transform.Calculate_Results(b.transformed1.calculated_reactions);

    for(int i = 0  ; i < b.myp.N1 ; i++) {
        for (int j = 0; j < b.myp.N2; j++) {
            perturb1(i, j) = b.reverse_transform.calculated_reactions[0][i * 256 + j];
            perturb2(i, j) = b.reverse_transform.calculated_reactions[1][i * 256 + j];
        }
    }

    cout << "fourier done" << endl;
    outfunc(perturb1, "temp1");
    outfunc(perturb2, "temp2");
    pausel();

    double rescale = 0.0;

    for (int i = 0; i < p.N1; i++)
    {

        for (int j = 0; j < p.N2; j++)
        {
            // double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            rescale += perturb1(i, j) * perturb1(i, j);
            rescale += perturb2(i, j) * perturb2(i, j);
        }
    }

    cout << rescale << endl;

    perturb1 /= sqrt(rescale);
    perturb2 /= sqrt(rescale);

    

    double eps = val;
    
    cout << eps << endl;

    matrix<double> perturbn1(perturb1);

    matrix<double> perturbn2(perturb2);

    perturbn1 *= eps;
    perturbn2 *= eps;
    perturbn1 += mat3;
    perturbn2 += mat4;

    b.set_field(perturbn1, 0);
    b.set_field(perturbn2, 1);

    double dis2 = 0;
    double dis3 = 0;
    for (int i1 = 0; i1 < p.N1; i1++)
    {
        for (int j1 = 0; j1 < p.N2; j1++)
        {
            dis2 += SQR(a.fields[0][i1 * p.N1 + j1] - b.fields[0][i1 * p.N1 + j1]);
            dis2 += SQR(a.fields[1][i1 * p.N1 + j1] - b.fields[1][i1 * p.N1 + j1]);
            dis3 += SQR(perturb1(i1, j1));
            dis3 += SQR(perturb2(i1, j1));
        }
    }
    cout << fixed;
    cout << std::scientific << dis2 << endl;
    cout << "original distance " << sqrt(dis2) << endl;
    cout << dis3 << endl;

    // pausel();


    a.calculate_initial_weight(SQR(512));
    b.calculate_initial_weight(SQR(512));

  
    cout << "calc" << endl;
    

    cout << "all fields set" << endl;

    // auto start = std::chrono::high_resolution_clock::now();
    int runtime = 100000;
    int every = 1;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);
    vector1<bool> ps2(nof,true);
    // auto start = std::chrono::high_resolution_clock::now();
    ofstream myfileLya;
    string s1 = importstring;
    string su = s1.substr(0, s1.size() - 4);
    stringstream ssv;
    ssv << eps;
    stringstream ssn;
    ssn << iter;
    stringstream ssx1,ssx2;
    ssx1 << cut_off1;
    ssx2 << cut_off2;
    string sh = string("cutoff1=")+ssx1.str()+string("_cutoff2=")+ssx2.str();
    myfileLya.open( (su+string("Lyapunov_diff=") +ssv.str()+string("_try=") +ssn.str()+sh+string(".csv")).c_str() );

    double diss = 0.;
    for (int i1 = 0; i1 < p.N1; i1++)
    {
        for (int j1 = 0; j1 < p.N2; j1++)
        {
            diss += SQR(a.fields[0][i1 * p.N1 + j1] - b.fields[0][i1 * p.N1 + j1]);
            diss += SQR(a.fields[1][i1 * p.N1 + j1] - b.fields[1][i1 * p.N1 + j1]);
        }
    }

    // double lambda1 = (1. / (every * simparams[2])) * log(sqrt(dis1) / eps);
    // myfileLya << lambda1 << endl;
    myfileLya << sqrt(diss) << endl;

    for (int i = 0; i < runtime; i++)
    {
        if(i % 10000000 ==0  && i > 0000) {

            string s1 = importstring;
            // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / 1000;
            string s2 = "_i=" + ss.str();
            string su = s1.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_some_results(string("a") + su + s2, ps2);
            b.print_some_results(string("b") + su + s2, ps2);
        }

        if (i % every == 0 && i > 0000 )
        {
            // stringstream strep1;
            // stringstream strep2;
            // stringstream strep3;
            // stringstream strep4;
            double dis1 = 0.;
            for(int i1 = 0 ; i1 < p.N1 ; i1++) {
                for(int j1 = 0 ; j1 < p.N2 ; j1++) {
                    dis1 += SQR(a.fields[0][i1*p.N1+j1]-b.fields[0][i1*p.N1+j1]);
                    dis1 += SQR(a.fields[1][i1 * p.N1 + j1] - b.fields[1][i1 * p.N1 + j1]);
                }
            }
            cout << sqrt(dis1) << endl;
            cout << eps << endl;
            double lambda1 = (1. / (every * simparams[2])) * log(sqrt(dis1) / eps);
            // myfileLya << lambda1 << endl;
            myfileLya << sqrt(dis1) << endl;
            
            //  matrix<double> orig1(p.N1,p.N2);
            // matrix<double> orig2(p.N1, p.N2);
            // for (int i1 = 0; i1 < p.N1; i1++)
            // {
            //     for (int j1 = 0; j1 < p.N2; j1++)
            //     {
            //         orig1(i1,j1) = a.fields[0][i1 * p.N1 + j1];
            //         orig2(i1, j1) = a.fields[1][i1 * p.N1 + j1];
            //     }
            // }
            
            // b.set_field(orig1 + eps * perturb1, 0);
            // b.set_field(orig2 + eps * perturb2, 1);

            // diss = 0.;
            // for (int i1 = 0; i1 < p.N1; i1++)
            // {
            //     for (int j1 = 0; j1 < p.N2; j1++)
            //     {
            //         diss += SQR(a.fields[0][i1 * p.N1 + j1] - b.fields[0][i1 * p.N1 + j1]);
            //         diss += SQR(a.fields[1][i1 * p.N1 + j1] - b.fields[1][i1 * p.N1 + j1]);
            //     }
            // }

            // // double lambda1 = (1. / (every * simparams[2])) * log(sqrt(dis1) / eps);
            // // myfileLya << lambda1 << endl;
            // cout << sqrt(diss) << endl;
            // pausel();

            // strep1 << dens;
            // strep4 << c0;
            // strep2 << c1;
            // strep3 <<  surf;
            
            // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
            // stringstream ss;
            // ss << setw(number_of_digits) << setfill('0') << i / every;
            // string s2 = "_i=" + ss.str();
            
            // cout << su + s2 << endl;
            // string bs = "b";
            // a.print_some_results(su + s2,ps2);
            // b.print_some_results(bs + su + s2, ps2);
            
        }
        cout << i << endl;
        cout << "begin" << endl;
        a.Update();
        b.Update();
        bool chck = true;
        a.check_field(chck);
        if (!chck)
            break;
    }
}