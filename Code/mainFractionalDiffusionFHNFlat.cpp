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
    // string importstring;

    // if (argc == 2)
    // {
    //     stringstream ss;
    //     ss << argv[1];
    //     importstring = ss.str();
    // }
    // else
    // {
    //     error("no");
    // }

    double T;
    bool err1;
   // matrix<double> mat1 = importcsv(importstring, T, err1);


    int n = 2;

    vector1<bool> ps(n);

    for(int i = 0 ; i < n ; i++) {
        ps[i] = 0;
    }


    // vector1<double> simparams(4);
    // for (int i = 0; i < 5; i++)
    //     simparams[i] = mat1(2, i);


    // matrix<double> phasesepsparams(n,5);

    // for(int j = 0 ; j < n ; j++ )
    //     for (int i = 0; i < 5; i++)
    //         phasesepsparams(j,i) = mat1(3+j, i);


    // vector1<double> epsi((n) * (n - 1) / 2);

    // for (int i = 0; i < (n) * (n - 1) / 2; i++)
    // {
    //     epsi[i] = mat1(3+n, i);
    // }

    int k = 0;
    matrix<double> epsa(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            epsa(i, j) = 0.;
            epsa(j, i) = 0.;
            k++;
        }
    }

    vector1<double> init(2);
    init[0]=0.;
    init[1]=0.;

    double L=770.0;
    double deltat=0.01;

    vector1<double> diffusion(n);
    diffusion[0] = 1.0;
    diffusion[1] = 3.01;

    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc, myc> FWCC;
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    CH_builder p;
    int nof = n;
    p.number_of_fields = nof;
    p.N1 = 512;
    p.N2 = 512;

    CHC<myc> a(p);

    for(int i = 0 ; i < nof ; i++) {
        for(int j = i+1  ; j < nof ; j++) {
            //cout << i << " " << j << endl;
             a.set_interaction(epsa(i,j), i, j);
        }
    }

    //ps[1] = true;


    a.set_phase_separators(ps);


    // double nuc = 2.77778;
    for(int i = 0 ; i < n ; i++) {
        
        // if(ps[i]) {
        // a.set_diffusion(phasesepsparams(i,0), i);
        // a.set_epsilon(phasesepsparams(i,1), i);
        // a.set_c0_c1(phasesepsparams(i,2), phasesepsparams(i,3), i, phasesepsparams(i,4));
        // }
        
        a.set_diffusion(diffusion[i], i);
        
    }




    // nuc = 1*0.563958;
    // a.set_diffusion(phasesepsparams[0], 1);
    // a.set_epsilon(0.664,1);
    // a.set_c0_c1(0.395,3.,1,nuc);

    // a.set_diffusion(phasesepsparams[0], 1);

    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);

    a.set_alpha(1.);
    a.set_dt(deltat);


    double rate_multiplier = 1.0;

    

    FWCC my_chemsitry(p);
    NoWeight<myc,myc> nw;
    double kf =  3.4;
    double alpha =0.186;
    double delta = 0.296;

    double kg =  3.23;
    double epsilon =  0.308;

    FitzHughNagumo<myc> fhn(kf, alpha,delta, 0);
    InvasionDecay<myc> u1(1., 1);

    MultipleReactions<myc> c6(2);
    c6.add_chemical_reaction(fhn, 0);
    c6.add_chemical_reaction(u1, 1);

    my_chemsitry.add_method(c6,0);

    InvasionDecay<myc> u2(-kg * epsilon, 0);
    InvasionDecay<myc> u3(epsilon, 1);

    MultipleReactions<myc> c7(2);
    c7.add_chemical_reaction(u2, 0);
    c7.add_chemical_reaction(u3, 1);

    my_chemsitry.add_method(c6, 0);

    my_chemsitry.add_method(c7, 1);
    // for (int j = 5+n; j < mat1.getnrows(); j++)
    // {
    //     int no_chem = mat1(j, 0);


    //     if (no_chem == 0)
    //     {
    //         my_chemsitry.add_method(nw, j - (5+n));
    //     }
    //     else
    //     {
    //         MultipleReactions<double> c6(no_chem);

    //         double tot = 0.0;

    //         for (int i = 0; i < no_chem; i++)
    //         {
    //             // cout << i << endl;
    //             vector1<int> jpow(nof);
    //             for (int k = i * (nof + 1) + 2; k < i * (nof + 1) + 2 + nof; k++)
    //             {

    //                 jpow[k - (i * (nof + 1) + 2)] = (int)mat1(j, k);
    //             }
    //             // cout << jpow << endl;

    //             GenericChemistry<double> c6_0(rate_multiplier*mat1(j, i * (nof + 1) + 1), jpow);
    //             c6.add_chemical_reaction(c6_0, i);
    //             double tot1 = 1.0;
    //             for (int k = 0; k < nof; k++)
    //             {
    //                 tot1 *= Power(init[k], jpow[k]);
                    
    //             }

    //             // cout << mat1(j, i * (nof + 1) + 1) <<","<< tot1 << "," <<mat1(j, i * (nof + 1) + 1) * tot1 << endl;

    //             tot += mat1(j, i * (nof + 1) + 1) * tot1;
                
                
    //         }
    //         // cout << tot << endl;
            

    //         my_chemsitry.add_method(c6, j - (5+n));
    //     }
    // }
    // pausel();
    a.set_chems(my_chemsitry);




    


    a.setup_matrices();

    vector<matrix<myc>> v;

    for (int j = 0; j < nof; j++)
    {
        matrix<myc> field1(p.N1, p.N2);
        v.push_back(field1);
    }

    double gt=0.6;

    for (int lk = 0; lk < nof; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1; i++)
        {
            //double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);

            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);

                v[lk](i, j) = x1 + gt  * r1;
            }
        }
    }
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
    // matrix<double> mat3 = importcsv(importstring1, T, err1);

    // string importstring2="/home/dino/External/Waves/chemistry19/field1res_i=3078_real.csv";
    // matrix<double> mat4 = importcsv(importstring2, T, err1);

    for (int lk = 0; lk < nof; lk++)
    {
        a.set_field(v[lk], lk);
    }
    // a.set_field(mat3, 0);
    // a.set_field(mat4, 1);

    a.calculate_initial_weight(SQR(10));


    cout << "calc" << endl;
    

    cout << "all fields set" << endl;

    // auto start = std::chrono::high_resolution_clock::now();
    int runtime = 20000;
    int every = 10;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);
    vector1<bool> ps2(nof, true);
    
    
    // // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runtime; i++)
    {

        if (i % every == 0 && i > 0 )
        {
            // stringstream strep1;
            // stringstream strep2;
            // stringstream strep3;
            // stringstream strep4;

            // strep1 << dens;
            // strep4 << c0;
            // strep2 << c1;
            // strep3 <<  surf;
            string s1 = "FHN";
            // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();
            string su = s1.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_some_results(su + s2,ps2);
        }
        cout << i << endl;
        cout << "begin" << endl;
        a.Update();
        bool chck = true;
        a.check_field(chck);
        if (!chck)
            break;
    }

    
    // ofstream myfile;
    // ofstream myfile2;
    // string fil1 = string("field0") + importstring;
    // string fil2 = string("field1") + importstring;
    // myfile.open(fil1.c_str());
    // myfile2.open(fil2.c_str());
    // for (int i = 0; i < runtime; i++)
    // {

    //     if (i % every == 0 && i > 00000)
    //     {

    //         for(int j = 0 ; j < p.N1-1 ; j++) {
    //             // cout << a.fields[0][j * p.N1] << endl;
    //             myfile << a.fields[0][j*p.N1] <<",";
    //         }
    //         myfile << a.fields[0][p.N1*p.N1-1] << endl;

    //         for (int j = 0; j < p.N1 - 1; j++) {
    //             myfile2 << a.fields[1][j * p.N1] << ",";
    //         }
    //         myfile2 << a.fields[1][p.N1 * p.N1 - 1] << endl;

    //         // pausel();
    //         // stringstream strep1;
    //         // stringstream strep2;
    //         // stringstream strep3;
    //         // stringstream strep4;

    //         // strep1 << dens;
    //         // strep4 << c0;
    //         // strep2 << c1;
    //         // strep3 <<  surf;
    //         // string s1 = importstring;
    //         // // // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
    //         // stringstream ss;
    //         // ss << setw(number_of_digits) << setfill('0') << i / every;
    //         // string s2 = "_i=" + ss.str();
    //         // string su = s1.substr(0, s1.size() - 4);
    //         // cout << su + s2 << endl;
    //         // a.print_some_results(su + s2,ps2);
    //     }
    //     cout << i << endl;
    //     cout << "begin" << endl;
    //     a.Update();
    //     bool chck = true;
    //     a.check_field(chck);
    //     if (!chck)
    //         break;
    // }
    
}