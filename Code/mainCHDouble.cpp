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


    double x12;
    double x13;
    double x23;
    double den1;
    double den2;
    int cutoff;
    if (argc == 7)
    {
        x12=atof(argv[1]);
        x13=atof(argv[2]);
        x23=atof(argv[3]);
        den1=atof(argv[4]);
        den2=atof(argv[5]);
        cutoff=atof(argv[6]);
    }
    else
    {
        error("no");
    }

    CH_builder p;
    int nof = 2;
    p.number_of_fields = nof;
    p.N1 = 512;
    p.N2 = 512;

    CHD a(p);
    matrix<double> Dx(2,2);
    Dx(0,0)=1.;
    Dx(1, 1) = 1.;

    a.set_interaction_and_diffusion(x12,x13,x23,Dx);
    
    a.set_dt(0.005);

    double L = 100.;
    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);
    a.set_epsilon(0.4);

    a.setup_matrices();

    typedef complex<double> myc;
    vector<matrix<double>> v;

    for (int j = 0; j < nof; j++)
    {
        matrix<double> field1(p.N1, p.N2);
        v.push_back(field1);
    }
    vector1<double> init(2);
    init[0]=0.0;
    init[1]=0.0;

    double gt = 0.1;

    for (int lk = 0; lk < nof; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = x1 + gt * r1;
            }
        }
    }



    vector1<double> tots(nof);

    for (int lk = 0; lk < nof; lk++)
    {
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                
                tots[lk] += v[lk](i,j);
            }
        }
    }


    tots /= double(p.N1*p.N2);

    //double dens = -0.1;
    vector1<double> dens(2);
    dens[0]=den1;
    dens[1]=den2;
    for (int lk = 0; lk < nof; lk++)
    {
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {

                v[lk](i, j) -= tots[lk];
                v[lk](i,j) += dens[lk];
            }
        }
    }


    //we can cut off the high frequency modes

    

    for (int lk = 0; lk < nof; lk++)
    {
        a.set_field(v[lk], lk);
    }


    vector1<int> cutoff2(2,cutoff);

    a.setupInitial(cutoff2);

    int runtime = 50001;
    int every = 10000;

    string importstring = "twofields";
    stringstream ss1,ss2,ss3,ss4,ss5;
    ss1 << x12;
    ss2 << x23;
    ss3 << x13;
    ss4 << den1;
    ss5 << den2;
    importstring = importstring + string("x12=") + ss1.str() + string("_x13=") + ss2.str() + string("_x23=") + ss3.str() + string("_den1=")+ss4.str() + string("_den2=") + ss5.str();
    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runtime; i++)
    {

        if (i % every == 0)
        {
            // stringstream strep1;
            // stringstream strep2;
            // stringstream strep3;
            // stringstream strep4;

            // strep1 << dens;
            // strep4 << c0;
            // strep2 << c1;
            // strep3 <<  surf;
            string s1 = importstring;
            // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();
            string su = s1;//.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_all_results(su + s2);
        }
        cout << i << endl;
        cout << "begin" << endl;
        a.Update();
        bool chck = true;
        a.check_field(chck);
        if (!chck)
            break;
    }
}