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
    double c0 = 0.2;
    double c1 = 0.8;
    double eps = 0.4;

    // line 1 = all parameters;

    // line 2-5 = no. of terms -> rate, pow1, pow2, pow3, pow4 etc
    string importstring;

    if (argc == 2)
    {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();
    }
    else
    {
        error("no");
    }

    double T;
    bool err1;
    matrix<double> mat1 = importcsv(importstring, T, err1);

    int n = mat1(0, 0);

    vector1<int> phaseseps(n);

    for (int i = 0; i < n; i++)
        phaseseps[i] = mat1(1, i);

    vector1<double> epsi((n) * (n - 1) / 2);

    for (int i = 0; i < (n) * (n - 1) / 2; i++)
    {
        epsi[i] = -mat1(2, i);
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
        init[i] = mat1(3, i);

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
    int nof = 4;
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;



    //SurfBundle A(p,dt,diffusion_constant,temp1,eps,alpha);

    

    CHFracDt a(p);

    FWCC my_chemsitry(p);
    NoWeight<myc, myc> nw;
    for (int j = 4; j < mat1.getnrows(); j++)
    {
        int no_chem = mat1(j, 0);
        // cout << no_chem << endl;

        if (no_chem == 0)
        {
            my_chemsitry.add_method(nw, j - 4);
        }
        else
        {
            MultipleReactions<myc> c6(no_chem);

            double tot = 0.0;

            for (int i = 0; i < no_chem; i++)
            {
                // cout << i << endl;
                vector1<int> jpow(nof);
                for (int k = i * (nof + 1) + 2; k < i * (nof + 1) + 2 + nof; k++)
                {

                    jpow[k - (i * (nof + 1) + 2)] = (int)mat1(j, k);
                }
                // cout << mat1(j,i*(nof+1)+1) << endl;
                // cout << jpow << endl;
                // pausel();
                GenericChemistry<myc> c6_0(20. * mat1(j, i * (nof + 1) + 1), jpow);
                c6.add_chemical_reaction(c6_0, i);
                double tot1 = 1.0;
                for (int k = 0; k < nof; k++)
                {
                    tot1 *= Power(init[k], jpow[k]);
                    // cout << tot1 << ",";
                }
                // cout << endl;

                tot += mat1(j, i * (nof + 1) + 1) * tot1;

                // cout << endl;
            }
            cout << tot << endl;

            my_chemsitry.add_method(c6, j - 4);
        }
    }

    cout << "chemistries added" << endl;

    FWCC my_weights(p);

    double nu = 1.0;

    {
        int k = 0;
        vector1<int> which1(nof - 1);
        vector1<double> myeps(nof - 1);
        for (int i = 0; i < nof; i++)
        {
            if (i == 0)
            {
            }
            else
            {
                which1[k] = i;
                myeps(k) = epsa(0, i);
                k++;
            }
        }
        CahnHilliardWithCouplingWeightGenericN<myc> w0(c0, c1, nu, myeps, which1);
        my_weights.add_method(w0, 0);
    }

    for (int lk = 1; lk < nof; lk++)
    {
        int k = 0;
        vector1<int> which1(nof - 1);
        vector1<double> myeps(nof - 1);
        for (int i = 0; i < nof; i++)
        {
            if (i == lk)
            {
            }
            else
            {
                which1[k] = i;
                myeps(k) = epsa(lk, i);
                k++;
            }
        }

        DiffDiffusiveWeightGenericN<myc> w1(myeps, which1);

        my_weights.add_method(w1, lk);
    }

    // my_weights.add_method(w0, 0);

    // my_weights.add_method(w1, 1);

    // my_weights.add_method(w2, 2);

    // my_weights.add_method(w3, 3);

    cout << "weights added" << endl;

    a.set_chems(my_chemsitry);
    a.set_weights(my_weights);

    double dt = 0.05;
    // double dx = 0.05;

    // double temp1 = (1. / (dx * p.N1)) ;

    double surf = eps;
    double L = 50.0;
    double temp1 = SQR(2. * pii / L);

    // //cout << (1. / (dx * p.N1)) << endl;
    // cout << temp1 << endl;

    double diffusion_constant = 1.;
    double alpha = 0.9;

    cout << "creating surf bundle" << endl;
    SurfBundle A(p, dt, diffusion_constant, temp1, eps, alpha);
    a.AddBundleMethod(A,0);
    cout << "Bundle 1 added" << endl;
    for (int lk = 1; lk < nof; lk++)
    {
        cout << lk << endl;
        IntBundle B(p, dt, diffusion_constant, temp1, alpha);
        a.AddBundleMethod(B,lk);
    }

    cout << "bundles added" << endl;
    vector<matrix<myc>> v;

    for (int j = 0; j < nof; j++)
    {
        matrix<myc> field1(p.N1, p.N2);
        v.push_back(field1);
    }

    double gt = 0.6;

    for (int lk = 0; lk < nof; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = x1 + gt * x1 * r1;
            }
        }
    }

    for (int lk = 0; lk < nof; lk++)
    {
        a.set_field(v[lk], lk);
    }

    a.setupInitial();

    // cout << "done" << endl;

    // a.Update();
    cout << "all fields set" << endl;

    // auto start = std::chrono::high_resolution_clock::now();
    int runtime = 50001;
    int every = 10;

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

        if (i % every == 0 && i > 0)
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
            string su = s1.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_all_results(su + s2);
        }
        cout << i << endl;
        cout << "begin" << endl;
        a.Update();
    }
}