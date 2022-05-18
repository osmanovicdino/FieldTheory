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

    vector1<double> phasesepsparams(8);

    for (int i = 0; i < 8; i++)
        phasesepsparams[i] = mat1(1, i);

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
    int nof = n;
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;

    CHC a(p);

    for(int i = 0 ; i < nof ; i++) {
        for(int j = i+1  ; j < nof ; j++) {
            //cout << i << " " << j << endl;
             a.set_interaction(epsa(i,j), i, j);
        }
    }

    a.set_diffusion(phasesepsparams[0]);
    a.set_epsilon(phasesepsparams[1]);
    a.set_c0_c1(phasesepsparams[2],phasesepsparams[3]);
    double L = phasesepsparams[4];
    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);

    a.set_alpha(phasesepsparams[5]);
    a.set_dt(phasesepsparams[6]);

    double rate_multiplier = phasesepsparams[7];
    

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
                GenericChemistry<myc> c6_0(rate_multiplier*mat1(j, i * (nof + 1) + 1), jpow);
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
    a.set_chems(my_chemsitry);

    // pausel();

    


    a.setup_matrices();

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

    a.calculate_initial_weight();

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
        bool chck = true;
        a.check_field(chck);
        if (!chck)
            break;
    }
}