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
    uint64_t microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // cout << microseconds_since_epoch << endl;
    // cout << time(NULL) << endl;
    int seed = microseconds_since_epoch % time(NULL);
    srand(seed);


    // if (argc == 8)
    // {
    //     x12=atof(argv[1]);
    //     x13=atof(argv[2]);
    //     x23=atof(argv[3]);
    //     den1=atof(argv[4]);
    //     den2=atof(argv[5]);
    //     rr=atof(argv[6]);
    //     wh = atof(argv[7]);
    // }
    // else
    // {
    //     error("no");
    // }
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

    // double x12 = 4.;
    // double x13 = 4.;
    // double x23 = 4.;
    // double den1 = 0.;
    // double den2 = 0.;
    // int rr = 1000;
    // int wh = 1 ;

    double T;
    bool err1;
    matrix<double> mat1 = importcsv(importstring, T, err1);//define the chemistry;

    CH_builder p;
    int nof = mat1(0,0);
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;

    //vector1<double> phasesepparams(8);
    double x12 = mat1(1,0);
    double x13 = mat1(1, 1);
    double x23 = mat1(1, 2);
    double den1 = mat1(1, 3);
    double den2 = mat1(1, 4);
    int rr = mat1(1, 5);
    int wh = mat1(1, 6);

    matrix<double> dm(nof,nof);
    for(int i = 0  ; i < nof ; i++) {
        dm(i,i) = 1.;
    }

    CHD a(p);

    a.set_interaction_and_diffusion(x12,x13,x23,dm);
    
    double dt = mat1(1,7);
    double L = mat1(1,8);
    double epsilon = mat1(1,9);

    a.set_dt(dt);

    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);
    a.set_epsilon(epsilon);

    for(int i = 0 ; i < nof ; i++) {
        for(int j = 2 ; j < nof ; j++) {
        a.set_interaction(0.,i,j);
        }
    }

    typedef complex<double> myc;
    typedef Field_Wrapper<myc, myc> FWCC;

    vector1<double> useshifts(nof);
    useshifts[0] = 1./3.;
    useshifts[1] = 1./3.;

    double rm =0.1;

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
                GenericChemistry<myc> c6_0(rm*mat1(j, i * (nof + 1) + 1), jpow);
                c6_0.set_shifts(useshifts);
                c6.add_chemical_reaction(c6_0, i);

                // cout << endl;
            }

            my_chemsitry.add_method(c6, j - 4);
        }
    }
    a.set_chems(my_chemsitry);

    a.setup_matrices();

    typedef complex<double> myc;
    vector<matrix<myc>> v;

    for (int j = 0; j < nof; j++)
    {
        matrix<myc> field1(p.N1, p.N2);
        v.push_back(field1);
    }
    vector1<double> init(2);
    init[0]=0.0;
    init[1]=0.0;

    int reduced_fieldno = 2;

    double gt = 0.1;

    for (int lk = 0; lk < reduced_fieldno; lk++) //we only change the first 2
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

    for (int lk = 0; lk < reduced_fieldno; lk++)
    {
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                tots[lk] += v[lk](i,j).real();
            }
        }
    }


    tots /= double(p.N1*p.N2);

    //double dens = -0.1;
    int which_first = wh;
    vector1<double> init2(reduced_fieldno);
    init2[0] = den1;
    init2[1] = den2;
    double zeroden=-1./3.;

    vector1<double> dens(2);
    int wh2 = wh == 1 ? 0 : 1;


    dens[wh2]=init2[wh2];
    dens[wh] = zeroden;

    for (int lk = 0; lk < reduced_fieldno; lk++)
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

    

    for (int lk = 0; lk < reduced_fieldno; lk++)
    {
        a.set_field(v[lk], lk);
    }


    matrix<myc> fieldtemp(p.N1, p.N2);

    vector1<int> cutoff(4,100);
    // cutoff[0]=100;
    // cutoff[1]=100;
    a.setupInitial(cutoff);

    // do the initial updating without chemistry

    int runtime = rr;
    int every = 100;

    string exportstring = "twofields";
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
            string s1 = exportstring;
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


    //inject the other material and turn on the chemistry



    //check the field variable 

    bool inject_away = true;
    double x1 = init2[wh];
    int totplaces=0;
    for (int i = 0; i < p.N1; i++)
    {
        for (int j = 0; j < p.N2; j++)
        {
            double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            if(inject_away && a.fields[wh2][i*p.N1+j] < 0.0) {
            fieldtemp(i, j) = x1 + gt * r1;
            totplaces++;
            }
            else{
                fieldtemp(i,j) =zeroden;
            }
        }
    }

    double profactor = (double)(p.N1*p.N2-totplaces)/((double)totplaces);

    for (int i = 0; i < p.N1; i++)
    {
        for (int j = 0; j < p.N2; j++)
        {
            double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            if (inject_away && a.fields[wh2][i * p.N1 + j] < 0.0)
            {
                fieldtemp(i, j) = x1 * (1 + profactor) - profactor * (zeroden) + gt * r1;
            }
            else if(!inject_away) {
                fieldtemp(i, j) = x1 + gt * r1;
            }
            else
            {
                fieldtemp(i, j) = zeroden;
            }
        }
    }

    // cout << wh << " " << wh2 << endl;
    
    a.set_field(fieldtemp, wh);
    // outfunc2D(a.fields[0],p.N1,p.N2,"b0");
    // outfunc2D(a.fields[1],p.N1,p.N2 ,"b1");

    // pausel();

    vector1<double> otherfields(nof);
    for(int i = 2 ; i < nof ; i++ ) {
        otherfields[i] = mat1(3,i-2);
    }

    // now we add the other material as well:
    for (int k = 2; k < nof; k++)
    {
        double x1 = otherfields[k];
        for (int i = 0; i < p.N1; i++)
        {

            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                fieldtemp(i, j) = x1 + gt * x1 * r1;
            }
        }
        a.set_field(fieldtemp, k);
    }

    // vector1<double> otherfields(nof);
    // otherfields[2] =0.4;
    // otherfields[3] =0.4;

    // //now we add the other material as well:
    // for(int k = 2 ; k < nof ; k++) {
    //     double x1 = otherfields[k];
    //     for (int i = 0; i < p.N1; i++)
    //     {

    //         for (int j = 0; j < p.N2; j++)
    //         {
    //             double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
    //             fieldtemp(i, j) = x1 + gt * x1 * r1;
    //         }
    // }
    //     a.set_field(fieldtemp, k);
    // }

    int newruntime = runtime+5000;
    tf = ceil((double)newruntime / (double)every);
    number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    for (int i = runtime; i < newruntime; i++)
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
            string su = s1.substr(0, s1.size() - 4);
            cout << su + s2 << endl;
            a.print_all_results(su + s2);
        }
        cout << i << endl;
        cout << "begin" << endl;
        a.Update_With_Chem();
        bool chck = true;
        a.check_field(chck);
        if (!chck)
            break;
    }
}