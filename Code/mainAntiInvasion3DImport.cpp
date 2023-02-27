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
    double c0;
    double c1;
    double eps;
    double densi;
    string str;
    int sp;
    c0 = 0.1;
    c1 = 1.0;
    eps = 0.2;
    double Di;
    string str1,str2;
    string dirstr;

    // string dir = "/home/dino/Documents/Chemistry/Invasion/Resubmission/DataGood/";
    // str1 = "field0original=0.1_densi=1_eps1=0_D2=10_i=";
    // str2 = "field1original=0.1_densi=1_eps1=0_D2=10_i=";

    if (argc == 6)
    {
        //sp = atof(argv[1]);
        densi = atof(argv[1]);
        Di = atof(argv[2]);
        dirstr = argv[3];
        str1 = argv[4];
        str2 = argv[5];
    }
    else {
        error("not enough inputs");
    }



    string import1 = str1;
    string import2 = str2;
    // stringstream ss;
    // ss << std::setw(3) << std::setfill('0') << sp;
    // string choice = ss.str();
    // string ext= "_real.csv";

    string total1 = dirstr + import1;
    string total2 = dirstr + import2;


    cout << total1 << endl;

    cout << total2 << endl;

    double Td;
    bool err1;
    matrix<double> mat1 = importcsv(total1, Td, err1);

    matrix<double> mat2 = importcsv(total2, Td, err1);

    //count the total number of  sites that are less than the limit
    // int tp =0;
    // for(int i = 0  ; i < mat1.getnrows() ; i++) {
    //     for(int j = 0  ; j < mat1.getncols() ; j++) {
    //         if(mat1(i,j) < 0.55) tp++;
    //     }
    // }

    // int px = 128*128*128;

    // double frac = (double)px/(double)(tp);

    // double densi2 = densi*frac;


    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc, myc> FWCC;
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    //double p0 = 0.4;
    double rate1 = 1.0;
    // double A0 = atof(argv[2]);

    // double eps1 = -0.8;
    // double eps2 = 0.8;

    CH_builder p;
    p.number_of_fields = 3;
    p.dimension = 3;
    p.N1 = 128;
    p.N2 = 128;
    p.N3 = 128;

    CH<myc> a(p);

    //double rate1 = 0.15;
    // InvasionSelfRelease c1(rate1, rate2, 0, 1, 2);
    // InvasionSelfRelease c2(-rate1, -rate2, 0, 1, 2);
    // InvasionSelfRelease c3(rate1, rate2, 0, 1, 2);
    //NoWeight c4;
    //2 is the anti-invader
    InvasionLinear<myc> ch1(-rate1, 1, 2);
    InvasionLinear<myc> ch2(rate1, 1, 2);
    InvasionLinear<myc> ch3(rate1, 1, 2);

    FWCC my_chemsitry(p);

    my_chemsitry.add_method(ch1, 0);
    my_chemsitry.add_method(ch2, 1);
    my_chemsitry.add_method(ch3, 2);

    FWCC my_weights(p);

    double c00 = c0;
    double c11 = c1;
    double nu = 1.0;

    double eps1 = 0.0;
    double eps2 = 0.0;

    CahnHilliardWithCouplingWeightSQR<myc> d1(c00, c11, nu, eps1, eps2, 1, 0); //activated
    DiffDiffusiveWeightSQR<myc> d2(eps1, 0); //inactivated
    DiffusiveWeight<myc> d3(1.0);// invader

    // DiffusiveWeight d2(1.0);
    // DiffDiffusiveWeightSQR d3(eps1, 0);


    my_weights.add_method(d1, 0);
    my_weights.add_method(d2, 1);
    my_weights.add_method(d3, 2);

    a.set_chems(my_chemsitry);
    a.set_weights(my_weights);

    cout << "weights added" << endl;

    RWC my_rules(p);
    cout << "rule wrapper created" << endl;

    double dt = 0.005;
    double L = 40.0;
    double D = 10.0;
    double temp1 = SQR(2*pi / L);
    
    double D2 = Di;
    

    
    DiffusionWithSurfaceTension3D<myc> e1(p, dt, D, temp1,eps);
    DiffusionWithInteraction3D<myc> e2(p, dt, D, temp1);
    NormalDiffusion3D<myc> e3(p, dt, D2, temp1);
    

    cout << "created diffusion" << endl;

    my_rules.add_method(e1, 0);
    cout << "method1" << endl;
    my_rules.add_method(e2, 1);
    cout << "method2" << endl;
    my_rules.add_method(e3, 2);
    cout << "method3" << endl;

    cout << "added methods" << endl;

    a.set_rules(my_rules);

    cout << "done" << endl;

    // matrix<double> field1 = importcsv("ci.csv", T, err1);
    // matrix<double> field2 = importcsv("Ii.csv", T, err2);
    // matrix<double> field3 = importcsv("ρi.csv", T, err3);
    // matrix<double> field4 = importcsv("Ai.csv", T, err4);
    double T;
    bool err;
    //matrix<double> field1 = importcsv("/u/home/d/dinoo/FieldTheory/Code/InitialConditions/data174.csv", T, err);
    //matrix<double> field1 = importcsv("./InitialConditions/data174.csv", T, err);

    cd *a1 = new cd[p.N1 * p.N2 * p.N3];
    cd *a2 = new cd[p.N1 * p.N2 * p.N3];
    cd *a3 = new cd[p.N1 * p.N2 * p.N3];

        // double bc = 0.166;
    // double p0 = 0.596;

    double meanofinv = 0.0;
    for (int i = 0; i < p.N1*p.N2*p.N3; i++)
    {

            a1[i] = mat1[i];
            a3[i] = densi;
            a2[i] = mat2[i] * (1 + (0.2 * ((double)rand() / (double)RAND_MAX) - 0.1) );
            //meanofinv += field3(i,j).real();
            //field4(i, j) = A0 + (0.2 * ((double)rand() / (double)RAND_MAX) - 0.1);
        
    }

    // meanofinv/=SQR((double)(p.N1));

    // //cout << meanofinv << endl;
    // myc fac = (densi / meanofinv);
    // field3 *= fac;


    a.set_field(a1, 0);
    a.set_field(a2, 1);
    a.set_field(a3, 2);

    cout << "set fields" << endl;

    int runtime = 10000;
    int every = 100;

    cout << "diffusion: " << D2 << endl;


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
            stringstream strep1;
            stringstream strep2;
            stringstream strep3;
            stringstream strep4;

            strep2 << densi;
            strep3 << 0.0;
            strep4 << D2;

            string s1 = "_densi=" + strep2.str() + "_D2=" + strep4.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();

            a.print_all_results(s1 + s2);
        }
        cout << i << endl;
        a.Update();
    }


}