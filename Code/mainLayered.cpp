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
#include "DensityFunctionalTheory/layered.h"

#include "fftw3.h"

using namespace std;

int main(int argc, char **argv)
{
    srand(time(NULL));

    double c1;
    double c2;
    double kT;
    string str;
    int sp;

    if (argc == 4)
    {
        c1 = atof(argv[1]);
        c2 = atof(argv[2]);
        kT = atof(argv[3]);

        // str = string(argv[7]);
        // sp = atof(argv[8]);
    }

    layered a;



    double w11 = -12.09 / kT;
    double w12 = -8.38 / kT;
    double w22 = -7.39 / kT;

    a.setchis(w12,w11,w22);

    cout << a.chi12 << endl;
    cout << a.chi13 << endl;
    cout << a.chi23 << endl;


    double dens1 = c1;
    double dens2 = c2;

    for(int i = 0  ; i < 128*128*128; i++) {
        double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
        double r2 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
        a.density1[i] = dens1 + 0.9* dens1*r1;
        a.density2[i] = dens2 + 0.9 * dens2 * r2;
    }

    vector1<double> prev1 = a.density1;
    vector1<double> prev2 = a.density2;

    for(int i = 0 ; i < 10000 ; i++) {
    cout << i << endl;
    a.update3D();
    
    }

    stringstream strep1;
    stringstream strep2;
    stringstream strep3;

    strep1 << c1;
    strep2 << c2;
    strep3 << kT;

    string s1 = "c1=" + strep1.str() + "_c2=" + strep2.str() + "_kT=" + strep3.str() + "_f1";
    string s2 = "c1=" + strep1.str() + "_c2=" + strep2.str() + "_kT=" + strep3.str() + "_f2";

    outfunc(a.density1,s1);

    outfunc(a.density2,s2);

    // for (int i = 0; i < 9900; i++)
    // {
    //     cout << i << endl;
    //     a.update3D();
    // }

    // string s3 = "c1=" + strep1.str() + "_c2=" + strep2.str() + "_kT=" + strep3.str() + "_f1f";
    // string s4 = "c1=" + strep1.str() + "_c2=" + strep2.str() + "_kT=" + strep3.str() + "_f2f";

    // outfunc(a.density1, s3);

    // outfunc(a.density2, s4);
}