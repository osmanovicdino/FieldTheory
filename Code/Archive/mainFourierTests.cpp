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


matrix<double> make_sphere_inside(double avd, int N1, int N2) {
    matrix<double> mat(N1,N2);


    double radialsize = sqrt(N1*N2/(2*pi ) );
    double centrex = (double)(N1) / 2.;
    double centrey = (double)(N2) / 2.;
    double l = 2.0;
    for(int i = 0  ; i < N1 ; i++) {
        for(int j = 0 ; j < N2 ; j++) {
            // if(SQR(i-centrex) + SQR(j-centrey) < SQR(radialsize) ) {
            //     mat(i,j) =  2*avd;

            // }
            // else{
            //     mat(i,j) =  1E-10;
            // }
            double r = (double)rand()/(double)RAND_MAX;
            double r2 = 1.0 + 0.1*(2*r-1);
            mat(i,j) = 2.*avd*(tanh(- (sqrt(SQR(i - centrex) + SQR(j - centrey)) -radialsize) /l ) + 1) * r2;
        }
    }
    return mat;
}

matrix<double> make_sphere_outside(double avd, int N1, int N2)
{
    matrix<double> mat(N1, N2);

    double radialsize = sqrt(N1 * N2 / (2 * pi));
    double centrex = (double)(N1) / 2.;
    double centrey = (double)(N2) / 2.;
    double l = 2.0;
    for (int i = 0; i < N1; i++)
    {
        for (int j = 0; j < N2; j++)
        {
            // if (SQR(i - centrex) + SQR(j - centrey) < SQR(radialsize))
            // {
            //     mat(i, j) = 1E-10;
            // }
            // else
            // {
            //     mat(i, j) = 2 * avd;
            // }
            
            mat(i,j) = 2. * avd *(tanh((sqrt(SQR(i - centrex) + SQR(j - centrey)) - radialsize) / l) + 1);
        }
    }

    return mat;
}

// struct modelBnoise {

//     double operator()(double k1, double k2) {
//         return k1*k1+k2*k2;
//     }

// };

int main(int argc, char **argv)
{

    CH_builder p;
    int nof = 4;
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;

    typedef complex<double> myc;
    GenNoise<myc> a(p); //create a field

    double* str = new double [4];
    str[0] = 1.;
    str[1] = 1.;
    str[2] = 1.;
    str[3] = 1.;

    modelBnoise A;
    //a.GenFields(A,str,1./128.);

    CHWithNoise<modelBnoise> B(p,A,str);
 

    // srand(time(NULL));
    // int N = 128;
    // fftw_complex *in, *out;
    // fftw_plan p;
    // in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    // out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    // for(int i = 0  ; i < N ; i++) {
    //     in[i][0] = (double)rand()/(double)(RAND_MAX);
    //     in[i][1] = 0.0;
    // }

    // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(p); /* repeat as needed */

    // for(int i = 0  ; i < N ; i++) {
    //     cout << out[i][0] << "+" << out[i][1] << ",";
    // }
    // cout << endl;
    // fftw_destroy_plan(p);
    // fftw_free(in);
    // fftw_free(out);
}