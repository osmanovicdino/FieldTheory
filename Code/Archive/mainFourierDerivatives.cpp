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

//#include "complex.h"
#include <complex>
#include "fftw3.h" //before the rest of my code

#include "DataStructures/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "CahnHilliard/cahnhilliard.h"



using namespace std;
using namespace std::chrono;

int main(int argc, char **argv)
{
    srand(time(NULL));

    complex<double> a = 1.0;

    cout << log(a) << endl;

    /*
    CH_builder pg;
    int nof = 1;
    pg.number_of_fields = nof;
    pg.N1 = 1024;
    pg.N2 = 1024;
    //my 2D transform
    int N = 1024;
    complex<double>  * in;
    complex<double>  *out;
    fftw_plan p, p2;

    in = (complex<double> *)fftw_malloc(N * N * sizeof(complex<double>));

    out = (complex<double> *)fftw_malloc(N * N * sizeof(complex<double>));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            in[i * N + j] = {cos(6. * pii * (double)i / 1024.) * sin(6. * pii * (double)j / 1024.), 0.};
            out[i * N + j] = {0., 0.};
        }
    }



    p = fftw_plan_dft_2d(N, N, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

    // cout << "created plan" << endl;
    fftw_execute(p);
    // cout << "done 2" << endl;
    outfunc(out, "StartData", pg);
    cout << "done" << endl;


    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double k1, k2;
            if (i <= 512)
            {
                k1 = (2 * pii) * i;
            }
            else
            {
                k1 = (2 * pii) * (i - 0*N);
            }
            if (j <= 512)
            {
                k2 = (2 * pii) * j;
            }
            else
            {
                k2 = (2 * pii) * (j - 0*N);
            }
            complex<double> myfac = {-(SQR(k1) + SQR(k2)), 0.};

            out[i * N + j] *= myfac;
        }
    }

        p2 = fftw_plan_dft_2d(N, N, reinterpret_cast<fftw_complex *>(out), reinterpret_cast<fftw_complex *>(in), FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p2);

        outfunc(in,"Laplacian",pg);

        fftw_destroy_plan(p);
        fftw_destroy_plan(p2);
        fftw_free(in);
        fftw_free(out);
    */
        //     fftw_execute(p2);

        //     auto start = high_resolution_clock::now();
        //     for(int ti = 0 ; ti < 10000 ; ti++) {

        //     int N = 1024;
        //     complex<double> *in, *out;
        //     fftw_plan p, p2;

        //     complex<double> bg;

        //performing derivatives in k space

        //     in = (complex<double> *)fftw_malloc(sizeof(complex<double>) * N);

        //     out = (complex<double> *)fftw_malloc(sizeof(complex<double>) * N);
        //     for (int i = 0; i < N; i++)
        //     {
        //         in[i] = { cos(6 * pii * i / 1024.), 0. };
        //     }

        //     for (int i = 0; i < N; i++)
        //     {
        //         out[i] = { 0., 0. };

        //     }

        //     p = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

        //     fftw_execute(p); /* repeat as needed */
        //     for (int i = 0; i < 1024; i++)
        //     {
        //         if (i < 511)
        //         {
        //             double re = out[i].real();
        //             double im = out[i].imag();

        //             double newre = -im * 2 * pii * i;
        //             double newim = re * 2 * pii * i;

        //             out[i] = { newre, newim};

        //         }
        //         else if (i > 511)
        //         {
        //             double re = out[i].real();
        //             double im = out[i].imag();

        //             double newre = -im * 2 * pii * (i - 1024);
        //             double newim = re * 2 * pii * (i - 1024);

        //             out[i] = {newre, newim};
        //         }
        //         else
        //         {
        //             out[i] = {0, 0};
        //         }
        //     }

        //     p2 = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex *>(out), reinterpret_cast<fftw_complex *> (in), FFTW_BACKWARD, FFTW_ESTIMATE);

        //     fftw_execute(p2);

        // /*     cout << "real={";
        //     for (int i = 0; i < 1024; i++)
        //     {
        //         cout << (1. / 1024.) * in[i][0] << ",";
        //     }
        //     cout << "};" << endl;

        //     cout << "imag={";
        //     for (int i = 0; i < 1024; i++)
        //     {
        //         cout << (1. / 1024.) * in[i][1] << ",";
        //     }
        //     cout << "};" << endl;
        //     pausel(); */

        //     fftw_destroy_plan(p);
        //     fftw_free(in);
        //     fftw_free(out);
        //     }
        //     auto stop = high_resolution_clock::now();

        //     auto duration = duration_cast<microseconds>(stop - start);

        //     cout << "Time taken by function: "
        //          << duration.count() << " microseconds" << endl;

        //     auto start = high_resolution_clock::now();
        //     for (int ti = 0; ti < 10000; ti++)
        //     {
                int N = 1024;
                fftw_complex *in, *out;
                fftw_plan p,p2;

                in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

                out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
                for(int i = 0  ; i < N ; i++) {
                    in[i][0] = cos(6*pii*i/1024.);
                    in[i][1] = 0.;
                }

                for (int i = 0; i < N; i++)
                {
                    out[i][0] = 0.;
                    out[i][1] = 0.;
                }

                p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

                fftw_execute(p); /* repeat as needed */
                for (int i = 0; i < 1024; i++)
                {
                    if(i< 511) {
                        double re = out[i][0];
                        double im = out[i][1];

                        double newre = -im*2*pii*i;
                        double newim = re * 2 * pii * i;

                        out[i][0] = newre;
                        out[i][1] = newim;
                    }
                    else if (i > 511) {
                        double re = out[i][0];
                        double im = out[i][1];

                        double newre = -im * 2 * pii * (i-1024);
                        double newim = re * 2 * pii * (i-1024);

                        out[i][0] = newre;
                        out[i][1] = newim;
                    }
                    else {
                        out[i][0] = 0.0;
                        out[i][1] = 0.0;
                    }
                }

                p2 = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

                fftw_execute(p2);

                 cout << "real={";
                for(int i = 0  ; i < 1024 ; i++) {
                    cout << (1. / 1024. )* in[i][0] << ",";
                }
                cout << "};" << endl;

                cout << "imag={";
                for (int i = 0; i < 1024; i++)
                {
                    cout << (1. / 1024.) * in[i][1] << ",";
                }
                cout << "};" << endl;
                pausel(); 

        //         fftw_destroy_plan(p);
        //         fftw_free(in); fftw_free(out);
        //         }
        //     auto stop = high_resolution_clock::now();

        //     auto duration = duration_cast<microseconds>(stop - start);

        //     cout << "Time taken by function: "
        //          << duration.count() << " microseconds" << endl;

        // int N1 = 512;
        // int N2 = 512;
        // double **an_in_array = new double*[4];
        // double **an_out_array = new double*[4];
        // an_in_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));
        // an_out_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));

        // fftw_plan p;

        // p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

        // fftw_execute(p);

        /*
int N1 = 512;
int N2 = 512;



double *an_array;
an_array = (double *)fftw_malloc(N1 * N2 * sizeof(double));

double *out;
out = (double *)fftw_malloc(N1 * N2 * sizeof(double));
// in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
// out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

for (int i = 0; i < N1; i++)
{
    for (int j = 0; j < N2; j++)
    {
        an_array[i * N1 + j] = cos(pi*i*j/10.);
    }
}

fftw_plan p;

p = fftw_plan_r2r_2d(N1, N2, an_array, out, FFTW_REDFT10, FFTW_REDFT10, 1);

fftw_execute(p);
cout << an_array[0] << endl;
cout << out[0] / (4 * 512) << endl;




outfunc2D(out,N1,N2,"cosine");

fftw_plan p2;
p2 = fftw_plan_r2r_2d(N1, N2, out, an_array, FFTW_REDFT01, FFTW_REDFT01, 1);

fftw_execute(p2);

cout << an_array[0]/SQR(2*512) << endl;
cout << an_array[1] << endl;
cout << out[0] / (4*512) << endl;
*/

        // cout << a <<endl;
    }