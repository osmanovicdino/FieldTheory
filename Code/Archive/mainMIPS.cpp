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
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "./DataStructures/basic.h"
#include "./DataStructures/vector1.h"
#include "./DataStructures/matrix2.h"
#include "./DataStructures/matrix2.cpp"
#include "./MD/potential.h"
//#include "intmatrix.h"
#include "./DensityFunctionalTheory/mips.h"

using namespace std;

int main(int argc, char **argv)
{
    // int N = 1000;
    // int nr = 128;
    // int nq = 40;
    // mips_mf a(nr, nq, N, 39.63, -0.03);
    // double tot = 0.0;


    // a.dat[0](0, 0) = 1;

    // vecmat b = a.calc_int2();

    // outfunc(b,"norm");

    // vecmat c(nq);

    // for(int q1 =  0 ; q1 < nq ; q1++ ) {
    //     c[q1]=matrix<double>(nr,nr);
    //     for(int j = 0 ; j < nr ; j++) {
    //         for(int k = 0  ; k < nr ; k++) {
    //             c[q1](j,k)  = a.potential(a.xpos[j], a.ypos[k], a.qpos[q1], 0, 0, 0);
    //         }
    //     }
    // }
    
    // outfunc(c,"norm2");

    
int N = 1000;
int nr =128;
int nq = 40;
mips_mf a(nr,nq,N,39.63,1.00);

double tot = 0.0;
for(int i = 0  ; i < nq ; i++) {
    for(int j = 0  ; j < nr ; j++) {
        for(int k = 0  ; k < nr ; k++) {
            double val = 1. + 0.1 * cos(2 * pi * j / (double)nr)+ 0.1 * cos(2 * pi * k / (double)nr);
            a.dat[i](j,k) = val;
            tot += (a.dq)*SQR(a.dx)*val; 
        }
    }
}


double norm = (double)((double)N / tot);
//cout << norm << endl;
(a.dat) *= norm;

//vecmat b = a.calc_int();


// outfunc(a.qintegral(),"int");


for(int i = 0; i < 19999; i++) {
pausel();

a.updatedatstore();


}


    /* 
matrix<COMPLEX> res2(128,128);
for(int i = 0  ; i < 128 ; i++) {
    for(int j  = 0 ; j < 128 ; j++) {
        COMPLEX b;
        b.real = cos(2*pi*i/128.);
        b.imag = 0.0;
        res2(i,j) = b;
    }
}

short checkl = FFT2D(res2, 128,128,1);



matrix<COMPLEX> mat(128,128);
for(int i = 0  ; i < 128 ; i++) {
    for(int j = 0 ; j < 128 ; j++) {
        mat(i,j).real = exp(-SQR(i-32)-SQR(j-32));
    }
}


FFT2D(mat,128,128,1);

matrix<COMPLEX> newm = mat&res2;

FFT2D(newm,128,128,-1);

matrix<double> final(128,128);
for(int i  = 0 ; i < 128 ; i++) {
    for(int j = 0  ; j < 128 ; j++) {
        final(i,j) = newm(i,j).getreal();
    }
}

cout << (128.*128.)*final << endl;

cout << "built" << endl;
pausel();
 */
}