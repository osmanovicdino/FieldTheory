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
#include "./MD/MD.h"
#include "./MD/Langevin.h"

// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

#include "./MeanField/IsingPol.h"
#include "./MeanField/IsingPolPartial.h"

using namespace std;

int main(int argc, char **argv)
{
    srand(time(NULL));
    int ns = 1;

    IsingPolymer a(100, ns, 1);

    // for(int i = 0 ; i < 100 ; i++) {
    // 	field[i]=0.01*i;
    // }
    // vector1<double> den = a.radialsolver(field);

    // 	for(int j = 0  ; j < 2 ; j++) {
    // 		field(i,j)=sqrt((i+1)*(j+1));
    // 	}
    // }

    // cout << a.radialsolver(fieldv) << endl;
    // cout << a.radialsolver2(fieldv) << endl;

    //double u0=10.0;
    double v0 = 1.0;

    //a.setv0(0.0);
    a.setv0(v0);

    a.setR(10.0);

    double stepsize = 0.001;

    a.setstot(10000.0);
    a.settotpolymers(1.0);

    vector1<double> dat(100);

    for (int i = 0; i < 100; i++)
    {
        dat[i] = 0.0;
    }

    //matrix<double> den =a.densityfromfield(field);
    //outfunc(den,"den");

    int runtime = 1000;
    double dv = 0.0;
    
    //a.setu0(10.0);

    //a.radialsolver_copolymer()
    vector1<double> field1(100);
    vector1<double> field2(100);

    int s = 100;


        int size = 10;
        vector1<bool> seq(10001);

        // for (int i = 0; i < 100-1; i++)
        // {
        //     if ((i > s && i < s + size) || (i > 1000 - s - size && i < 1000 - s))
        //         seq[i] = true;
        //     //if(i%10==0) seq[i]=true;

        //     else
        //         seq[i] = false;
        // }
        for(int i = 0 ; i < 10001 ; i++) {
            int j = (i) % s;
            if(j <= size || j >= s-size) {
                seq[i] = true;
                //cout << i << endl;
            }
        }
        //cout << seq << endl;
       // pausel();

        //cout << seq << endl;
        a.setseq(seq);

 
        double startv = -1.0;

        matrix<double> field(100, ns + 1, 0.0);
        for (int i = 0; i < 100; i++)
        {
            field(i, 0) = 0.+0.001*(double)rand()/(double)RAND_MAX;
            field1[i] = field(i,0);
            field(i, 1) = 0.+0.001*(double)rand() / (double)RAND_MAX;

            field2[i] =  field(i,1);
        }

        // double norm;
        // matrix<double> testden = a.radialsolver_copolymer(field1,field2,norm);

        // cout << norm << endl;
        // cout << testden << endl;
        

        //pausel();
        matrix<double> fixed_field(100, ns + 1);

        double u0 = 10.0;
        double b = 0.02;

        a.setu0(u0);
        a.setb(b);

        //double dv = 1.00;

        double v = 0.;
        //double jaa = -c*u0;

        matrix<double> V2(ns + 1, ns + 1);
        double jaa = 9.0;
        V2(0, 0) = -jaa;

        V2(0, 1) = jaa;
        V2(1, 0) = jaa;


        V2(1, 1) = -jaa;
 

        //cout << V2 << endl;

        //double dv = 10.0;

        for (int i = 0; i < 100; i++)
        {
            fixed_field(i, 0) = v - dv / 2.;
        }
        for (int i = 0; i < 100; i++)
        {
            fixed_field(i, 1) = v + dv / 2.;
        }


        // for(int i = 0 ; i < ns+1 ; i++) {
        // 	for(int j = i ; j < ns+1 ; j++) {
        // 		if(i == j) {
        // 			if(i==0) V2(i,j) = 0.0;
        // 			if(i==1) V2(i,j) = jaa;
        // 		}
        // 		else {
        // 			V2(i,j) =  0.;
        // 			V2(j,i) =  0.;
        // 		}
        // 	}
        // }
        a.setV(V2);
        a.setfixed_field(fixed_field);



        //iterate fields

        
        stringstream ss1, ss2, ss3, ss4, ss5;

        ss1 << b;
        ss2 << u0;
        ss3 << jaa;
        ss4 << dv;
        ss5 << size;

        string bstring = "_b=" + ss1.str();
        string uistring = "_ui=" + ss2.str();
        string jstring = "_jaa=" + ss3.str();
        string dvstring = "_dv=" + ss4.str();
        string sizestring = "_size=" + ss5.str();

        string string1 = "field" + bstring + uistring + jstring + dvstring + sizestring;

        string string2 = "den" + bstring + uistring + jstring + dvstring + sizestring;

        string string3 = "fe" + bstring + uistring + jstring + dvstring + sizestring;

        // a.setfixed_field(fixed_field2);

        // cout << a.freeenergy(field) << endl;
        // cout << a.freeenergy(fixed_field2) << endl;

        a.run_copolymer(field, runtime, stepsize);

        matrix<double> ft(1, 1);

        double fe = a.freeenergy_copolymer(field);
        ft(0, 0) = fe;
        matrix<double> den = a.densityfromfield_copolymer(field);
        outfunc(field, string1);
        outfunc(den, string2);
        outfunc(ft, string3);
        
}