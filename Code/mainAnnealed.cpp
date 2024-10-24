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
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "./DataStructures/basic.h"
#include "./DataStructures/vector1.h"
#include "./DataStructures/matrix2.h"
#include "./DataStructures/matrix2.cpp"
#include "./MD/potential.h"
// #include "intmatrix.h"
//  #include "./MD/MD.h"
//  #include "./MD/Langevin.h"

// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

#include "./MeanField/IsingPol.h"
#include "./MeanField/IsingPolPartial.h"

using namespace std;




int main(int argc, char** argv) {
srand (time(NULL));


int ns = 3;

IsingPolymer a(100,ns);


double target[100] = {2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917, 2.61917,2.61917, 2.61918, 2.61918, 2.61918, 2.61918, 2.61918, 2.61918,2.61918, 2.61917, 2.61917, 2.61915, 2.61913, 2.6191, 2.61906,2.61901, 2.61896, 2.61892, 2.6189, 2.61894, 2.61908, 2.61938, 2.6199,2.62075, 2.62199, 2.62367, 2.62577, 2.62814, 2.6304, 2.6319, 2.63154,2.62773, 2.61826, 2.60026, 2.57009, 2.52332, 2.45463, 2.35763,2.22448, 2.04515, 1.80596, 1.48673, 1.05576, 0.474853, -0.00264318};

double u0 = 10.0;
double b = 1.0;

double c,c2;

double j11,j12,j13,j22,j23,j33;

if (argc == 7)
{
	j11 = atof(argv[1]);
	j12 = atof(argv[2]);
	j13 = atof(argv[3]);
	j22 = atof(argv[4]);
	j23 = atof(argv[5]);
	j33 = atof(argv[6]);
}
else
{
	error("no");
}

matrix<double> field(100,ns);
matrix<double> fixed_field(100,ns);

for(int i = 0 ; i < 100 ; i++) {
	for(int j = 0  ; j < ns ; j++) {
		field(i,j) = 0.;
	}
}

// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=2.71828+10./(1. + exp(-(50. - i)/1.));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-(target[i])) - exp(-field(i,0)));
// }


//double u0 =  pow(10.0,(double)ui);

a.setu0(u0);
a.setb(b);

double dv = 0.8;
//double c = 1.;

double v = log(2.*cosh(dv/2.));
// double jaa = -c*u0;
// double jab = c2*u0;

stringstream ss1,ss2,ss3,ss4,ss5;


ss1 << b;
ss2 << u0;
// ss3 << jaa;
// ss5 << jab;
ss4 << dv;

// Iterate through the matrix and build the string



for(int i = 0 ; i < 100 ; i++) {
	for(int j = 0  ; j < ns ; j++)
		fixed_field(i,j)=v-dv/2.+((double)j/(double)(ns-1))*dv;
}
// }
// for(int i = 0 ; i < 100 ; i++) {
// fixed_field(i,1)=v+dv/2.;
// }

a.setfixed_field(fixed_field);



matrix<double> V2(ns,ns);

V2(0,0) = j11;
V2(1, 1) = j22;
V2(2, 2) = j33;

V2(0, 1) = j12;
V2(1, 0) = j12;

V2(0, 2) = j13;
V2(2, 0) = j13;

V2(1, 2) = j23;
V2(2, 1) = j23;

string result;
for (int i = 0; i < V2.getnrows(); i++)
{
	for (int j = 0; j < V2.getncols(); j++)
	{
		result += std::to_string(V2(i,j));
		// Add a comma after each element except the last one
		if (i != V2.getnrows() - 1 || j != V2.getncols() - 1)
		{
			result += ",";
		}
	}
}

string bstring = "_b=" + ss1.str();
string uistring = "_ui=" + ss2.str();
string jstring = "_j=" + result;
// string jbstring = "_jab=" + ss5.str();
string dvstring = "_dv=" + ss4.str();

string string1 = "field" + bstring + jstring + uistring + dvstring;

string string2 = "den" + bstring + jstring + uistring + dvstring;

string string3 = "fe" + bstring + jstring + uistring + dvstring;

// for(int i = 0 ; i < ns ; i++) {
// 	for(int j = i ; j < ns ; j++) {
// 		// if(i == j) {
// 		// 	if(i==0) V2(i,j) = jaa*;
// 		// 	if(i==1) V2(i,j) = jaa;
// 		// }
// 		// else {
// 		// 	V2(i,j) =  0.;
// 		// 	V2(j,i) =  0.;
// 		// }
// 		// double sigma = ((double)i / (double)(ns - 1));
// 		// double sigma1 = ((double)j / (double)(ns - 1));

// 		// double gigma = sigma *sigma1 *jaa + sigma * (1 - sigma1) * jab + (1 - sigma) * sigma1 *jab ;
// 		// V2(i, j) = gigma;
// 		// V2(j, i) = gigma;
// 	}
// }
a.setV(V2);

// a.setfixed_field(fixed_field2);

// cout << a.freeenergy(field) << endl;
// cout << a.freeenergy(fixed_field2) << endl;
int runtime = 10000;
double stepsize = 0.01;
a.run(field,runtime,stepsize);

matrix<double> ft(1,1);

double fe = a.freeenergy(field);
ft(0,0)=fe;
matrix<double> den =a.densityfromfield(field);
outfunc(field,string1);
outfunc(den,string2);
outfunc(ft,string3);









return 0;
}
