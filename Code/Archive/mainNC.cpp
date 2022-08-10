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

#include "basic.h"
#include "vector1.h"
#include "matrix2.h"
#include "matrix2.cpp"
#include "potential.h"
//#include "intmatrix.h"
#include "MD.h"
#include "Langevin.h"


// #include "BrownianGel.cpp"
// #include "BrownianGel2.cpp"
// #include "LangevinGel.cpp"
// #include "LangevinGelFixed.cpp"

#include "IsingPol.h"
#include "IsingPolPartial.h"


using namespace std;




int main(int argc, char** argv) {
srand (time(NULL));
/*
IsingPolymerPartial a(100,2);

//vector1<double> field(100);

double u0=10.0;
double v0=1.0;

a.setu0(u0);
a.setv0(v0);



a.setR(10.0);

double stepsize =  0.01;

a.setstot(1000.0);
a.settotpolymers(1.0);


int s = 250;
int size = 100;

vector1<bool> seq(1001);

// for(int i = 0; i < 1001 ; i++) {
// 	if((i>s&&i<s+size)||(i>1000-s-size&&i<1000-s)) seq[i]=true;
// 	//if(i%10==0) seq[i]=true;

	
// 	else seq[i]=false;
// }

for(int i = s-size/2  ; i < s + size/2 ; i++) {
	seq[i] =  true;
}
for(int i = 1001-s-size/2  ; i < 1001-s + size/2 ; i++) {
	seq[i] =  true;
}

a.setseq(seq);

double b = 1.0;

a.setb(b);

//double dv = 1.00;


double jaa = -2.0;
double dv1 = 1.0;
int ns = 1;
double v1 = log(2.*cosh(dv1/2.));


//double jaa = -c*u0;

matrix<double> fixed_field(100,ns+ns);


matrix<double> V2(ns+ns,ns+ns);


//0 and 2 are the bookmarks on
V2(1,1) = jaa;

//cout << V2 << endl;

//double dv = 10.0;

for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,0)=v1-dv1/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,1)=v1+dv1/2.;
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

stringstream ss1,ss2,ss3,ss4,ss5;


ss1 << b;
ss2 << u0;
ss3 << jaa;
ss4 << dv1;
ss5 << size;

string bstring = "_b="+ss1.str();
string uistring = "_ui="+ss2.str();
string jstring  = "_jaa="+ss3.str();
string dvstring = "_dv="+ss4.str();
string sizestring = "_size="+ss5.str();

string string1 = "field"+bstring+uistring+jstring+dvstring+sizestring;

string string2 = "den"+bstring+uistring+jstring+dvstring+sizestring;

string string3 = "fe"+bstring+uistring+jstring+dvstring+sizestring;

// a.setfixed_field(fixed_field2);

// cout << a.freeenergy(field) << endl;
// cout << a.freeenergy(fixed_field2) << endl;
int runtime = 100000;
vector1<double> fixed(100,-log(2.));

// for(int i = 0 ; i < 100 ; i++) {
// 	field(i,0)=1.*(double)rand()/(double)RAND_MAX;
// }
vector1<double> f2(100,0.5);


a.run_mixed();
*/




// matrix<double> fr(100,2);

// a.run_fixedseq(fr,runtime,stepsize);

/* WEAK BOOKMARKS
int ns = 2;

IsingPolymer a(100,ns,2);

double u0=10.0;
double v0=1.0;

a.setu0(u0);
a.setv0(v0);



a.setR(10.0);

double stepsize =  0.01;

a.setstot(1000.0);
a.settotpolymers(1.0);


int s = 150;
int size = 100;

vector1<bool> seq(1001);

for(int i = 0; i < 1001 ; i++) {
	if((i>s&&i<s+size)||(i>1000-s-size&&i<1000-s)) seq[i]=true;
	//if(i%10==0) seq[i]=true;

	
	else seq[i]=false;
}



a.setseq(seq);

double b = 1.0;

a.setb(b);

//double dv = 1.00;


double jaa = 0.0;
double dv1 = 1.0;
double dv2 = 1.0;

double v1 = log(2.*cosh(dv1/2.));

double v2 = log(2.*cosh(dv2/2.));
//double jaa = -c*u0;

matrix<double> fixed_field(100,ns+ns);


matrix<double> V2(ns+ns,ns+ns);


//0 and 2 are the bookmarks on
V2(1,1) = jaa;

V2(3,3) = jaa;

V2(1,3) = jaa;

V2(3,1) = jaa;

//cout << V2 << endl;

//double dv = 10.0;

for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,0)=v1-dv1/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,1)=v1+dv1/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,2)=v2-dv2/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,3)=v2+dv2/2.;
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

stringstream ss1,ss2,ss3,ss4,ss5;


ss1 << b;
ss2 << u0;
ss3 << jaa;
ss4 << dv1;
ss5 << size;

string bstring = "_b="+ss1.str();
string uistring = "_ui="+ss2.str();
string jstring  = "_jaa="+ss3.str();
string dvstring = "_dv="+ss4.str();
string sizestring = "_size="+ss5.str();

string string1 = "field"+bstring+uistring+jstring+dvstring+sizestring;

string string2 = "den"+bstring+uistring+jstring+dvstring+sizestring;

string string3 = "fe"+bstring+uistring+jstring+dvstring+sizestring;

// a.setfixed_field(fixed_field2);

// cout << a.freeenergy(field) << endl;
// cout << a.freeenergy(fixed_field2) << endl;
int runtime = 10000;
matrix<double> field(100,ns+ns,0.0);
// for(int i = 0 ; i < 100 ; i++) {
// 	field(i,0)=1.*(double)rand()/(double)RAND_MAX;
// }


a.run_dcopolymer(field,runtime,stepsize);

matrix<double> ft(1,1);

double fe = a.freeenergy_dcopolymer(field);
ft(0,0)=fe;
matrix<double> den = a.densityfromfield_dcopolymer(field);
outfunc(field,string1);
outfunc(den,string2);
outfunc(ft,string3);
*/


/* COPOLYMER SIMULATIONS
int ns = 2;

IsingPolymer a(100,ns,1);

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
double v0=1.0;

//a.setv0(0.0);
a.setv0(v0);



a.setR(10.0);

double stepsize =  0.01;

a.setstot(1000.0);
a.settotpolymers(1.0);


vector1<double> dat(100);

for(int i = 0 ; i < 100 ; i++) {
	dat[i]=0.0;
}








//matrix<double> den =a.densityfromfield(field);
//outfunc(den,"den");



int runtime = 1000;


//a.setu0(10.0);


//a.radialsolver_copolymer()
vector1<double> field1(100);
vector1<double> field2(100);




int s = 150;

for(int size = 70 ; size < 71 ; size += 20){
	for(double dv = 0.8 ; dv < 1.01 ; dv += 0.1) {
		for(double jaa = 0.0 ; jaa > -9.01 ; jaa -= 1.0) {




vector1<bool> seq(1001);

for(int i = 0; i < 1001 ; i++) {
	if((i>s&&i<s+size)||(i>1000-s-size&&i<1000-s)) seq[i]=true;
	//if(i%10==0) seq[i]=true;

	
	else seq[i]=false;
}



//cout << seq << endl;
a.setseq(seq);

matrix<double> field(100,ns+1,0.0);
for(int i = 0 ; i < 100 ; i++) {
	field(i,0)=1.*(double)rand()/(double)RAND_MAX;
}


matrix<double> fixed_field(100,ns+1);

double u0 =  10.0;
double b = 2.0;

a.setu0(u0);
a.setb(b);

//double dv = 1.00;

double v = log(2.*cosh(dv/2.));
//double jaa = -c*u0;

matrix<double> V2(ns+1,ns+1);



V2(0,0)=jaa;

V2(0,1)=0.0;
V2(1,0)=0.0;

V2(0,2)=jaa;
V2(2,0)=jaa;


V2(1,1)=0.0;
V2(2,1)=0.0;
V2(1,2)=0.0;


V2(2,2)=0.0;

//cout << V2 << endl;

//double dv = 10.0;

for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,1)=v-dv/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,2)=v+dv/2.;
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

stringstream ss1,ss2,ss3,ss4,ss5;


ss1 << b;
ss2 << u0;
ss3 << jaa;
ss4 << dv;
ss5 << size;

string bstring = "_b="+ss1.str();
string uistring = "_ui="+ss2.str();
string jstring  = "_jaa="+ss3.str();
string dvstring = "_dv="+ss4.str();
string sizestring = "_size="+ss5.str();

string string1 = "field"+bstring+uistring+jstring+dvstring+sizestring;

string string2 = "den"+bstring+uistring+jstring+dvstring+sizestring;

string string3 = "fe"+bstring+uistring+jstring+dvstring+sizestring;

// a.setfixed_field(fixed_field2);

// cout << a.freeenergy(field) << endl;
// cout << a.freeenergy(fixed_field2) << endl;



a.run_copolymer(field,runtime,stepsize);

matrix<double> ft(1,1);

double fe = a.freeenergy_copolymer(field);
ft(0,0)=fe;
matrix<double> den = a.densityfromfield_copolymer(field);
outfunc(field,string1);
outfunc(den,string2);
outfunc(ft,string3);

}
}
}

*/


//a.run(field,runtime,0.1);

//a.run(field,runtime,0.001);


// cout << a.freeenergy(fixed_field2) << endl;


// for(double b = 1.0 ; b < 3.1 ; b+=1.0)
// 	for(int ui = 0 ; ui < 2 ; ui++)
// 		for(double c = 0.0 ; c < 0.91 ; c+= 0.1 )
// 			for( double dv  = 0.0 ; dv < 1.01 ; dv += 0.1) {



int ns = 2;

IsingPolymer a(100,ns);


double u0 = 10.0;
double b = 1.0;


matrix<double> field(100,ns);
matrix<double> fixed_field(100,ns);


for(int i = 0 ; i < 100 ; i++) {
field(i,0)=1./(1. + exp((50. - i)/2.));
}
for(int i = 0 ; i < 100 ; i++) {
field(i,1)=-log(1. - exp(-field(i,0)));
}

//double u0 =  pow(10.0,(double)ui);

a.setu0(u0);
a.setb(b);

double dv = 0.5;
double c = 0.5;

double v = log(2.*cosh(dv/2.));
double jaa = -c*u0;

stringstream ss1,ss2,ss3,ss4;


ss1 << b;
ss2 << u0;
ss3 << jaa;
ss4 << dv;

string bstring = "_b="+ss1.str();
string uistring = "_ui="+ss2.str();
string jstring  = "_jaa="+ss3.str();
string dvstring = "_dv="+ss4.str();

string string1 = "field"+bstring+uistring+jstring+dvstring;

string string2 = "den"+bstring+uistring+jstring+dvstring;

string string3 = "fe"+bstring+uistring+jstring+dvstring;


for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,0)=v-dv/2.;
}
for(int i = 0 ; i < 100 ; i++) {
fixed_field(i,1)=v+dv/2.;
}

a.setfixed_field(fixed_field);



matrix<double> V2(ns,ns);
for(int i = 0 ; i < ns ; i++) {
	for(int j = i ; j < ns ; j++) {
		if(i == j) {
			if(i==0) V2(i,j) = 0.0;
			if(i==1) V2(i,j) = jaa;
		}
		else {
			V2(i,j) =  0.;
			V2(j,i) =  0.;
		}
	}
}
a.setV(V2);

// a.setfixed_field(fixed_field2);

// cout << a.freeenergy(field) << endl;
// cout << a.freeenergy(fixed_field2) << endl;
double stepsize = 0.01;
a.run(field,runtime,stepsize);

matrix<double> ft(1,1);

double fe = a.freeenergy(field);
ft(0,0)=fe;
matrix<double> den =a.densityfromfield(field);
outfunc(field,string1);
outfunc(den,string2);
outfunc(ft,string3);








// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=0.01+5./(1. + exp((90. - i)/1.));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-0.0) - exp(-field(i,0)));
// }

// matrix<double> den=a.densityfromfield(field);
// outfunc(den,"den");
// pausel();
// int runtime = 10000;
// a.run(field,runtime,0.1);
// den=a.densityfromfield(field);
// outfunc(den,"den");

// cout << a.freeenergy(field) << endl;
// pausel();
// a.setstot(100.0);
// a.settotpolymers(10.0);
// matrix<double> field2(100,2);
// for(int i = 0 ; i < 100 ; i++) {
// field2(i,0)=0.01+5./(1. + exp((90. - i)/1.));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field2(i,1)=-log(exp(-0.0) - exp(-field2(i,0)));
// }
// den=a.densityfromfield(field2);
// outfunc(den,"den2");
// pausel();
// a.run(field2,runtime,0.1);
// den=a.densityfromfield(field2);
// cout << a.freeenergy(field2) << endl;
// outfunc(den,"den2");
// pausel();


// a.setR(20.0);
// a.setstot(8000.);

// vector1<double> fieldv(100,0.0);
// double norm1 = 0.0;
// vector1<double> den = a.radialsolver(fieldv,norm1);

// vector1<double> fac = den&fieldv;
// cout << norm1 << endl;
// cout << a.integrate(fac) << endl;

// double unmixed = a.freeenergy(field);
// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=1./(1. + exp((50. - i)/2.));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(1. - exp(-field(i,0)));
// }
//cout << a.freeenergy(field) <<endl;

//cout << unmixed << endl;

// int runtime = 10000;
// a.setv0(0.0);
// a.run(field,runtime,0.001);
// cout << field << endl;

// int ju = 40;
// double vfac = ju*0.02*10.0;
// a.setv0(vfac);
// stringstream ss;
// ss<<ju;
// string endf=ss.str();
// endf+=".csv";

// string fi1 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/field_mixed"+endf;
// string de1 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/den_mixed"+endf;
// string fe1 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/freeenergy_mixed"+endf;

// string fi2 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/field_demixed"+endf;
// string de2 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/den_demixed"+endf;
// string fe2 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/freeenergy_demixed"+endf;

// string fi3 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/field_homo"+endf;
// string de3 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/den_homo"+endf;
// string fe3 = "/home/dino/External/Code/PolycombCode/sim-20-04-02-23:12:18/freeenergy_homo"+endf;

// double T;
// bool err1, err2, err3;

// matrix<double> field_mixed = importcsv(fi1,T,err1);
// matrix<double> field_demixed = importcsv(fi2,T,err2);
// matrix<double> field_homo = importcsv(fi3,T,err3);

// cout << a.freeenergy(field_mixed) << endl;
// cout << a.freeenergy(field_demixed) << endl;
// cout << a.freeenergy(field_homo) << endl;
// a.setv0(5.4);
// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=dat[i]+log(2.);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// double temp = 0.0;
// matrix<double> den = a.densityfromfield(field,temp);
// cout << a.freeenergy(field) << endl;
// outfunc(den,"den");

// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=5./(1. + exp((90. - i)/1.));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }

// den = a.densityfromfield(field,temp);
// cout << a.freeenergy(field) << endl;
// outfunc(den,"den");
// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=5.;
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// den = a.densityfromfield(field,temp);
// cout << a.freeenergy(field) << endl;
// outfunc(den,"den");
// pausel();

//double dat[100]={-0.000058, -0.000058, -0.000058, -0.000058, -0.000058, -0.000058, -0.000058, -0.000058, -0.000057, -0.000057, -0.000057, -0.000056, -0.000056, -0.000056, -0.000055, -0.000055, -0.000054, -0.000053, -0.000053, -0.000052, -0.000051, -0.00005, -0.000049, -0.000048, -0.000047, -0.000046, -0.000044, -0.000043, -0.000042, -0.00004, -0.000038, -0.000037, -0.000035, -0.000033, -0.000031, -0.000029, -0.000026, -0.000024, -0.000022, -0.000019, -0.000017, -0.000015, -0.000012, -0.00001, -8.e-6, -6.e-6, -4.e-6, -2.e-6, -1.e-6, -1.e-6, 0., -1.e-6, -2.e-6, -4.e-6, -8.e-6, -0.000012, -0.000019, -0.000027, -0.000037, -0.000049, -0.000064, -0.000083, -0.000105, -0.00013, -0.000161, -0.000196, -0.000237, -0.000285, -0.000339, -0.000401, -0.000472, -0.000552, -0.000642, -0.000743, -0.000856, -0.000982, -0.00112, -0.001273, -0.00144, -0.001622, -0.001818, -0.002029, -0.002255, -0.002493, -0.002744, -0.003004, -0.003273, -0.003546, -0.003822, -0.004095, -0.004362, -0.004618, -0.004858, -0.005078, -0.005274, -0.00544, -0.005573, -0.005671, -0.00573, -0.00575};
//double dat[100]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//double dat[100]={-0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428461805599453, -0.5428461805599453, -0.5428461805599453, -0.5428461805599453, -0.5428461805599453, -0.5428461805599453, -0.5428451805599454, -0.5428451805599454, -0.5428451805599454, -0.5428451805599454, -0.5428451805599454, -0.5428451805599454, -0.5428451805599454, -0.5428461805599453, -0.5428461805599453, -0.5428461805599453, -0.5428471805599453, -0.5428471805599453, -0.5428471805599453, -0.5428481805599453, -0.5428481805599453, -0.5428471805599453, -0.5428471805599453, -0.5428461805599453, -0.5428451805599454, -0.5428421805599454, -0.5428401805599453, -0.5428361805599453, -0.5428321805599453, -0.5428271805599453, -0.5428221805599454, -0.5428171805599453, -0.5428141805599453, -0.5428121805599453, -0.5428121805599453, -0.5428171805599453, -0.5428271805599453, -0.5428441805599453, -0.5428681805599452, -0.5429011805599454, -0.5429411805599452, -0.5429861805599453, -0.5430311805599454, -0.5430701805599454, -0.5430881805599453, -0.5430721805599453, -0.5429981805599453, -0.5428441805599453, -0.5425831805599454, -0.5421941805599453, -0.5416721805599454, -0.5410441805599453, -0.5403991805599453, -0.5399321805599454, -0.5400171805599453, -0.5413081805599453, -0.5448711805599453, -0.5523081805599452, -0.5657271805599453, -0.5872111805599454, -0.6172408805599453, -0.6519749805599453, -0.6814217805599454, -0.6931471805599453};
//int ju = 0;
// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=dat[i]+log(2.);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// a.setu0(1.0);
// a.setv0(7.0);

// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=dat[i]+log(2.);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// cout << a.freeenergy(field) << endl;

// // for(int i = 0 ; i < 100 ; i++) {
// // field(i,0)=-0.542+5./(1. + exp((90. - i)/1/));
// // }
// // for(int i = 0 ; i < 100 ; i++) {
// // field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// // }

// for(double po = 0. ; po < 90. ; po+=1.0) {
// for(int i = 0 ; i < 100 ; i++) {
// 	double poi =po;
// field(i,0)=0.01+5./(1. + exp((poi - i)/0.1));
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// cout << a.freeenergy(field) << ",";
// }
// cout << endl;
// cout << endl;

// double fg =0.0;

// matrix<double> den = a.densityfromfield(field,fg);

// outfunc(field,"field");
// outfunc(den,"den");


// for(double v = 0.0 ; v < 6.0 ; v+=0.1) {
// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=dat[i]+v+log(2.);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }

// cout << a.freeenergy(field) <<", ";
// }

// cout << endl;
// pausel();

// double dt =1.;

// // cout << field << endl;
// // pausel();
// int runtime = 100000;

// // a.run(field,runtime,0.001);

// // cout << field << endl;

// for(double vfac =0.0 ; vfac < 5.0 ; vfac += 0.1) {


// cout << ju << " " << vfac << endl;
// a.setv0(vfac);



// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=dat[i]+log(2.);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// a.run(field,runtime,dt);

// stringstream ss;
// ss<<ju;
// string endf=ss.str();

// string fi1 = "field_mixed"+endf;
// string de1 = "den_mixed"+endf;
// string fe1 = "freeenergy_mixed"+endf;

// double fe=0.0;
// matrix<double> den = a.densityfromfield(field,fe);
// fe=a.freeenergy(den,field,fe);
// matrix<double> fi(1,1);
// fi(0,0) = fe;

// outfunc(field,fi1);
// outfunc(den,de1);
// outfunc(fi,fe1);

// // as(ju,0) = a.freeenergy(field);


// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=0.01+0.1*rand()/(RAND_MAX);
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// a.run(field,runtime,dt);

// string fi2 = "field_demixed"+endf;
// string de2 = "den_demixed"+endf;
// string fe2 = "freeenergy_demixed"+endf;

// fe=0.0;
// den = a.densityfromfield(field,fe);
// fe=a.freeenergy(den,field,fe);

// fi(0,0) = fe;

// outfunc(field,fi2);
// outfunc(den,de2);
// outfunc(fi,fe2);

// for(int i = 0 ; i < 100 ; i++) {
// field(i,0)=5.;
// }
// for(int i = 0 ; i < 100 ; i++) {
// field(i,1)=-log(exp(-dat[i]) - exp(-field(i,0)));
// }
// a.run(field,runtime,dt);
// string fi3 = "field_homo"+endf;
// string de3 = "den_homo"+endf;
// string fe3 = "freeenergy_homo"+endf;
// // as(ju,2) = a.freeenergy(field);


// fe=0.0;
// den = a.densityfromfield(field,fe);	
// fe=a.freeenergy(den,field,fe);

// fi(0,0) = fe;

// outfunc(field,fi3);
// outfunc(den,de3);
// outfunc(fi,fe3);

// ju++;
// }


//cout << as << endl;
// a.run(field,100000,0.001);

// double nn = 0.0;
// matrix<double> b = a.densityfromfield(field,nn);

// outfunc(b,"den");

//cout << den << endl;


return 0;
}
