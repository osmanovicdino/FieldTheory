/* 
 * File:   basic.h
 * Author: Dino
 *
 * Created on 28 November 2010, 00:14
 */


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>


using namespace std;

/* This file contains many basic functions for error checking and output,
 * as well as some common mathemtical operations */

#ifndef BASIC_H
#define	BASIC_H



const double pii = 3.1415926535897932384626433832795028841971693993751;
const long double pi = 3.1415926535897932384626433832795028841971693993751;
const long double ee=2.7182818284590452353;
const double eee=2.7182818284590452353;

template <class T>
inline void output(char *s, T x ) { cout << s << x; }
template <class T>
inline void input (char *s, T &x) { cout << s; cin >> x; }


//definitions, mainly inline functions, and mathematica plotting tools
#define check cout << "fine up to here" << endl; //basic check
#define mathematica myfile << "\nll=" << //input into mathematica
#define convert(dr) myfile << "\nTable[ll[[i]]={(i-1)*" << dr << ",ll[[i]]},{i,1,Length[ll]}];" << endl;
#define listplot3 myfile << "\nr1=ListPlot[ll,PlotJoined->True,PlotRange->All]" << endl;
#define listplot(a) myfile << "\nr1=ListPlot[ll,PlotJoined->True,PlotRange->All,PlotLabel->"<<a<<"]" << endl;
#define listplot2(a,b) a << "\nr1=ListPlot["<<b<<",PlotJoined->True,PlotRange->All]" << endl;
#define listplot3d(q) q << "\nListPlot3D[ll,PlotRange->All]"<<endl;
#define G(r0,R,t,k) myfile << "\nr2=Plot[G[r," <<r0<<","<<R<<","<<t<<","<<k<<"], {r, 0, "<<R<<"}, PlotRange -> All]" << endl;
#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define SIGN(a,b)((b)>= 0 ? fabs(a) : -fabs(a))
#define ABS(a)( (a>0? (a) :(-a)) )

int inline fasterfloor( const double &x ) { return x > 0 ? (int) x : (int) x - 1; }


inline void error(const char *errmsg)
{
	cerr << "\n Runtime error ... \n";
	cerr << errmsg;
	cerr << "\n...Quitting program! \n";


        double a;
        cin >> a; //pause the program
	exit(1);
}

template <class T>
void declare(ofstream &s, T a, char *ss, bool c = false ) {
    if (c) {
    s.open(ss,ios::app);
    }
    else {
        s.open(ss);
    }
    s.precision(20);
    s << "ll=" << fixed << a << endl;
    listplot3d(s);
}

template <class T>
void mathe(ofstream &s, T a) { //output a templated object T to a file s
    if ( !s.is_open() ) error("open a file first");
    s << "ll=" << fixed << a << endl;
    listplot3d(s);
}



inline void pausel() { //pause the program (this is faster than system("pause")
    cout << "\nenter an integer: " << endl;
    int pause;
    cin >> pause;
}
#define Stop pausel();


template <class T>
void SWAP( T &a, T &b ) { //swap the values of two objects a and b
	T c = a;
	a = b;
	b = c;
}

template <class T>
int sign(T y) { //return -1 if negative and 1 otherwise
	if ( y>0 ) return 1;
	return -1;
}

template <class T>
T absolute(T x) { //retrun the absolute value of a function
    if (x<0) return x*T(-1.0);
    else return x;
}


template <class T, class Q>
T MAX(Q x, T y) { //return the largest of the objects x and y (operator> must be defined)
if ( x>y ) return T(x);
else return y;

}

template <class T, class Q>
T MIN(Q x, T y) { //same as function above but vice versa
if ( x<y ) return T(x);
else return y;

}

inline double arctan(double y, double x) { //arctan(y/x) allowing for singularity
    if ( x == 0 ) return sign(y)*pi;
    if ( x > 0 ) return atan(y/x);
    if( x < 0 && y >= 0 ) return atan(y/x)+pi;
    if( x < 0 && y < 0 ) return atan(y/x)-pi;
            
    else return atan(y/x);
}

inline long double Heaviside(long double x) { //Heaviside step function
    if ( x < 0) return 0;
    else return 1;

}

inline long double arccos(long double x) {
    
    if (x > 1) return 0.;
    if ( x < -1) return pi;
    else return acos(x);
}

inline double KroneckerDelta(const int &i, const int &j) {
    if(i==j) {
        return 1.;
    }
    else{
        return 0.;
    }
}

template <class T, class Q, class F>
struct combin {
        T *func1;
        Q *func2;
        combin(T *funcc1, Q *funcc2 ) : func1(funcc1), func2(funcc2) {}
        F operator()( F x1) {
            return (*func1)(x1)*(*func2)(x1);
        }  
};


template<class H, class T,class Q>
struct combine {
    T &func1;
    Q &func2;
    combine(T &funcc1, Q &funcc2) : func1(funcc1),func2(funcc2) {}
    H operator()(H x, H y) {
        return (func1)(x,y)*(func2)(x,y);
    }
    
    
    
};


    template <class T, class Q, class F>
    struct combinefuncs {
        T *func1;
        Q *func2;
        combinefuncs(T *funcc1, Q *funcc2 ) : func1(funcc1), func2(funcc2) {}
        F operator()( F x1, F x2) {
            return (*func1)(x1,x2)*(*func2)(x1,x2);
        }
    };

    template<class T, class Q>
     struct d1a {
        Q &func;
        T R;
        d1a(Q &funcc, T RR) : func(funcc),R(RR) {}               
        long double operator()(long double x0) {
            return func(x0,R-x0);
                       }
                   };
                   
    template<class T, class Q>
     struct d1a1 {
        Q &func;
        T R;
        d1a1(Q &funcc, T RR) : func(funcc),R(RR) {}               
        long double operator()(long double x0) {
            return (R-x0)*func(x0,R-x0);
                       }
                   };
                   
    template<class T, class Q>
     struct d1a2 {
        Q &func;
        T R,z;
        d1a2(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        long double operator()(long double x0) {
            return (x0-R)*func(x0,z+R-x0);
                       }
                   };
    template<class T, class Q>
     struct d1a3 {
        Q &func;
        T R,z;
        d1a3(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        long double operator()(long double x0) {
            return (x0)*func(x0,z-sqrt(SQR(R)-SQR(x0)))-(x0)*func(x0,z+sqrt(SQR(R)-SQR(x0)));
                       }
                   };
    template<class T, class Q>
     struct d1a4 {
        Q &func;
        T R,z;
        d1a4(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        T operator()(T x0) {
            return ((x0/(sqrt(SQR(R)-SQR(x0))))*func(x0,z-sqrt(SQR(R)-SQR(x0))))+((x0/(sqrt(SQR(R)-SQR(x0))))*func(x0,z+sqrt(SQR(R)-SQR(x0))));
                       }
                   };

       template <class T, class Q>
       struct d1a5 {
           Q &func;
           T R,z;
           d1a5(Q &funcc,T RR, T zz) : func(funcc),R(RR),z(zz) {}
           T operator()(T x0) {
               return -((SQR(x0)/(sqrt(SQR(R)-SQR(x0))))*func(x0,z-sqrt(SQR(R)-SQR(x0))))-((SQR(x0)/(sqrt(SQR(R)-SQR(x0))))*func(x0,z+sqrt(SQR(R)-SQR(x0))));
           }
           
           
       };

template <class T>
T Power(T a, int b) { //number to an integer power
	if (b<0) { cout << "b is: " << b << endl; error("b must be >= 0");  return 0; }
    if (b==0) return 1;
    if (a==0) return 0;
    if (b%2==0) {
        return Power(a*a, b/2);
    } else if (b%2==1) {
        return a*Power(a*a,b/2);
    }
    return 0;
}


	




#endif	/* BASIC_H */




