#ifndef ISINGPOL_H
#define ISINGPOL_H



struct IsingPolymer {

double ds;
double stot;
double b;
double R;
int Nr;
double dr;

vector1<double> radialpoints;
matrix<double> amat;
matrix<double> bmat;
matrix<double> bigD;

int no_types;
matrix<double> V;


vector1<double> propd;

double beta;
double u0;
double v0;


IsingPolymer(int,int);

void setds(double);
void setstot(double);
void setb(double);
void setR(double);
void setV(matrix<double>&);
void setu0(double);
void setv0(double);

double integrate(vector1<double>&);

double determine_stepsize2(matrix<double> &x, matrix<double> &a, double param, double fe1);

double freeenergy(matrix<double> &field);

double freeenergy(matrix<double> &field, matrix<double> &den, double &fe);

vector1<double> radialsolver(vector1<double> &field,double&);

vector1<double> radialsolver2(vector1<double> &field);

matrix<double> densityfromfield(matrix<double>&,double&);

matrix<double> genconvolocal(matrix<double>&);

void iteratefield(matrix<double> &field, double dt, double &fen, double&);
//vector1<double> forcecalc1(double,double,double,double,double,double);
//void forcecalc2(int, int, vector1<double>&);
//vector1<double> forcecalc3(int);

void run(matrix<double> &mf, int tot, double dt);

};

#include "IsingPol.cpp"

#endif