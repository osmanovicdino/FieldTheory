#ifndef ISINGPOL_H
#define ISINGPOL_H



struct IsingPolymer {

double ds;
double stot;
double totpolymers;
double b;
double R;
int Nr;
double dr;
bool polymer_or_fluid;

int cop;

vector1<double> radialpoints;
matrix<double> amat;
matrix<double> bmat;
matrix<double> bigD;
matrix<double> fixed_field;
vector1<bool> *seq;


int no_types;
matrix<double> V;


vector1<double> propd;

double beta;
double u0;
double v0;


IsingPolymer(int,int,int);
~IsingPolymer();

void setds(double);
void setstot(double);
void setb(double);
void setR(double);
void setV(matrix<double>&);
void setu0(double);
void setv0(double);
void setfixed_field(matrix<double>&);
void settotpolymers(double);
void setseq(vector1<bool>&);

template <class T>
T integrate(vector1<T>&);

double determine_stepsize2(matrix<double> &x, matrix<double> &a, double param, double fe1);

double determine_stepsize2_copolymer(matrix<double> &x, matrix<double> &a, double param, double fe1);

double freeenergy(matrix<double> &field);

double freeenergy_copolymer(vector1<double>&, matrix<double>&);

double freeenergy_copolymer(matrix<double>&,bool);

double freeenergy(matrix<double>&, matrix<double>&, double &);

double freeenergy_copolymer(matrix<double>&, matrix<double>&, double&);

double freeenergy_copolymer(matrix<double>&, vector1<double>&, matrix<double>&, double&);

vector1<double> radialsolver(vector1<double>&,double&);

vector1<double> radialsolver2(vector1<double>&);

matrix<double> radialsolver_copolymer(vector1<double>&, vector1<double>&,double&);

matrix<double> densityfromfield(matrix<double>&);

matrix<double> densityfromfield_copolymer(matrix<double>&);

matrix<double> densityfromfield(matrix<double>&,double&);


matrix<double> genconvolocal(matrix<double>&);

matrix<double> genconvolocal_copolymer(matrix<double>&);

void iteratefield(matrix<double> &field, double dt, double &fen, double&);

void iteratefield_copolymer(matrix<double> &field, double dt, double &fen, double&);
//vector1<double> forcecalc1(double,double,double,double,double,double);
//void forcecalc2(int, int, vector1<double>&);
//vector1<double> forcecalc3(int);

void run(matrix<double> &mf, int tot, double dt);
void run_copolymer(matrix<double> &mf, int tot, double dt);

matrix<double> densityfromfield_fluid(matrix<double>&);
matrix<double> densityfromfield_fluid(matrix<double>&,double&);

matrix<double> densityfromfield_copolymer(vector1<double>&,matrix<double>&,double&);


};

#include "IsingPol.cpp"

#endif