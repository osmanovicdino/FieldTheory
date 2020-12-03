#ifndef MIPS_H
#define MIPS_H

#include "../DataStructures/vector1.h"

#include "../DataStructures/matrix2.h"

typedef vector1<matrix<double>> vecmat;

struct mips_mf {
int Nr;
int Nq;

int N;

vector1<matrix<double> > dat; //where our data is stored


vector1<matrix<COMPLEX> > potm;

double l;
vector1<double> xpos;
vector1<double> ypos;
vector1<double> qpos;

double sigma;
double epsilon;
double v0;

// vector1<matrix<double> > gen

double dx;
double dq;



mips_mf(int Nr, int Nq, int, double l, double vo);
inline double potential(double, double, double, double, double, double);

void set_eps(double a) { epsilon = a; }
void set_sig(double a) { sigma = a; }
void set_v0(double a) { v0 = a; }

double min(vector1<matrix<double> > &);
double max(vector1<matrix<double> > &);
// void gen_fm();
matrix<double> qintegral();

double integrate();

double integrate(vecmat&);

matrix<double> calc_HS();

vector1< matrix<double> > calc_int();

vector1<matrix<double>> calc_int2();

void updatedat(double);
void updatedatstore();
};

#include "mips.cpp"

#endif /* MIPS_H */
