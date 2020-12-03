#ifndef ISINGPOLPARTIAL_H
#define ISINGPOLPARTIAL_H



struct IsingPolymerPartial : IsingPolymer {

IsingPolymerPartial(int,int,int);

vector1<double> get_other_field(vector1<double> &field1,vector1<double> &fixed);

matrix<double> density_fixedden(vector1<double> &field, vector1<double> &field1, double &fe);
void determinesequence(matrix<double> &density, int s );
matrix<double> density_fixedseq(matrix<double> &field, double &fe);


double freeenergy_fixedseq(matrix<double> &field);
double freeenergy_fixedseq(matrix<double> &den, matrix<double> &field, double &fe);
matrix<double> genconvolocal_fixedseq(matrix<double> &density) ;
void iteratefield_fixedseq(matrix<double> &field, double dt, double &fen, double &stepsize2);
double determine_stepsize2_fixedseq(matrix<double> &x, matrix<double> &a, double param, double fe1); //determine stepsize for annealed copolymer with fixed bookmarks


// double freeenergy_fixedden(vector1<double>&,vector1<double>&);
// vector1<double> genconvolocal_fixedden(matrix<double> &density) ;
// void iteratefield_fixedden(vector1<double> &field, vector1<double> &fixed, double dt, double &fen, double &stepsize2);
// double determine_stepsize2_fixedden(vector1<double> &x, vector1<double> &a, vector1<double> &fixed, double param, double fe1); //determine stepsize for annealed copolymer with fixed bookmarks
double freeenergy_fixedden(vector1<double> &density,vector1<double> &field,vector1<double> &prop,double &fe);
vector1<double> genconvolocal_fixedden(vector1<double> &density,vector1<double> &field,vector1<double> &prop);
void iteratefield_fixedden(vector1<double> &density, vector1<double> &field,vector1<double> &prop, double dt, double &fen, double fe, double &stepsize2);

void corrfunc(vector1<double>&);

void run_fixedseq(matrix<double> &mf, int tot, double dt);
void run_fixedden(vector1<double> &mf, vector1<double> &den, int tot, double dt);

void run_mixed();

};


#include "IsingPolPartial.cpp"
#endif 
