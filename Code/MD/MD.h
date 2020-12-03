#ifndef MD_H
#define MD_H

#include "../DataStructures/vector1.h"
#include "../DataStructures/vector1.cpp"

#include "../DataStructures/matrix2.h"
#include "../DataStructures/matrix2.cpp"
//#include "matrix.cpp"

// #include "intmatrix.h"
#include "geometry.h"
#include "potential.h"

//#include "vec_vec.h"

class MD { //virtual base class of all calculations 3 D
protected:

geometry *geo; //specify the geometry of the system

matrix<double> *dat; //specifies the positions of all the particles this is the baseline of all molecular dynamics

potential *ints; //interactions between all the elements in the system

int dimension;

public:


MD();
MD(const MD&);
~MD();

void setdat(matrix<double>&);

void setinteractions(potential&);
void setgeometry(geometry&);

double getcoordinate(int,int);

int getN(); // get the number of particles
int getdimension();
matrix<double>& getdat();
geometry& getgeo();
potential& getints();




// double distance(const int&,const int&); //the distance between particles i and j; 
// bool distance_less_than(const int&,const int&, double);
// inline vector1<double> unit_vector(const int&, const int&, double&);
// inline void unit_vector(const int&, const int&, vector1<double>&, double&);
// inline double getcoordinate(int,int);
// void setcoordinate(int,int,double);
//complicated interactions
//void setinteractions();
// double distance(int&, int&); //the distance between particles i and j; 
// bool distance_less_than(const int&,const int&, double);



matrix<int>* calculatepairs(matrix<int>&,double); //calculate the pairs over which forces exist

matrix<int>* calculatepairs_sorted(matrix<int>&,double); //calculate the pairs over which forces exist

//virtual matrix<int>* calculatepairs(matrix<int>&);//calculate pairs using the cell list algorithm

//virtual matrix<int>* calculatepairs_global(matrix<int>&,double);//calculate pairs using the cell list algorithm with a global cut off distance

//virtual matrix<int>* calculatepairs_global_onetype(matrix<int>&,double);//calculate pairs using the cell list algorithm with a global cut off distance

//s_matrix<double> calculatepaircorrelation(double,double);

//s_matrix<double> calculateforces_truncated(matrix<int>&,double); //calculate the forces using the pairs as an input, with the forces truncated

//s_matrix<double> calculateforces(matrix<int>&); //calculate the forces using the pairs as an input

matrix<double> calculateforces(matrix<int>&); //calculate the forces using the pairs as an input
matrix<double> calculatestress(matrix<int>&); //calculate the stress;
matrix<double> calculateforces_truncated(matrix<int>&, double); //calculate the forces using the pairs as an input
matrix<double> calculateforces_old(matrix<int>&);
//s_matrix<double> calculateforces(matrix<int>&,vector1<int>&,vector1<bool>&,intmatrix&); //calculate the forces using the pairs as an input

//virtual double calculate_energy(vec_vec<int>&);

//virtual double calculate_energy(matrix<int>&);
//virtual void adv(vec_vec<int>&); //advance the system one timestep with the pairs as an input

virtual void adv(matrix<int>&)=0; //advance the system one timestep with the pairs as an input

};

struct mdpair {
	int a;
	int b;
	mdpair() : a(0), b(0) {}
	mdpair(const int &aa,const int &bb) : a(aa), b(bb) {}
	bool operator<(mdpair const &rhs) const {
		return a < rhs.a || (a==rhs.a && b < rhs.b);
	} 
};

#include "MD.cpp"

#endif