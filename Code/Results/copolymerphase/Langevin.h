#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "MD.h"

class LangevinNVT : public MD {
protected:
	double gamma;
	double dt;
	double kT;
	double m;
	matrix<double> *mom; //momenta of all the particles;

	double d;//((gamma*dt/2.));
	double q;//((dt)/2.);
	double r;//(sqrt(0.5*kT*(gamma)*(m)*(dt)));

	double c1;// = (dt/m);
	double c2;// = (1.0/(1.0+(d)));
	double c3;// = (1.0/(1.0+(d)))*q;
	double c4;// = (1.0/(1.0+(d)))*r;
	double c5;// = (1-(d));	
public:
	LangevinNVT();
	LangevinNVT(geometry &a); //define a Langevin system with a given geometry
	void setgamma(double); //set the damping factor
	void setdt(double); //set the timestep
	void setkT(double); //set the temperature
	void setm(double); //set the mass of the particles
	void setmom(matrix<double>&); //set the momentum of the particles

	double getkT() {return kT;}

	matrix<double>& getmom(); //return a matrix of the particles momenta
	vector1<double> avmom();
	double getdt() { return dt;}
	void printparams();

	void advance_pos(); //advance the position of the particles
	virtual void advance_mom(matrix<double>&,matrix<double>&);//advance the momenta a full step with the arguments being the deterministic and random forces

	virtual void initialize(matrix<int>&);
	virtual void adv(matrix<int>&);

	//virtual void advm(matrix<int>&,vector1<int>&,vector1<bool>&,intmatrix&); //advance with forces that depend on the internal state of the particle
};

#include "Langevin.cpp"

#endif