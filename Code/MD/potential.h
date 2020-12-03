#ifndef POTENTIAL_H
#define POTENTIAL_H

struct potential { //this is the base constructor for all potentials.
	bool dl; //does this interaction decay at large distances?
	double interaction_distance;

	virtual potential* clone() const = 0;

	virtual double energy(double)=0;
	virtual double force(double)=0;
	virtual double force2mdx(double)=0; //the factor as a function of distance*squared which when multiplied with dx will give the total force

	virtual void setparameters(const vector1<double>&)=0;
	virtual void printparameters()=0;
	virtual void printparameters(ofstream&)=0;


};

struct LJPotential : potential { //Lennard Jones potential
	double epsilon;
	double sigma;
	double sigma6,sigma12;
	LJPotential(double a, double b) : epsilon(a),sigma(b) { dl=true; interaction_distance=2.5*sigma; sigma6= SQR(CUB(sigma)); sigma12 = SQR(sigma6);}
	double energy(double x) {
        if(x<interaction_distance) {
        double a;
        a=x*x*x;
        
        a=SQR(a);


        return 4*epsilon*(sigma12/SQR(a)-sigma6/a);	
        }
        else return 0;	
	}
	double force(double x) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        double a;
        a=x*x*x;
        
        a=SQR(a);

		return 4*epsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));
	}
	double force2mdx(double x2) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        double a;
        // a=x*x*x;
        
        a=x2*x2*x2;

		return 4*epsilon*(12*sigma12/(x2*SQR(a))-6*sigma6/(x2*a));
	}	

	void setparameters(const vector1<double> &x) {
		if(x.getsize() != 2) error("parameters vector for LJ not of the right size");
		epsilon = x.gpcons(0);
		sigma = x.gpcons(1);
		dl=true;
		interaction_distance=2.5*sigma;
		sigma6= SQR(CUB(sigma));
		sigma12 = SQR(sigma6);
	}

	void printparameters() {
		cout << "LJ epsilon: " << epsilon << " sigma: " << sigma << endl; 
	}
	void printparameters(ofstream &s) {
		s << epsilon << "," << sigma; 
	}	

	LJPotential* clone() const {
		return new LJPotential(*this);
	}
};

struct HSPotential: potential {
	double paramlimit;
	double epsilon;
	double sigma;
	double sigma6,sigma12;	
	HSPotential(double a, double b) : epsilon(a), sigma(b) {paramlimit=1.12246*sigma; dl=true; interaction_distance=paramlimit+sigma; sigma6= SQR(CUB(sigma)); sigma12 = SQR(sigma6); }

	double energy(double x) {
      	if(x < paramlimit) {
        	double a;
        	a=x*x*x;
       		a=SQR(a);
        return 4*epsilon*(sigma12/SQR(a)-sigma6/a)+epsilon;	
        }
        else {
        	return 0;
        }
	}
	double force(double x) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        if(x < paramlimit) {
        double a;
        a=x*x*x;
        
        a=SQR(a);

		return 4*epsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));
		}
		else {
			return 0;
		}
	}
	double force2mdx(double x2) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        if(x2 < SQR(paramlimit)) {
        double a;
        a=x2*x2*x2;
        
        //a=SQR(a);

		return 4*epsilon*(12*sigma12/(x2*SQR(a))-6*sigma6/(x2*a));
		}
		else {
			return 0;
		}
	}	

	void setparameters(const vector1<double> &x) {
		if(x.getsize() != 2) error("parameters vector for LJ not of the right size");
		epsilon = x.gpcons(0);
		sigma = x.gpcons(1);
		dl=true;
		interaction_distance=2.5*sigma;
		sigma6= SQR(CUB(sigma));
		sigma12 = SQR(sigma6);
	}

	void printparameters() {
		cout << "HS epsilon: " << epsilon << " sigma: " << sigma << endl; 
	}
	void printparameters(ofstream &s) {
		s << epsilon << "," << sigma; 
	}

	HSPotential* clone() const {
		return new HSPotential(*this);
	}

};

struct WCAPotential : potential {
	double epsilon; //epsilon
	double sigma; //sigma
	double attrepsilon; //epsilon modifier
	double paramlimit;
	double sigma6;
	double sigma12;

	WCAPotential(double a, double b, double c) : epsilon(a),sigma(b),attrepsilon(c) {
		sigma6= SQR(CUB(sigma)); sigma12 = SQR(sigma6);
		interaction_distance=2.5*sigma;
		dl = true;
		paramlimit=1.12246*sigma;
	}

	double energy(double x) {
      	if(x < paramlimit) {
        	double a;
        	a=x*x*x;
       		a=SQR(a);
        return 4*epsilon*(sigma12/SQR(a)-sigma6/a)+epsilon-attrepsilon;	
        }
        else if(x >= paramlimit && x< interaction_distance) {
            double a;
        	a=x*x*x;
       		a=SQR(a);
       		return 4*attrepsilon*(sigma12/SQR(a)-sigma6/a);	
        }
        else {
        	return 0;
        }
	}

	double force(double x) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        if(x < paramlimit) {
        double a;
        a=x*x*x;
        
        a=SQR(a);

		return 4*epsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));
		}
		else if(x >= paramlimit && x< interaction_distance) {
        double a;
        a=x*x*x;
        
        a=SQR(a);

		return 4*attrepsilon*(12*sigma12/(x*SQR(a))-6*sigma6/(x*a));			
		}
		else {
			return 0;
		}
	}
	double force2mdx(double x2) {
		//where x is the distance, this is the magnitude of the force, it is multiplied by r/|r| to give vectorial
        if(x2 < SQR(paramlimit)) {
        double a;
        a=x2*x2*x2;
        
        //a=SQR(a);

		return 4*epsilon*(12*sigma12/(x2*SQR(a))-6*sigma6/(x2*a));
		}
		else if(x2 >= SQR(paramlimit) && x2< SQR(interaction_distance)) {
        double a;
        a=x2*x2*x2;
        
        //a=SQR(a);

		return 4*attrepsilon*(12*sigma12/(x2*SQR(a))-6*sigma6/(x2*a));			
		}
		else {
			return 0;
		}
	}	

	void setparameters(const vector1<double> &x) {
		if(x.getsize() != 3) error("parameters vector for LJ not of the right size");
		epsilon = x.gpcons(0);
		sigma = x.gpcons(1);
		attrepsilon = x.gpcons(2);
		dl=true;
		interaction_distance=2.5*sigma;
		sigma6= SQR(CUB(sigma));
		sigma12 = SQR(sigma6);
	}

	void printparameters() {
		cout << "WCA epsilon: " << epsilon << " sigma: " << sigma << " att epsilon: " << attrepsilon <<endl; 
	}
	void printparameters(ofstream &s) {
		s << epsilon << "," << sigma << "," << attrepsilon; 
	}

	WCAPotential* clone() const {
		return new WCAPotential(*this);
	}		
};

struct HarmonicPotential : potential {
	double k;
	double x0;

	HarmonicPotential(double a, double b) : k(a),x0(b) { dl=false; interaction_distance = 0.0;}

	double energy(double x) {
		return 0.5*k*SQR(x-x0);
	}
	double force(double x) {
		return -k*(x-x0);
	}
	double force2mdx(double x2) {
		return -k;
	}

	void setparameters(const vector1<double> &x) {
		if(x.getsize() != 2) error("parameters vector for HarmonicPotential not of the right size");
		this->k = x.gpcons(0);
		this->x0 = x.gpcons(1);
	}
	void printparameters() {
		cout << "harmonic spring constant: " << k << " minima: " << x0 << endl; 
	}
	void printparameters(ofstream &s) {
		s << k << "," << x0; 
	}	

	HarmonicPotential* clone() const {
		return new HarmonicPotential(*this);
	}

};

struct FENEPotential : potential {
	double kbond;
	double R_0;

	FENEPotential(double a, double b) : kbond(a),R_0(b) { dl = false; interaction_distance =0.0; }

	double energy(double x) { 
		if(x<R_0) return -(kbond/2.)*SQR(R_0)*log(1.0-SQR(x)/SQR(R_0));
		else return 1.E5*x;
	}

	double force(double x) { return -kbond*x/(1-SQR(x)/SQR(R_0));}

	double force2mdx(double x2) { return -kbond/(1-x2/SQR(R_0));}

	void setparameters(const vector1<double> &x) {
		if(x.getsize() != 2) error("parameters vector for FENEPotential not of the right size");
		this->kbond =x.gpcons(0);
		this->R_0 = x.gpcons(1);
	}
	void printparameters() {
		cout << "FENE spring constant: " << kbond << " R_0: " << R_0 << endl; 
	}
	void printparameters(ofstream &s) {
		s << kbond << "," << R_0; 
	}	

	FENEPotential* clone() const {
		return new FENEPotential(*this);
	}

};

struct ColoumbPotential : potential {
	double k;

	ColoumbPotential(double a) : k(a) {dl = true; interaction_distance = 50.0*k;}

	double energy(double x) { return k/x; }

	double force(double x) { return -k/SQR(x); }

	double force2mdx(double x2) { return -k/(x2*sqrt(x2)); }

	void setparameters(const vector1<double> &x) {
	if(x.getsize() != 1) error("parameters vector for ColoumbPotential not of the right size");
	this->k = x.gpcons(0);
	}
	void printparameters() {
		cout << "electrostatic constant: " << k << endl; 
	}
	void printparameters(ofstream &s) {
		s << k; 
	}		
	ColoumbPotential* clone() const {
		return new ColoumbPotential(*this);
	}
};

// struct ScreenedColoumbPotential : potential {
// 	double k;
// 	double la;

// 	ScreenedColoumbPotential(double a, double b): k(a),la(b) {dl = true; interaction_distance = 2*la;}

// 	double energy(double x) { return (k/x)*(exp(-x/la)); }

// 	double force(double x) {  return -(k*(exp(-x/la)/SQR(x))) - k*exp(-x/la)/la*x; }

// 	void setparameters(const vector1<double> &x) {
// 	if(x.getsize() != 2) error("parameters vector for ScreenedColoumbPotential not of the right size");
// 	this->k = x.gpcons(0);
// 	this->la = x.gpcons(1);
// 	}
// 	void printparameters() {
// 		cout << "spring constant: " << k << " screening length: " << la << endl; 
// 	}
// 	void printparameters(ofstream &s) {
// 		s << k << "," << la; 
// 	}
// 	ScreenedColoumbPotential* clone() const {
// 		return new ScreenedColoumbPotential(*this);
// 	}			
// };

// struct ExponentialPotential : potential {
// 	double k;
// 	double la;

// 	ExponentialPotential(double a, double b): k(a),la(b) {dl = true; interaction_distance = 2*la;}

// 	double energy(double x) { return k*(exp(-x/la)); }

// 	double force(double x) {  return k*exp(-x/la)/la; }

// 	void setparameters(const vector1<double> &x) {
// 	if(x.getsize() != 2) error("parameters vector for ExponentialPotential not of the right size");
// 	this->k = x.gpcons(0);
// 	this->la = x.gpcons(1);
// 	}
// 	void printparameters() {
// 		cout << "spring constant: " << k << " screening length: " << la << endl; 
// 	}
// 	void printparameters(ofstream &s) {
// 		s << k << "," << la; 
// 	}
// 	ExponentialPotential* clone() const {
// 		return new ExponentialPotential(*this);
// 	}			
// };

#endif