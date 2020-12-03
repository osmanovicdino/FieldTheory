double sphericalbessel0(double r) {
	if(r<1E-10) return 1.0;
	else return sin(r)/r;
}

IsingPolymer::IsingPolymer(int Nrr, int no_typess) : radialpoints(vector1<double>(Nrr)) , propd(vector1<double>(Nrr)), amat(matrix<double>(Nrr,Nrr)),bmat(matrix<double>(Nrr,Nrr)),bigD(matrix<double>(Nrr,Nrr)),V(matrix<double>(no_typess,no_typess)) {
	ds = 1.0;
	stot =  1000.0;
	b = 1.0;
	R = 10.0;
	Nr = Nrr;
	no_types = no_typess;


	dr =  R/((double)(Nr-1));
	for(int i = 0 ; i < Nr ; i++) {
	radialpoints[i]=i*dr;
	}

	for(int i = 0  ; i < Nr ; i++) {
		propd[i]=exp(-SQR(b)*SQR(i)*SQR(pi)*ds/(6.*SQR(R)));
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
			amat(i,j)=(2./(double)Nr)*SQR(j*pi/R)*sphericalbessel0(radialpoints[i]*j*pi/R)*SQR(radialpoints[i]);
			bmat(i,j)=sphericalbessel0(radialpoints[i]*j*pi/R);
		}
	}

	matrix<double> tempd(Nr,Nr);
	for(int i = 0 ; i < Nr ; i++) {
		tempd(i,i)=propd[i];
	}
	tempd =  tempd*amat;
	bigD = bmat*tempd;

	//bigD = bmat*(tempd*amat);

	beta = 1.0;
	u0 = 1.0;
	v0 = 1.0;
}

void IsingPolymer::setds(double dss) {
	ds = dss;
	for(int i = 0  ; i < Nr ; i++) {
		propd[i]=exp(-SQR(b)*SQR(i)*SQR(pi)*ds/(6.*SQR(R)));
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
			amat(i,j)=(2./(double)Nr)*SQR(j*pi/R)*sphericalbessel0(radialpoints[i]*j*pi/R)*SQR(radialpoints[i]);
			bmat(i,j)=sphericalbessel0(radialpoints[i]*j*pi/R);
		}
	}	
	matrix<double> tempd(Nr,Nr);
	for(int i = 0 ; i < Nr ; i++) {
		tempd(i,i)=propd[i];
	}
	tempd =  tempd*amat;
	bigD = bmat*tempd;



};
void IsingPolymer::setstot(double stott) {
	stot = stott;

}
void IsingPolymer::setb(double bb) {
	b = bb;
	for(int i = 0  ; i < Nr ; i++) {
		propd[i]=exp(-SQR(b)*SQR(i)*SQR(pi)*ds/(6.*SQR(R)));
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
			amat(i,j)=(2./(double)Nr)*SQR(j*pi/R)*sphericalbessel0(radialpoints[i]*j*pi/R)*SQR(radialpoints[i]);
			bmat(i,j)=sphericalbessel0(radialpoints[i]*j*pi/R);
		}
	}	
	matrix<double> tempd(Nr,Nr);
	for(int i = 0 ; i < Nr ; i++) {
		tempd(i,i)=propd[i];
	}
	tempd =  tempd*amat;
	bigD = bmat*tempd;

}
void IsingPolymer::setR(double RR) {
	R = RR;
	dr =  R/((double)(Nr-1));
	for(int i = 0 ; i < Nr ; i++) {
	radialpoints[i]=i*dr;
	}

	for(int i = 0  ; i < Nr ; i++) {
		propd[i]=exp(-SQR(b)*SQR(i)*SQR(pi)*ds/(6.*SQR(R)));
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
			amat(i,j)=(2./(double)Nr)*SQR(j*pi/R)*sphericalbessel0(radialpoints[i]*j*pi/R)*SQR(radialpoints[i]);
			bmat(i,j)=sphericalbessel0(radialpoints[i]*j*pi/R);
		}
	}
	matrix<double> tempd(Nr,Nr);
	for(int i = 0 ; i < Nr ; i++) {
		tempd(i,i)=propd[i];
	}
	tempd =  tempd*amat;
	bigD = bmat*tempd;


}
void IsingPolymer::setV(matrix<double> &V2) {
	if(V2.getnrows() != no_types || V2.getncols() != no_types) error("incorrect size for matrix of interactions");
	for(int i = 0  ; i < no_types ; i++) {
		for(int j = 0  ; j < no_types ; j++) {
			V(i,j)=V2(i,j);
		}
	}

}
void IsingPolymer::setu0(double u00) {
	u0 = u00;
}
void IsingPolymer::setv0(double v00) {
	v0 = v00;
}

double IsingPolymer::integrate(vector1<double> &m) {
	double tot = 0.0;
	for(int i  = 0 ; i < Nr ; i++)
	tot+=m[i]*4*pi*SQR(radialpoints[i])*dr;

	return tot;
}


vector1<double> IsingPolymer::radialsolver(vector1<double> &field, double &norm) {
	vector1<double> den(Nr);
	// vector1<double> propw(Nr);

	// for(int i  = 0 ; i < Nr ; i++) {
	// 	propw[i]=exp(-field[i]*ds/2.);
	// }

	matrix<double> propw(Nr,Nr);
		for(int i  = 0 ; i < Nr ; i++) {
		propw(i,i)=exp(-field[i]*ds/2.);
	}
	matrix<long double> itmat = (propw)*(bigD*propw);

	int totsize = 1+(int)(stot/ds);

	matrix<long double> store(totsize,Nr);
	for(int j  = 0 ; j < Nr-1 ; j++) {
		store(0,j)=1.0;
	}
	store(0,Nr-1)=0.0;

	vector1<long double> q13(Nr);
	for(int j = 0  ; j < Nr ; j++) {
		q13[j] = store(0,j);
	}
	for(int i = 0  ; i < totsize-1 ; i++) {

		// q13=propw&q13;
		// q13=amat*q13;
		// q13=propd&q13;
		// q13=bmat*q13;
		// q13=propw&q13;

		q13=itmat*q13;

		for(int j = 0  ; j < Nr ; j++) {
			store(i+1,j)=q13[j];
		}
		
	}

	//cout << store(100,'r') << endl;
	//vector1<double> den(Nr);
	vector1<long double> den2(Nr);
	for(int i = 0 ; i < totsize ; i++) {
		for(int j = 0 ; j < Nr ; j++) {		
		den2[j]+=(store(i,j)*store(totsize-1-i,j));
		}
	}
	long double tot = 0.0;
	for(int i  = 0 ; i < Nr ; i++)
	tot+=den2[i]*4.*pi*SQR(radialpoints[i])*dr;

	//cout << "tot: " << tot << endl;	

	double V=(4./3.)*pi*R*R*R;
	// double qint = 0.0;
	// for(int i  = 0 ; i < Nr ; i++)
	// qint+=store(totsize-1,i)*4*pi*SQR(radialpoints[i])*dr;

	
	norm=-log((tot/(long double)stot)/(long double)V);

	//cout << norm << endl;

	den2*=((long double)stot/tot);
	for(int i = 0  ; i < Nr ; i++) {
		den[i]=den2[i];
	}

	//return stot*(den/tot);
	return den;

}

vector1<double> IsingPolymer::radialsolver2(vector1<double> &field) {
	vector1<double> den(Nr);
	vector1<double> propw(Nr);

	for(int i  = 0 ; i < Nr ; i++) {
		propw[i]=exp(-field[i]*ds/2.);
	}

	// matrix<double> propw(Nr,Nr);
	// 	for(int i  = 0 ; i < Nr ; i++) {
	// 	propw(i,i)=exp(-field[i]*ds/2.);
	// }
	// matrix<double> itmat = (propw)*(bigD*propw);

	int totsize = (int)(stot/ds);

	matrix<double> store(totsize,Nr);
	for(int j  = 0 ; j < Nr-1 ; j++) {
		store(0,j)=1.0;
	}
	store(0,Nr-1)=0.0;

	vector1<double> q13(Nr);
	for(int j = 0  ; j < Nr ; j++) {
		q13[j] = store(0,j);
	}
	for(int i = 0  ; i < totsize-1 ; i++) {

		q13=propw&q13;
		q13=amat*q13;
		q13=propd&q13;
		q13=bmat*q13;
		q13=propw&q13;

		//q13*=itmat;
		
		for(int j = 0  ; j < Nr ; j++) {
			store(i+1,j)=q13[j];
		}
		
	}

	//vector1<double> den(Nr);

	for(int i = 0 ; i < totsize ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
		den[j]+=(store(i,j)*store(totsize-1-i,j));
		}
	}
	double tot = 0.0;
	for(int i  = 0 ; i < Nr ; i++)
	tot+=den[i]*4*pi*SQR(radialpoints[i])*dr;


	return stot*(den/tot);
}

matrix<double> IsingPolymer::densityfromfield(matrix<double> &field, double &fi) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field(i,j));
		}
		totalfield[i]=-(1/beta)*log(norm[i]);
	}
	//cout << totalfield << endl; 

	vector1<double> density = this->radialsolver(totalfield,fi);

	// cout << fi << endl;
	// vector1<double> fac = density&totalfield;
	// cout << fi-(this->integrate(fac)) << endl;

	matrix<double> partialden(Nr,no_types);
	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density[i]*exp(-beta*field(i,j))/norm[i];		
		}
	}

	// vector1<double> fro = partialden(0,'c');
	// vector1<double> fro1 = partialden(1,'c');

	// vector1<double> gro = field(0,'c');
	// vector1<double> gro1 = field(1,'c');

	// vector1<double> temp = fro&gro;
	// vector1<double> temp2 = fro1&gro1;

	// cout << this->integrate(density) << endl;
	// cout << this->integrate(fro) << endl;
	// cout << this->integrate(fro1) << endl;
	// cout << this->integrate(temp) << endl;
	// cout << this->integrate(temp2) << endl;
	// cout << "horse" << endl;
	// cout << endl;
	return partialden;
}

matrix<double> IsingPolymer::genconvolocal(matrix<double> &density) {
	matrix<double> convo(Nr,no_types);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < no_types ; j++) {
		
			double fac1 = 0.0;
			double fac2 = 0.0;

			for(int k1  = 0  ; k1 < no_types ; k1++)
				fac1+=density(i,k1);
			for(int k2 = 0 ; k2 < no_types ; k2++)
				fac2+=V(j,k2)*density(i,k2);


		convo(i,j)=u0*fac1+v0*fac2;
		}
	}
	return convo;
}

double IsingPolymer::freeenergy(matrix<double> &field) {
double fe=0.0;
//cout << field << endl;
matrix<double> den = densityfromfield(field,fe);
//cout << den << endl;

double tot =0.0;

for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types ; k2++)
for(int j = 0  ; j <no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3;
}


double IsingPolymer::freeenergy(matrix<double> &den, matrix<double> &field, double &fe) {

double tot =0.0;

for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;

for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types ; k2++)
for(int j = 0  ; j <no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
//cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;
//cout << tot << endl;
// cout << field << endl;

// cout << tot2 << endl;
// cout << tot3 << endl;

return fe-tot+u0*tot2+v0*tot3;

}

double IsingPolymer::determine_stepsize2(matrix<double> &x, matrix<double> &a, double param, double fe1) {

	// x is our initial mean field;
	// a is our convoled meanfield
	//fe1 is our initial free energy
	//param is some parameter that determines step size
    double stepsize = 0.45*param;        
    
    matrix<double> temp2 = x+stepsize*(a-x);
    
    double fe2 = this->freeenergy(temp2);
    double param2;
    if(fe2-fe1<-1E-5) {
        double steptemp = stepsize*2;

        matrix<double> gemp2 = x+steptemp*(a-x);                
        double ftemp = this->freeenergy(gemp2);
        
        if(ftemp < fe2) param2 = steptemp;
        else {
            double a = -((-(fe1*stepsize) + ftemp*stepsize + fe1*steptemp - fe2*steptemp)/
 (stepsize*(stepsize - steptemp)*steptemp));
            double b = -((fe1*SQR(stepsize) - ftemp*SQR(stepsize) - 
   fe1*SQR(steptemp) + fe2*SQR(steptemp))/
 (stepsize*(stepsize - steptemp)*steptemp));
            param2 = -b/(2.*a);                
        }
//            cout << fe1 << "," << fe2 << "," << ftemp << endl;
//            cout << 0 << "," << stepsize << "," << steptemp << endl;
        
    }
    else if (fe1-fe2 < -1E-5) {
        int iter = 0;
        bool lesser = true;
        double steptemp = stepsize;
        double ftemp;
        while(fe1 < fe2) {
            if(iter > 6) {lesser = false; break; }
            ftemp = fe2;
            steptemp = stepsize;
            stepsize/=2;
            matrix<double> gemp2 = x+stepsize*(a-x);
            fe2 = this->freeenergy(gemp2);

        }
        
        if(lesser) {
            double a = -((-(fe1*stepsize) + ftemp*stepsize + fe1*steptemp - fe2*steptemp)/
 (stepsize*(stepsize - steptemp)*steptemp));
            double b = -((fe1*SQR(stepsize) - ftemp*SQR(stepsize) - 
   fe1*SQR(steptemp) + fe2*SQR(steptemp))/
 (stepsize*(stepsize - steptemp)*steptemp));
            param2 = -b/(2.*a);
        }
        else {
            param2 = stepsize/2.;
        }
//            cout << fe1 << "," << fe2 << "," << ftemp << endl;
//            cout << 0 << "," << stepsize << "," << steptemp << endl;
        // have fe1, fe2, ftemp
        // with 0  , stepsize, steptemp
        
    }
    else {
        cout << "fe1=fe2" << endl;
        param2 = stepsize;
        
    }

    return param2;        
}    

void IsingPolymer::iteratefield(matrix<double> &field, double dt, double &fen, double &stepsize2) {
	double fe=0.0;
	matrix<double> den = this->densityfromfield(field, fe);
	// outfunc(den,"den");
	// outfunc(field,"field");
	double fen1  = this->freeenergy(den,field,fe);
	cout << "free energy: " << fen1 << endl;
	 
	fen = fen1;
	matrix<double> convolocal = this->genconvolocal(den);
	double stepsize = this->determine_stepsize2(field,convolocal,dt,fen1);
	cout << "steps: " << dt << " " << stepsize << endl;
	stepsize2 = stepsize;
	field = field+stepsize*(-field+convolocal);
}

void IsingPolymer::run(matrix<double> &field, int tot, double dt) {
	//cout << field << endl;
	double fen = 0.0;
	for(int i = 0 ; i < tot ; i++)
	{
		double fi = fen;
		cout << i << endl;
		double stepsize2 = 0.0;
		this->iteratefield(field,dt,fen,stepsize2);
		double ff = fen;
		//cout << abs(fi-ff) << endl;
		if(abs(fi-ff)<0.1*stepsize2) { cout << "convergence test made" << endl; break; }
	}
	//cout << field << endl;
}