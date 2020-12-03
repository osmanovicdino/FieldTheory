IsingPolymerPartial::IsingPolymerPartial(int Nrr, int no_typess, int copp=0) : IsingPolymer(Nrr,no_typess,copp) {




}

vector1<double> IsingPolymerPartial::get_other_field(vector1<double> &field1,vector1<double> &fixed) {
	vector1<double> field2(Nr);
	for(int i = 0  ; i < Nr ; i++) {
		field2[i] = -(1./beta)*log(exp(-beta*fixed[i])-exp(-beta*field1[i]));
	}
	return field2;
}

// matrix<double> IsingPolymerPartial::density_fixedden(vector1<double> &field, vector1<double> &field1, double &fe) { // field 1 and field 2 equals field
// 	double fi = 0.0;
// 	vector1<double> density = this->radialsolver(field,fi);
// 	fe= fi;
// 	//vector1<double> field2(Nr);
// 	vector1<double> norm(Nr);

// 	vector1<double> field2 = get_other_field(field1,field);


// 	for(int i = 0  ; i < Nr ; i++) {
		
// 			norm[i]=exp(-beta*field1[i])+exp(-beta*field2[i]);
		
// 		//totalfield[i]=-(1./beta)*log(norm[i]);
// 	}


// 	matrix<double> partialden(Nr,no_types);
// 	for(int i = 0  ; i < Nr ; i++) {
		
// 			partialden(i,0)=density[i]*exp(-beta*field1[i])/norm[i];	
// 			partialden(i,1)=density[i]*exp(-beta*field2[i])/norm[i];		
		
// 	}
// 	return partialden;
// }

void IsingPolymerPartial::determinesequence(matrix<double> &den, int s) { // s is the centre of the sequence
	vector1<double> denA(Nr);

	for(int i = 0 ; i < Nr ; i++) {
		denA[i]=den(i,0);
	}

	double totA = integrate(denA);

	int ch = (int)totA/4.;
	int totsize = 1+(int)(stot/ds);
	vector1<bool> seq(totsize,false);

	for(int i = s-ch  ; i < s + ch ; i++) {
		seq[i] =  true;
	}
	for(int i = totsize-s-ch  ; i < totsize-s + ch ; i++) {
		seq[i] =  true;
	}

	int sef = 0;
	for(int i = 0  ; i < totsize ; i++) {
		sef += seq[i];
	}
	cout << sef << endl;
	this->setseq(seq);


}

matrix<double> IsingPolymerPartial::density_fixedseq(matrix<double> &field, double &fe) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);
	vector1<double> field1(Nr);
	vector1<double> field2(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		field1[i]=field(i,0);
		field2[i]=field(i,1);
	}

	double fi=0.0;
	matrix<double> density = this->radialsolver_copolymer(field1,field2,fi);

	fe = fi;

	matrix<double> partialden(Nr,2);

	for(int i = 0 ; i < Nr ; i++) {
		partialden(i,0) = density(i,0);
		partialden(i,1) = density(i,1);		
	}


	return partialden;
}

double IsingPolymerPartial::freeenergy_fixedden(vector1<double> &density,vector1<double> &field,vector1<double> &prop,double &fe) {
	matrix<double> den(Nr,2);
	for(int i = 0 ; i < Nr ; i++) {
		den(i,0) = density[i]*prop[i];
		den(i,1) = density[i]*(1-prop[i]);
	}

	matrix<double> tfield(Nr,2);
	for(int i = 0 ; i < Nr ; i++) {
		tfield(i,0) =  -(1./beta)*log(exp(-beta*field[i])*prop[i]);
		tfield(i,1) =  -(1./beta)*log(exp(-beta*field[i])*(1-prop[i]));
	}

	return this->freeenergy_fixedseq(den, tfield,fe);
}

// double IsingPolymerPartial::freeenergy_fixedden(vector1<double> &field1, vector1<double> &fixed) {

// double fe=0.0;
// //cout << field << endl;
// matrix<double> den(Nr,2);
// den = density_fixedden(fixed,field1,fe);

// 	matrix<double> field(Nr,2);

// for(int i = 0 ; i < Nr ; i++) {
// 	field(i,0) = field1[i];
// 	field(i,1) = -(1./beta)*log(exp(-beta*fixed[i])-exp(-beta*field1[i]));		
// }


// double tot =0.0;

// for(int j = 0  ; j < 2 ; j++) {
// for(int i  = 0 ; i < Nr ; i++) {
// 	//cout << field(j,i) << endl;
// tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
// }
// }
	
// double tot2 = 0.0;

// for(int i  = 0 ; i < Nr ; i++) {
// 	double denatpoint = 0.0;
// 	for(int j = 0  ; j <2 ;j++)
// 		denatpoint+=den(i,j);

// 	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
// }

// double tot3 = 0.0;
// for(int k2 = 0 ; k2 < 2 ; k2++)
// for(int j = 0  ; j < 2 ;j++) {
// for(int i  = 0 ; i < Nr ; i++) {
// 	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
// }
// }

// double tot4 = 0.0;
// for(int j = 0  ; j <2 ; j++) {
// for(int i  = 0 ; i < Nr ; i++) {
// 	//cout << field(j,i) << endl;
// tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
// }
// }

// // cout << "free energy from partition function: " << fe << endl;
// // cout << "free energy due to w p term: " << -tot << endl;
// // cout << "hard sphere free energy: " << u0*tot2 << endl;
// // cout << "different density field interactions: " << v0*tot3 << endl;
// // cout << "fixed field free energy" << tot4 << endl;
// //cout << tot << endl;
// // cout << field << endl;

// // cout << tot2 << endl;
// // cout << tot3 << endl;

// return fe-tot+u0*tot2+v0*tot3+tot4;
// }





double IsingPolymerPartial::freeenergy_fixedseq(matrix<double> &field) {

double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,2);
den = density_fixedseq(field,fe);

double tot =0.0;

for(int j = 0  ; j < 2 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;

for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <2 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < 2 ; k2++)
for(int j = 0  ; j < 2 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

double tot4 = 0.0;
for(int j = 0  ; j <2 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
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

return fe-tot+u0*tot2+v0*tot3+tot4;
}

double IsingPolymerPartial::freeenergy_fixedseq(matrix<double> &den, matrix<double> &field, double &fe) {
double tot =0.0;

for(int j = 0  ; j < 2 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;

for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <2 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < 2 ; k2++)
for(int j = 0  ; j < 2 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

double tot4 = 0.0;
for(int j = 0  ; j <2 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
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

return fe-tot+u0*tot2+v0*tot3+tot4;
}



// vector1<double> IsingPolymerPartial::genconvolocal_fixedden(matrix<double> &density) {
// 	matrix<double> convo(Nr,2);

// 	for(int i = 0  ; i < Nr ; i++) {
// 		for(int j = 0  ; j < 2 ; j++) {
		
// 			double fac1 = 0.0;
// 			double fac2 = 0.0;

// 			for(int k1  = 0  ; k1 < 2 ; k1++)
// 				fac1+=density(i,k1); //excluded volume effects

// 			for(int k2 = 0 ; k2 < 2 ; k2++)
// 				fac2+=V(j,k2)*density(i,k2); //matrix of interactions between all the different fields

// 		convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
// 		}
// 	}

// 	vector1<double> field1(Nr);

// 	for(int i = 0  ; i <  Nr ; i++)
// 		field1[i] =  convo(i,0);


// 	return field1;
// }

matrix<double> IsingPolymerPartial::genconvolocal_fixedseq(matrix<double> &density) {
	matrix<double> convo(Nr,2);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < 2 ; j++) {
		
			double fac1 = 0.0;
			double fac2 = 0.0;

			for(int k1  = 0  ; k1 < 2 ; k1++)
				fac1+=density(i,k1); //excluded volume effects

			for(int k2 = 0 ; k2 < 2 ; k2++)
				fac2+=V(j,k2)*density(i,k2); //matrix of interactions between all the different fields

		convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
		}
	}


	return convo;
}

double IsingPolymerPartial::determine_stepsize2_fixedseq(matrix<double> &x, matrix<double> &a, double param, double fe1) {

	// x is our initial mean field;
	// a is our convoled meanfield
	//fe1 is our initial free energy
	//param is some parameter that determines step size
    double stepsize = 0.45*param;        
    
    matrix<double> temp2 = x+stepsize*(a-x);
    
    double fe2 = this->freeenergy_fixedseq(temp2);
    double param2;
    if(fe2-fe1<-1E-5) {
        double steptemp = stepsize*2;

        matrix<double> gemp2 = x+steptemp*(a-x);                
        double ftemp = this->freeenergy_fixedseq(gemp2);
        
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
            fe2 = this->freeenergy_fixedseq(gemp2);
            iter++;
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

// double IsingPolymerPartial::determine_stepsize2_fixedden(vector1<double> &x, vector1<double> &a, vector1<double> &fixed, double param, double fe1) {

// 	// x is our initial mean field;
// 	// a is our convoled meanfield
// 	//fe1 is our initial free energy
// 	//param is some parameter that determines step size


//     double stepsize = 0.45*param;        
    
//     vector1<double> temp2 = x+stepsize*(a-x);
    
//     double fe2 = this->freeenergy_fixedden(temp2,fixed);
//     double param2;
//     if(fe2-fe1<-1E-5) {
//         double steptemp = stepsize*2;

//         vector1<double> gemp2 = x+steptemp*(a-x);                
//         double ftemp = this->freeenergy_fixedden(gemp2,fixed);

//         // cout << "branch 1 free energy: " << ftemp << endl;
        
//         if(ftemp < fe2) param2 = steptemp;
//         else {
//             double a = -((-(fe1*stepsize) + ftemp*stepsize + fe1*steptemp - fe2*steptemp)/
//  (stepsize*(stepsize - steptemp)*steptemp));
//             double b = -((fe1*SQR(stepsize) - ftemp*SQR(stepsize) - 
//    fe1*SQR(steptemp) + fe2*SQR(steptemp))/
//  (stepsize*(stepsize - steptemp)*steptemp));
//             param2 = -b/(2.*a);                
//         }
// //            cout << fe1 << "," << fe2 << "," << ftemp << endl;
// //            cout << 0 << "," << stepsize << "," << steptemp << endl;
        
//     }
//     else if (fe1-fe2 < -1E-5) {
//     	//cout << "branch 2" << endl;
//         int iter = 0;
//         bool lesser = true;
//         double steptemp = stepsize;
//         double ftemp;
//         while(fe1 < fe2) {
//             if(iter > 6) {lesser = false; break; }
//             ftemp = fe2;
//             steptemp = stepsize;
//             stepsize/=2;
//             vector1<double> gemp2 = x+stepsize*(a-x);
//             fe2 = this->freeenergy_fixedden(gemp2,fixed);
//             //cout << fe2 << ", ";
//             iter++;
//         }
        
//         if(lesser) {
//             double a = -((-(fe1*stepsize) + ftemp*stepsize + fe1*steptemp - fe2*steptemp)/
//  (stepsize*(stepsize - steptemp)*steptemp));
//             double b = -((fe1*SQR(stepsize) - ftemp*SQR(stepsize) - 
//    fe1*SQR(steptemp) + fe2*SQR(steptemp))/
//  (stepsize*(stepsize - steptemp)*steptemp));
//             param2 = -b/(2.*a);
//         }
//         else {
//             param2 = stepsize/2.;
//         }
// //            cout << fe1 << "," << fe2 << "," << ftemp << endl;
// //            cout << 0 << "," << stepsize << "," << steptemp << endl;
//         // have fe1, fe2, ftemp
//         // with 0  , stepsize, steptemp
        
//     }
//     else {
//         cout << "fe1=fe2" << endl;
//         param2 = stepsize;
        
//     }

//     return param2;        
// }

vector1<double> IsingPolymerPartial::genconvolocal_fixedden(vector1<double> &density,vector1<double> &field,vector1<double> &prop) {
	//matrix<double> convo(Nr,2);
	vector1<double> cv1(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		cv1[i]=-(1./beta)*density[i]*(1+log(exp(-beta*field[i])*prop[i]))+(1./beta)*density[i]*(1+log(exp(-beta*field[i])*(1.-prop[i])));
	}

	vector1<double> cv2(Nr);
	for(int i = 0 ; i < Nr ; i++) {
		cv2[i] = -(V(0,0)*SQR(density[i])*prop[i]+V(0,1)*SQR(density[i])*(1-2*prop[i])-V(1,1)*SQR(density[i])*(1-prop[i]));
	}

	vector1<double> cv3(Nr);
	for(int i = 0  ; i < Nr ; i++) {
		cv3[i] = -density[i]*fixed_field(i,0)+density[i]*fixed_field(i,1);
	}

	// for(int i = 0  ; i < Nr ; i++) { 
	// 	for(int j = 0  ; j < 2 ; j++) {
		
	// 		double fac1 = 0.0;
	// 		double fac2 = 0.0;

	// 		for(int k1  = 0  ; k1 < 2 ; k1++)
	// 			fac1+=density(i,k1); //excluded volume effects

	// 		for(int k2 = 0 ; k2 < 2 ; k2++)
	// 			fac2+=V(j,k2)*density(i,k2); //matrix of interactions between all the different fields

	// 	convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
	// 	}
	// }

	return cv1+v0*cv2+cv3;
}

void IsingPolymerPartial::corrfunc(vector1<double> &f) {
	for(int i = 0 ; i < f.getsize() ; i++ ) {
		if(f[i]<0) f[i]=1E-10;
		else if(f[i]>1) f[i]=1.-1E-10;
		else {}
	}
}

void IsingPolymerPartial::iteratefield_fixedden(vector1<double> &density, vector1<double> &field,vector1<double> &prop, double dt, double &fen, double fe, double &stepsize2) {
	//density is the total density rA+rB

	double fen1  = this->freeenergy_fixedden(density,field,prop,fe);
	//f2 = rA/rA+rB
	cout << "free energy: " << fen1 << endl;

	fen = fen1;
	vector1<double> convolocal = this->genconvolocal_fixedden(density,field,prop);
	// cout << prop << endl;
	// pausel();

	double stepsize = dt;
	stepsize2 = stepsize;
	// matrix<double> den(Nr,2);
	prop = 	prop+stepsize*(convolocal);

	this->corrfunc(prop);

	// cout << prop << endl;
	// pausel();
	// for(int i = 0 ; i < Nr ; i++) {
	// 	den(i,0) = density[i]*f2[i];
	// 	den(i,1) = density[i]*(1-f2[i]);
	// }

	// matrix<double> tfield(Nr,2);
	// for(int i = 0 ; i < Nr ; i++) {
	// 	field(i,0) =  -(1./beta)*log(exp(-beta*field[i])*prop[i]);
	// 	field(i,1) =  -(1./beta)*log(exp(-beta*field[i])*(1-prop[i]));
	// }
	//f2 is the proprtion field

}

// void IsingPolymerPartial::iteratefield_fixedden(vector1<double> &field, vector1<double> &field1,  double dt, double &fen, double &stepsize2) { //field is the fixed field
// 	double fe=0.0;
// 	//cout << field << endl;


// 	matrix<double> den(Nr,2);
// 	den = density_fixedden(field,field1,fe);

// 	//cout << "output fe: " << fe << endl;

// //	outfunc(den,"den");
// //	outfunc(field,"field");
// 	matrix<double> fieldm(Nr,2);
// 	vector1<double> field2(Nr);
// 	for(int i = 0  ; i < Nr ; i++) {
// 		field2[i] = -(1./beta)*log(exp(-beta*field[i])-exp(-beta*field1[i]));
// 	}
// 	for(int i = 0 ; i < Nr ; i++) {
// 		fieldm(i,0)=field1[i];
// 		fieldm(i,1)=field2[i];
// 	}

// 	double fen1  = this->freeenergy_fixedseq(den,fieldm,fe);
// 	cout << "free energy: " << fen1 << endl;
// 	//pausel();

// 	fen = fen1;
// 	vector1<double> convolocal = this->genconvolocal_fixedden(den);
// 	// outfunc(fixed_field,"ff");
// 	// outfunc(den,"den");
// 	// outfunc(field1,"field1");
// 	// outfunc(get_other_field(field1,field),"field2");
// 	// outfunc(convolocal,"convolocal");
// 	// outfunc(get_other_field(convolocal,field),"convolocalother");
// 	// outfunc(this->genconvolocal_fixedseq(den),"matconvolocal");

// 	// vector1<double> v2 = field1+0.0001*(convolocal-field1);
// 	// cout << this->freeenergy_fixedden(field1,field) << endl;
// 	// cout << this->freeenergy_fixedden(v2,field) << endl;

// //	matrix<double> den2 =  density_fixedden(field,v2,fe);
// 	// outfunc(den2,"den2");

// 	// outfunc(v2,"field1n");
// 	// outfunc(get_other_field(v2,field),"field2n");
// 	// pausel();
// //	outfunc(convolocal,"convolocal");
	
// 	double stepsize = this->determine_stepsize2_fixedden(field1,convolocal,field,dt,fen1);
// 	cout << "steps: " << dt << " " << stepsize << endl;
// 	stepsize2 = stepsize;
	
// 	// outfunc(field,"f1");
// 	// outfunc((-stepsize)*field,"t1");

// 	// outfunc(stepsize*convolocal,"t2");
// 	// pausel();
// 	field1 = field1+stepsize*(-field1+convolocal);
// }

void IsingPolymerPartial::iteratefield_fixedseq(matrix<double> &field, double dt, double &fen, double &stepsize2) {
	double fe=0.0;
	//cout << field << endl;


	matrix<double> den(Nr,2);
	den = density_fixedseq(field,fe);

	//cout << "output fe: " << fe << endl;

//	outfunc(den,"den");
//	outfunc(field,"field");
	double fen1  = this->freeenergy_fixedseq(den,field,fe);
	

	// cout <<  "fe1: " << this->freeenergy_copolymer(field) << endl;
	cout << "free energy: " << fen1 << endl;
	//pausel();

	fen = fen1;
	matrix<double> convolocal = this->genconvolocal_fixedseq(den);
//	outfunc(convolocal,"convolocal");
	
	double stepsize = this->determine_stepsize2_fixedseq(field,convolocal,dt,fen1);
	cout << "steps: " << dt << " " << stepsize << endl;
	stepsize2 = stepsize;
	
	// outfunc(field,"f1");
	// outfunc((-stepsize)*field,"t1");

	// outfunc(stepsize*convolocal,"t2");
	// pausel();
	field = field+stepsize*(-field+convolocal);
}


void IsingPolymerPartial::run_fixedseq(matrix<double> &field, int tot, double dt) {
	double fen = 0.0;
	for(int i = 0 ; i < tot ; i++)
	{
		double fi = fen;
		cout << i << endl;
		double stepsize2 = 0.0;
		this->iteratefield_fixedseq(field,dt,fen,stepsize2);

		double ff = fen;
		//cout << abs(fi-ff) << endl;
		if(abs(fi-ff)<stepsize2) { cout << "convergence test made" << endl; break; }
		else if(stepsize2<1E-7) {cout << "vanishing step" << endl; break; }
		else{ }
	}
}


void IsingPolymerPartial::run_fixedden(vector1<double> &f,vector1<double> &prop, int tot, double dt) {
	double fen = 0.0;
	double fe = 0.0;
	vector1<double> density = this->radialsolver(f,fe);


	//fe = fi;

	for(int i = 0 ; i < tot ; i++)	{
		//f is the original field;

		double fi = fen;
		cout << i << endl;
		double stepsize2 = 0.0;
		//f2 is THE PROPORTION FIELD
		this->iteratefield_fixedden(density,f,prop,dt,fen,fe,stepsize2);
		

		double ff = fen;
		//cout << abs(fi-ff) << endl;
		if(abs(fi-ff)<1E-6) { cout << "convergence test made" << endl; break; }
		else if(stepsize2<1E-6) {cout << "vanishing step" << endl; break; }
		else{ }
	}
//cout << density << endl;
}

void IsingPolymerPartial::run_mixed() {
	matrix<double> fixfie = fixed_field;
	matrix<double> start(Nr,2);
	matrix<double> den(Nr,2);
	double stepsize = 0.01;
	int runtime = 10000;
	
//	matrix<double> field_fixseq(Nr,2);
	int iter = 0;

	//this->setfixed_field(field_fixseq);
	for(;;){
	this->run_fixedseq(start,runtime,stepsize);

	vector1<double> field(Nr);
	for(int i = 0 ; i < Nr ; i++ ) {
		field[i]=-(1./beta)*log(exp(-beta*start(i,0))+exp(-beta*start(i,1)));
	}
	double fe=0;

	den = this->density_fixedseq(start,fe);

	string s1 = "denfixseq";
	stringstream ss1;
	ss1 << iter;
	string tf=s1+ss1.str();
	outfunc(den,tf);

	vector1<double> prop(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		prop[i] = den(i,0)/(den(i,0)+den(i,1));
	}
	//this->setfixed_field(fixfie);

	this->run_fixedden(field,prop,runtime,stepsize);
	
	for(int i = 0 ; i < Nr ; i++) {
		double fac  = den(i,0)+den(i,1);
		den(i,0) = fac*prop[i];
		den(i,1) = fac*(1-prop[i]);
	}

	string s2 = "denfixden";
	string tf2=s2+ss1.str();
	outfunc(den,tf2);

	for(int i = 0 ; i < Nr ; i++) {
		start(i,0) =  -(1./beta)*log(exp(-beta*field[i])*prop[i]);
		start(i,1) =  -(1./beta)*log(exp(-beta*field[i])*(1-prop[i]));
	}
	this->determinesequence(den,250);
	
	pausel();

	iter++;
	}




}