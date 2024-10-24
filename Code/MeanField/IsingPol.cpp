double sphericalbessel0(double r) {
	if(r<1E-10) return 1.0;
	else return sin(r)/r;
}


//if cop = 1, copolymer , if = 0 normal system
IsingPolymer::IsingPolymer(int Nrr, int no_typess, int copp = 0) : radialpoints(vector1<double>(Nrr)), propd(vector1<double>(Nrr)), amat(matrix<double>(Nrr, Nrr)), bmat(matrix<double>(Nrr, Nrr)), bigD(matrix<double>(Nrr, Nrr)), bigDL(matrix<LogSpace>(Nrr, Nrr)), V(matrix<double>(no_typess + copp, no_typess + copp)), fixed_field(matrix<double>(Nrr, no_typess + copp))
{
	cop = copp;
	ds = 1.0;
	stot =  1000.0;
	totpolymers = 1.0;
	b = 1.0;
	R = 10.0;
	Nr = Nrr;
	no_types = no_typess;
	polymer_or_fluid = true; //false is a fluid, true is a polymer
	int totsize = 1+(int)(stot/ds);
	seq = new vector1<bool>(totsize);

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


	for(int i = 0 ; i < Nr ; i++) {
		for(int j = 0 ; j < Nr ; j++) {
			bigDL(i,j)=LogSpace(bigD(i,j));
		}
	}
	


	//bigD = bmat*(tempd*amat);

	beta = 1.0;
	u0 = 1.0;
	v0 = 1.0;
}

IsingPolymer::~IsingPolymer() {
	delete seq;
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
	for (int i = 0; i < Nr; i++)
	{
		for (int j = 0; j < Nr; j++)
		{
			bigDL(i, j) = LogSpace(bigD(i, j));
		}
	}

	int totsize = 1+(int)(stot/ds);
	seq = new vector1<bool>(totsize);



};

void IsingPolymer::setstot(double stott) {
	stot = stott;
	int totsize = 1+(int)(stot/ds);
	seq = new vector1<bool>(totsize);
}

void IsingPolymer::setseq(vector1<bool> &seqq) {
	//stot = stott;
	int totsize = 1+(int)(stot/ds);
	if(totsize!=seqq.getsize() ) error("error in setting sequence");
	*seq = seqq;

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
	for (int i = 0; i < Nr; i++)
	{
		for (int j = 0; j < Nr; j++)
		{
			bigDL(i, j) = LogSpace(bigD(i, j));
		}
	}
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
	for (int i = 0; i < Nr; i++)
	{
		for (int j = 0; j < Nr; j++)
		{
			bigDL(i, j) = LogSpace(bigD(i, j));
		}
	}
}
void IsingPolymer::setV(matrix<double> &V2) {
	if(V2.getnrows() != no_types+cop || V2.getncols() != no_types+cop) error("incorrect size for matrix of interactions");
	for(int i = 0  ; i < no_types+cop ; i++) {
		for(int j = 0  ; j < no_types+cop ; j++) {
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

void IsingPolymer::setfixed_field(matrix<double> &fixed_fieldd) {
	if(fixed_fieldd.getnrows() !=  Nr) error("wrong dimension fixed field setting");

	if(fixed_fieldd.getncols() !=  no_types+cop) error("wrong dimension fixed field setting");

	fixed_field = fixed_fieldd;
}

void IsingPolymer::settotpolymers(double totpolymerss) {
	totpolymers = totpolymerss;
}

template <class T>
T IsingPolymer::integrate(vector1<T> &m) {
	T tot = 0.0;
	for(int i  = 0 ; i < Nr ; i++)
	tot+=m[i]*T(4*pi*SQR(radialpoints[i])*dr);

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

	
	norm=-totpolymers*log((tot/(long double)stot)/(long double)V);

	//cout << norm << endl;

	den2*=((long double)(totpolymers*stot)/tot);
	for(int i = 0  ; i < Nr ; i++) {
		den[i]=den2[i];
	}

	//return stot*(den/tot);
	return den;

}

matrix<double> IsingPolymer::radialsolver_copolymer(vector1<double> &field1, vector1<double> &field2, double &norm) { //MAKE SURE SEQUENCE IS SYMMETRIC
	//vector1<double> den(Nr);
	// vector1<double> propw(Nr);
	int totsize = 1+(int)(stot/ds);
	// for(int i  = 0 ; i < Nr ; i++) {
	// 	propw[i]=exp(-field[i]*ds/2.);
	// }
	if ((*seq).getsize()!=totsize) error("incorrect sequence definition");

	LogSpace aprop = 0.0;

	for(int i = 0  ; i < totsize ; i++) {
		aprop += (double)(*seq)[i];
	}

	LogSpace bprop = totsize - aprop;



	matrix<double> propw1(Nr,Nr);

	matrix<double> propw2(Nr,Nr);
	for(int i  = 0 ; i < Nr ; i++) {
	
		propw1(i,i)=exp(-field1[i]*ds/2.);
		propw2(i,i)=exp(-field2[i]*ds/2.);
	}
	matrix<long double> itmat1 = (propw1)*(bigD*propw1);

	matrix<long double> itmat2 = (propw2)*(bigD*propw2);

	matrix<LogSpace> itmat1LS(Nr,Nr);
	matrix<LogSpace> itmat2LS(Nr,Nr);

	for(int i =  0 ; i < Nr ; i++)
		for(int j = 0  ; j < Nr ; j++) {
			itmat1LS(i, j) << itmat1(i, j);
			itmat2LS(i, j) << itmat2(i, j);
		}


	//int totsize = 1+(int)(stot/ds);

	matrix<LogSpace> store(totsize,Nr);
	for(int j  = 0 ; j < Nr-1 ; j++) {
		store(0,j)=LogSpace(1.0);
	}
	store(0,Nr-1)=LogSpace(0.0);

	vector1<LogSpace> q13(Nr);
	for(int j = 0  ; j < Nr ; j++) {
		q13[j] = store(0,j);
	}

		// cout << itmat2LS << endl;

		// cout << itmat1LS << endl;



	for(int i = 0  ; i < totsize-1 ; i++) {

		// q13=propw&q13;
		// q13=amat*q13;
		// q13=propd&q13;
		// q13=bmat*q13;
		// q13=propw&q13;


		if((*seq)[i]) {
		q13=(itmat1LS*q13);
		}
		else {
		q13=(itmat2LS*q13);	
		}


		

		for(int j = 0  ; j < Nr ; j++) {
			store(i+1,j)=q13[j];
		}
		
	}






	//vector1<double> den(Nr);
	vector1<LogSpace> den2(Nr);
	for(int i = 0 ; i < totsize ; i++) {
		for(int j = 0 ; j < Nr ; j++) {		
		den2[j]+=(store(i,j)*store(totsize-1-i,j));
		}
	}
	LogSpace tot(0.0);
	for(int i  = 0 ; i < Nr ; i++)
	tot+=den2[i]*LogSpace(4.*pi*SQR(radialpoints[i])*dr);

	//cout << "tot: " << tot << endl;	

	double V=(4./3.)*pi*R*R*R;
	// LogSpace qint(0.0);
	// for(int i  = 0 ; i < Nr ; i++)
	// qint+=store(totsize-1,i)*LogSpace(4*pi*SQR(radialpoints[i])*dr);
	
	norm = -totpolymers *log((tot / LogSpace(stot)) / LogSpace(V));
	// cout << qint << endl;
	// norm=(-totpolymers/(long double)V)*(log(tot)-log((long double)stot));


	
	//cout << -totpolymers * log((tot / LogSpace(stot)) / (long double)V) << endl;
	// cout << totpolymers*log((tot)/LogSpace(stot)) / V << endl;
	// cout << log(tot) - log(stot) << endl;
	// cout << norm << endl;
	// pausel();
	// //cout << norm << endl;

	vector1<LogSpace> denA(Nr);
	vector1<LogSpace> denB(Nr);
	for(int i = 0 ; i < totsize ; i++) {
		if((*seq)[i]) {
		for(int j = 0 ; j < Nr ; j++) {		
		denA[j]+=(store(i,j)*store(totsize-1-i,j));
		}
		}
		else{
		for(int j = 0 ; j < Nr ; j++) {
		denB[j]+=(store(i,j)*store(totsize-1-i,j));
		}
		}
	}

	LogSpace na = this->integrate(denA);

	LogSpace nb = this->integrate(denB);


	// cout << na << endl;

	// cout << nb << endl;

	// cout << na + nb << endl;

	// cout << tot << endl;
	// cout << aprop << endl;
	// cout << bprop << endl;


	denA*=((totpolymers*aprop)/na);
	denB*=((totpolymers*bprop)/nb);
	
	matrix<double> den(Nr,2);
	for(int i = 0  ; i < Nr ; i++) {
		den(i,0)=denA[i].getdouble();
		den(i,1)=denB[i].getdouble();
	}
	// cout << aprop << endl;
	// cout << bprop << endl;
	// cout << this->integrate(denA) << endl;
	// cout << this->integrate(denB) << endl;
	// pausel();
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

matrix<double> IsingPolymer::densityfromfield(matrix<double> &field) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field(i,j));
		}
		totalfield[i]=-(1./beta)*log(norm[i]);
	}

	double fi=0.0;
	vector1<double> density = this->radialsolver(totalfield,fi);


	matrix<double> partialden(Nr,no_types);
	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density[i]*exp(-beta*field(i,j))/norm[i];		
		}
	}
	return partialden;
}

matrix<double> IsingPolymer::densityfromfield_copolymer(matrix<double> &field) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);
	vector1<double> field1(Nr);
	matrix<double> field2(Nr,no_types);

	for(int i = 0  ; i < Nr ; i++) {
		field1[i]=field(i,0);
		for(int j = 0  ; j < no_types ; j++) {
		field2(i,j)=field(i,j+1);
		}
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field2(i,j));
		}
		totalfield[i]=-(1./beta)*log(norm[i]);
	}


	double fi=0.0;
	matrix<double> density = this->radialsolver_copolymer(field1,totalfield,fi);


	matrix<double> partialden(Nr,no_types+1);

	for(int i = 0 ; i < Nr ; i++) {
		partialden(i,0) = density(i,0);
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j+1)=density(i,1)*exp(-beta*field2(i,j))/norm[i];		
		}
	}

	return partialden;
}

matrix<double> IsingPolymer::densityfromfield_dcopolymer(matrix<double> &field) {
	vector1<double> totalfield1(Nr),totalfield2(Nr);
	vector1<double> norm1(Nr),norm2(Nr);
	matrix<double> field1(Nr,no_types);
	matrix<double> field2(Nr,no_types);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < no_types ; j++) {
		field1(i,j)=field(i,j);
		}
		for(int j = 0  ; j < no_types ; j++) {
		field2(i,j)=field(i,j+no_types);
		}
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm1[i]+=exp(-beta*field1(i,j));
		}
		totalfield1[i]=-(1./beta)*log(norm1[i]);
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm2[i]+=exp(-beta*field2(i,j));
		}
		totalfield2[i]=-(1./beta)*log(norm2[i]);
	}	


	double fi=0.0;
	matrix<double> density = this->radialsolver_copolymer(totalfield1,totalfield2,fi);


	matrix<double> partialden(Nr,no_types+no_types);

	for(int i = 0 ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density(i,0)*exp(-beta*field1(i,j))/norm1[i];		
		}
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j+no_types)=density(i,1)*exp(-beta*field2(i,j))/norm2[i];		
		}
	}

	return partialden;
}




matrix<double> IsingPolymer::densityfromfield_fluid(matrix<double> &field) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field(i,j));
		}
		totalfield[i]=-(1./beta)*log(norm[i]);
	}

	double fi=0.0;
	vector1<double> density(Nr);
	for(int i = 0  ; i < Nr ; i++) {
		density[i]=exp(-beta*totalfield[i]);
	}
	double normal=this->integrate(density);

	density *= (stot*totpolymers)/normal;

	matrix<double> partialden(Nr,no_types);
	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density[i]*exp(-beta*field(i,j))/norm[i];		
		}
	}

	return partialden;
}

matrix<double> IsingPolymer::densityfromfield_fluid(matrix<double> &field, double &fe) {
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field(i,j));
		}
		totalfield[i]=-(1./beta)*log(norm[i]);
	}

	double fi=0.0;
	vector1<double> density(Nr);
	for(int i = 0  ; i < Nr ; i++) {
		density[i]=exp(-beta*totalfield[i]);
	}
	double normal=this->integrate(density);

	density *= (stot*totpolymers)/normal;

	fe = -stot*totpolymers*log(normal);

	matrix<double> partialden(Nr,no_types);
	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density[i]*exp(-beta*field(i,j))/norm[i];		
		}
	}

	return partialden;
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

	vector1<double> density = this->radialsolver(totalfield,fi);

	matrix<double> partialden(Nr,no_types);
	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density[i]*exp(-beta*field(i,j))/norm[i];		
		}
	}

	return partialden;
}

matrix<double> IsingPolymer::densityfromfield_copolymer(vector1<double> &field1,matrix<double> &field2, double &fe) { //copolymer field theory
	vector1<double> totalfield(Nr);
	vector1<double> norm(Nr);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm[i]+=exp(-beta*field2(i,j));
		}
		totalfield[i]=-(1./beta)*log(norm[i]);
	}



	double fi=0.0;
	matrix<double> density = this->radialsolver_copolymer(field1,totalfield,fi);

	fe=fi;

	matrix<double> partialden(Nr,no_types+1);

	for(int i = 0 ; i < Nr ; i++) {
		partialden(i,0) = density(i,0);
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j+1)=density(i,1)*exp(-beta*field2(i,j))/norm[i];		
		}
	}

	return partialden;
}

matrix<double> IsingPolymer::densityfromfield_dcopolymer(matrix<double> &field1,matrix<double> &field2, double &fe) { //copolymer field theory
	vector1<double> totalfield1(Nr),totalfield2(Nr);
	vector1<double> norm1(Nr),norm2(Nr);
	// matrix<double> field1(Nr,no_types);
	// matrix<double> field2(Nr,no_types);

	// for(int i = 0  ; i < Nr ; i++) {
	// 	for(int j = 0  ; j < no_types ; j++) {
	// 	field1(i,j)=field(i,j);
	// 	}
	// 	for(int j = 0  ; j < no_types ; j++) {
	// 	field2(i,j)=field(i,j+no_types);
	// 	}
	// }

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm1[i]+=exp(-beta*field1(i,j));
		}
		totalfield1[i]=-(1./beta)*log(norm1[i]);
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			norm2[i]+=exp(-beta*field2(i,j));
		}
		totalfield2[i]=-(1./beta)*log(norm2[i]);
	}	


	double fi=0.0;
	matrix<double> density = this->radialsolver_copolymer(totalfield1,totalfield2,fi);
	fe = fi;


	matrix<double> partialden(Nr,no_types+no_types);

	for(int i = 0 ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j)=density(i,0)*exp(-beta*field1(i,j))/norm1[i];		
		}
	}

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0 ; j < no_types ; j++) {
			partialden(i,j+no_types)=density(i,1)*exp(-beta*field2(i,j))/norm2[i];		
		}
	}

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


		convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
		}
	}
	return convo;
}

matrix<double> IsingPolymer::genconvolocal_copolymer(matrix<double> &density) {
	matrix<double> convo(Nr,no_types+1);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < no_types+1 ; j++) {
		
			double fac1 = 0.0;
			double fac2 = 0.0;

			for(int k1  = 0  ; k1 < no_types+1 ; k1++)
				fac1+=density(i,k1); //excluded volume effects

			for(int k2 = 0 ; k2 < no_types+1 ; k2++)
				fac2+=V(j,k2)*density(i,k2); //matrix of interactions between all the different fields

		convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
		}
	}


	return convo;
}

matrix<double> IsingPolymer::genconvolocal_dcopolymer(matrix<double> &density) {
	matrix<double> convo(Nr,no_types+no_types);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < no_types+no_types ; j++) {
		
			double fac1 = 0.0;
			double fac2 = 0.0;

			for(int k1  = 0  ; k1 < no_types+no_types ; k1++)
				fac1+=density(i,k1); //excluded volume effects

			for(int k2 = 0 ; k2 < no_types+no_types ; k2++)
				fac2+=V(j,k2)*density(i,k2); //matrix of interactions between all the different fields

		convo(i,j)=u0*fac1+v0*fac2+fixed_field(i,j);
		}
	}


	return convo;
}


double IsingPolymer::freeenergy(matrix<double> &field) {
double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,no_types);
if(polymer_or_fluid) {
den = densityfromfield(field,fe);
}
else{
den = densityfromfield_fluid(field,fe);	
}
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


double tot4 = 0.0;
for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3+tot4;
}

double IsingPolymer::freeenergy_copolymer(vector1<double> &field1,matrix<double> &field2) {
double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,no_types+1);

den = densityfromfield_copolymer(field1,field2,fe);

//cout << den << endl;

double tot =0.0;


for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,0)*field1[i]*4.*pi*SQR(radialpoints[i])*dr;
}


for(int j = 1  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field2(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+1 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+1 ; k2++)
for(int j = 0  ; j <no_types+1 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3+tot4;
}

double IsingPolymer::freeenergy_dcopolymer(matrix<double> &field1,matrix<double> &field2) {
double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,no_types+no_types);

den = densityfromfield_dcopolymer(field1,field2,fe);

//cout << den << endl;

double tot =0.0;


for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field1(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

for(int j = no_types  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field2(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+no_types ; k2++)
for(int j = 0  ; j <no_types+no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3+tot4;
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

double tot4 = 0.0;
for(int j = 0  ; j <no_types ; j++) {
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

double IsingPolymer::freeenergy_copolymer(matrix<double> &den, vector1<double> &field1,matrix<double> &field2, double &fe) {

double tot =0.0;


for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,0)*field1[i]*4.*pi*SQR(radialpoints[i])*dr;
}


for(int j = 1  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field2(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+1 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+1 ; k2++)
for(int j = 0  ; j <no_types+1 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3+tot4;

}

double IsingPolymer::freeenergy_dcopolymer(matrix<double> &den, matrix<double> &field1,matrix<double> &field2, double &fe) {

double tot =0.0;


for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field1(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

for(int j = no_types  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field2(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+no_types ; k2++)
for(int j = 0  ; j <no_types+no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;

return fe-tot+u0*tot2+v0*tot3+tot4;

}



double IsingPolymer::freeenergy_copolymer(matrix<double> &field, bool show = false) {

vector1<double> field1(Nr);
matrix<double> field2(Nr,no_types);

for(int i = 0  ; i < Nr ; i++) {
	field1[i]=field(i,0);
	for(int j = 0  ; j < no_types ; j++) {
	field2(i,j)=field(i,j+1);
	}
}

double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,no_types+1);

den = densityfromfield_copolymer(field1,field2,fe);

//cout << den << endl;

double tot =0.0;


for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,0)*field1[i]*4.*pi*SQR(radialpoints[i])*dr;
}


for(int j = 1  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field2(i,j-1)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+1 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+1 ; k2++)
for(int j = 0  ; j <no_types+1 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

if(show) {
cout << "free energy from partition function: " << fe << endl;
cout << "free energy due to w p term: " << -tot << endl;
cout << "hard sphere free energy: " << u0*tot2 << endl;
cout << "different density field interactions: " << v0*tot3 << endl;
cout << "fixed field contrib: " << tot4 << endl;
}

return fe-tot+u0*tot2+v0*tot3+tot4;
}


double IsingPolymer::freeenergy_dcopolymer(matrix<double> &field, bool show = false) {

matrix<double> field1(Nr,no_types);
matrix<double> field2(Nr,no_types);

for(int i = 0  ; i < Nr ; i++) {
	for(int j = 0  ; j < no_types ; j++) {
	field1(i,j)=field(i,j);
	}
	for(int j = 0  ; j < no_types ; j++) {
	field2(i,j)=field(i,j+no_types);
	}
}

double fe=0.0;
//cout << field << endl;
matrix<double> den(Nr,no_types+no_types);

den = densityfromfield_dcopolymer(field1,field2,fe);

//cout << den << endl;
double tot =0.0;


for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field1(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

for(int j = 0  ; j <no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j+no_types)*field2(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+no_types ; k2++)
for(int j = 0  ; j <no_types+no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

if(show) {
cout << "free energy from partition function: " << fe << endl;
cout << "free energy due to w p term: " << -tot << endl;
cout << "hard sphere free energy: " << u0*tot2 << endl;
cout << "different density field interactions: " << v0*tot3 << endl;
cout << "fixed field contrib: " << tot4 << endl;
}

return fe-tot+u0*tot2+v0*tot3+tot4;
}



double IsingPolymer::freeenergy_copolymer(matrix<double> &den, matrix<double> &field, double &fe) {


double tot =0.0;


for(int j = 0  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+1 ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+1 ; k2++)
for(int j = 0  ; j <no_types+1 ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+1 ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// cout << "free energy from partition function: " << fe << endl;
// cout << "free energy due to w p term: " << -tot << endl;
// cout << "hard sphere free energy: " << u0*tot2 << endl;
// cout << "different density field interactions: " << v0*tot3 << endl;
// cout << "fixed field contrib: " << tot4 << endl;


return fe-tot+u0*tot2+v0*tot3+tot4;

}


double IsingPolymer::freeenergy_dcopolymer(matrix<double> &den, matrix<double> &field, double &fe) {
// double fe=0.0;
// //cout << field << endl;
// matrix<double> den(Nr,no_types+1);

// den = densityfromfield_copolymer(field1,field2,fe);

//cout << den << endl;

//double tot1 =0.0;
double tot =0.0;


for(int j = 0  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}

// for(int j = no_types  ; j <no_types+no_types ; j++) {
// for(int i  = 0 ; i < Nr ; i++) {
// 	//cout << field(j,i) << endl;
// tot+=den(i,j)*field(i,j)*4*pi*SQR(radialpoints[i])*dr;
// }
// }
	
double tot2 = 0.0;
for(int i  = 0 ; i < Nr ; i++) {
	double denatpoint = 0.0;
	for(int j = 0  ; j <no_types+no_types ;j++)
		denatpoint+=den(i,j);

	tot2+=SQR(denatpoint)*4*pi*SQR(radialpoints[i])*dr;
	
}

double tot3 = 0.0;
for(int k2 = 0 ; k2 < no_types+no_types ; k2++)
for(int j = 0  ; j <no_types+no_types ;j++) {
for(int i  = 0 ; i < Nr ; i++) {
	tot3+=V(j,k2)*den(i,k2)*den(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


double tot4 = 0.0;
for(int j = 0  ; j <no_types+no_types ; j++) {
for(int i  = 0 ; i < Nr ; i++) {
	//cout << field(j,i) << endl;
tot4+=den(i,j)*fixed_field(i,j)*4*pi*SQR(radialpoints[i])*dr;
}
}


return fe-tot+u0*tot2+v0*tot3+tot4;


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

double IsingPolymer::determine_stepsize2_copolymer(matrix<double> &x, matrix<double> &a, double param, double fe1) {

	// x is our initial mean field;
	// a is our convoled meanfield
	//fe1 is our initial free energy
	//param is some parameter that determines step size
    double stepsize = 0.45*param;        
    
    matrix<double> temp2 = x+stepsize*(a-x);
    
    double fe2 = this->freeenergy_copolymer(temp2);

   
  //  pausel();

    double param2;
    if(fe2-fe1<-1E-5) {
    	//cout << "branch 1" << endl;
        double steptemp = stepsize*2;

        matrix<double> gemp2 = x+steptemp*(a-x); 
        //cout << gemp2 << endl;

        double ftemp = this->freeenergy_copolymer(gemp2);
        //cout << "branch 1 free energy: " << ftemp << endl;

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
    	//cout << "branch 2" << endl;
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
            fe2 = this->freeenergy_copolymer(gemp2);
            //cout << fe2 << ", ";
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
    	//cout << "branch 3" << endl;
        cout << "fe1=fe2" << endl;
        param2 = stepsize;
        
    }

    return param2;        
} 


double IsingPolymer::determine_stepsize2_dcopolymer(matrix<double> &x, matrix<double> &a, double param, double fe1) {

	// x is our initial mean field;
	// a is our convoled meanfield
	//fe1 is our initial free energy
	//param is some parameter that determines step size
    double stepsize = 0.45*param;        
    
    matrix<double> temp2 = x+stepsize*(a-x);
    
    double fe2 = this->freeenergy_dcopolymer(temp2);

   
  //  pausel();

    double param2;
    if(fe2-fe1<-1E-5) {
    	//cout << "branch 1" << endl;
        double steptemp = stepsize*2;

        matrix<double> gemp2 = x+steptemp*(a-x); 
        //cout << gemp2 << endl;

        double ftemp = this->freeenergy_dcopolymer(gemp2);
        cout << "branch 1 free energy: " << ftemp << endl;

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
    	cout << "branch 2" << endl;
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
            fe2 = this->freeenergy_dcopolymer(gemp2);
            cout << fe2 << ", ";
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
    	cout << "branch 3" << endl;
        cout << "fe1=fe2" << endl;
        param2 = stepsize;
        
    }

    return param2;        
} 


void IsingPolymer::iteratefield(matrix<double> &field, double dt, double &fen, double &stepsize2) {
	double fe=0.0;
	//cout << field << endl;
	matrix<double> den(Nr,no_types);
	if(polymer_or_fluid) {
	den = densityfromfield(field,fe);
	}
	else{
	den = densityfromfield_fluid(field,fe);	
	}
	// outfunc(den,"den");
	// outfunc(field,"field");
	double fen1  = this->freeenergy(den,field,fe);
	cout << "free energy: " << fen1 << endl;
	 
	fen = fen1;
	matrix<double> convolocal = this->genconvolocal(den);
	double stepsize = this->determine_stepsize2(field,convolocal,dt,fen1);
	cout << "steps: " << dt << " " << stepsize << endl;
	stepsize2 = stepsize;

	outfunc(field,"f1");
	outfunc(field,"t1");
	outfunc(convolocal,"t2");
	outfunc(den,"p1");
	//pausel();
	

	field = field+stepsize*(-field+convolocal);
}

void IsingPolymer::iteratefield_copolymer(matrix<double> &field, double dt, double &fen, double &stepsize2) {
	double fe=0.0;
	//cout << field << endl;
	matrix<double> den(Nr,no_types+1);
	vector1<double> field1(Nr);
	matrix<double> field2(Nr,no_types);

	for(int i = 0  ; i < Nr ; i++) {
		field1[i]=field(i,0);
		for(int j = 0  ; j < no_types ; j++) {
		field2(i,j)=field(i,j+1);
		}
	}	

	den = densityfromfield_copolymer(field1,field2,fe);

	//cout << "output fe: " << fe << endl;

	outfunc(den,"den");
	outfunc(field,"field");
	double fen1  = this->freeenergy_copolymer(den,field,fe);
	

	// cout <<  "fe1: " << this->freeenergy_copolymer(field) << endl;
	cout << "free energy: " << fen1 << endl;
	//pausel();

	fen = fen1;
	matrix<double> convolocal = this->genconvolocal_copolymer(den);
	outfunc(convolocal,"convolocal");
	
	double stepsize = this->determine_stepsize2_copolymer(field,convolocal,dt,fen1);
	cout << "steps: " << dt << " " << stepsize << endl;
	stepsize2 = stepsize;

	// outfunc(field,"f1");
	// outfunc((-stepsize)*field,"t1");

	// outfunc(stepsize*convolocal,"t2");
	// pausel();
	field = field+stepsize*(-field+convolocal);

}

void IsingPolymer::iteratefield_dcopolymer(matrix<double> &field, double dt, double &fen, double &stepsize2) {
	double fe=0.0;
	//cout << field << endl;
	matrix<double> den(Nr,no_types+no_types);
	matrix<double> field1(Nr,no_types);
	matrix<double> field2(Nr,no_types);

	for(int i = 0  ; i < Nr ; i++) {
		for(int j = 0  ; j < no_types ; j++) {
		field1(i,j)=field(i,j);
		}
		for(int j = 0  ; j < no_types ; j++) {
		field2(i,j)=field(i,j+no_types);
		}
	}	

	den = densityfromfield_dcopolymer(field1,field2,fe);

	//cout << "output fe: " << fe << endl;

	outfunc(den,"den");
	outfunc(field,"field");
	double fen1  = this->freeenergy_dcopolymer(den,field,fe);
	

	// cout <<  "fe1: " << this->freeenergy_copolymer(field) << endl;
	cout << "free energy: " << fen1 << endl;
	//pausel();

	fen = fen1;
	matrix<double> convolocal = this->genconvolocal_dcopolymer(den);
//	outfunc(convolocal,"convolocal");
	
	double stepsize = this->determine_stepsize2_dcopolymer(field,convolocal,dt,fen1);
	cout << "steps: " << dt << " " << stepsize << endl;
	stepsize2 = stepsize;

	// outfunc(field,"f1");
	// outfunc((-stepsize)*field,"t1");

	// outfunc(stepsize*convolocal,"t2");
	// pausel();
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
		if(abs(fi-ff)<0.01*stepsize2) { cout << "convergence test made" << endl; break; }
		else if(stepsize2<1E-7) {cout << "vanishing step" << endl; break; }
		else{ }
	}
	//cout << field << endl;
}


void IsingPolymer::run_copolymer(matrix<double> &field, int tot, double dt) {
	//cout << field << endl;
	double fen = 0.0;
	for(int i = 0 ; i < tot ; i++)
	{
		double fi = fen;
		cout << i << endl;
		double stepsize2 = 0.0;
		this->iteratefield_copolymer(field,dt,fen,stepsize2);

		double ff = fen;
		//cout << abs(fi-ff) << endl;
		if(abs(fi-ff)<stepsize2) { cout << "convergence test made" << endl; break; }
		else if(stepsize2<1E-7) {cout << "vanishing step" << endl; break; }
		else{ }
	}
	//cout << field << endl;
}

void IsingPolymer::run_dcopolymer(matrix<double> &field, int tot, double dt) {
	//cout << field << endl;
	double fen = 0.0;
	for(int i = 0 ; i < tot ; i++)
	{
		double fi = fen;
		cout << i << endl;
		double stepsize2 = 0.0;
		this->iteratefield_dcopolymer(field,dt,fen,stepsize2);

		double ff = fen;
		//cout << abs(fi-ff) << endl;
		if(abs(fi-ff)<stepsize2) { cout << "convergence test made" << endl; break; }
		else if(stepsize2<1E-7) {cout << "vanishing step" << endl; break; }
		else{ }
	}
	//cout << field << endl;
}

