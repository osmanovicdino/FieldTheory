LangevinNVT::LangevinNVT() {
gamma = 1.0;
dt = 1.0;
kT = 1.0;
m = 1.0;
vector1<bool> pb(1,true);
cube a(10.,pb,1);
this->setgeometry(a);
matrix<double> rs(1,a.dimension);
this->setdat(rs);

matrix<double> res(rs);

mom  = new matrix<double>;
*mom = res;


d = (gamma*dt/2.);
q = (dt)/2.;
r = sqrt(0.5*kT*(gamma)*(m)*(dt));

c1 = (dt/m);
c2 = (1.0/(1.0+(d)));
c3 = (1.0/(1.0+(d)))*q;
c4 = (1.0/(1.0+(d)))*r;
c5 = (1-(d));
}

LangevinNVT::LangevinNVT(geometry &a)  {
gamma = 1.0;
dt = 1.0;
kT = 1.0;
m = 1.0;
this->setgeometry(a);
matrix<double> rs(1,a.dimension);
this->setdat(rs);

matrix<double> res(rs);

mom  = new matrix<double>;
*mom = res;


d = (gamma*dt/2.);
q = (dt)/2.;
r = sqrt(0.5*kT*(gamma)*(m)*(dt));


c1 = (dt/m);
c2 = (1.0/(1.0+(d)));
c3 = (1.0/(1.0+(d)))*q;
c4 = (1.0/(1.0+(d)))*r;
c5 = (1-(d));
}




void LangevinNVT::setgamma(double g) {
	this->gamma = g;
	d = (gamma*dt/2.);
	q = (dt)/2.;
	r = sqrt(0.5*kT*(gamma)*(m)*(dt));

	c1 = (dt/m);
	c2 = (1.0/(1.0+(d)));
	c3 = (1.0/(1.0+(d)))*q;
	c4 = (1.0/(1.0+(d)))*r;
	c5 = (1-(d));	
}

void LangevinNVT::setdt(double d) {
	this->dt = d;
	d = (gamma*dt/2.);
	q = (dt)/2.;
	r = sqrt(0.5*kT*(gamma)*(m)*(dt));

	c1 = (dt/m);
	c2 = (1.0/(1.0+(d)));
	c3 = (1.0/(1.0+(d)))*q;
	c4 = (1.0/(1.0+(d)))*r;
	c5 = (1-(d));		
}

void LangevinNVT::setkT(double k) {
	this->kT = k;
	d = (gamma*dt/2.);
	q = (dt)/2.;
	r = sqrt(0.5*kT*(gamma)*(m)*(dt));

	c1 = (dt/m);
	c2 = (1.0/(1.0+(d)));
	c3 = (1.0/(1.0+(d)))*q;
	c4 = (1.0/(1.0+(d)))*r;
	c5 = (1-(d));		
}

void LangevinNVT::setm(double mm) {
	this->m = mm;
	d = (gamma*dt/2.);
	q = (dt)/2.;
	r = sqrt(0.5*kT*(gamma)*(m)*(dt));	

	c1 = (dt/m);
	c2 = (1.0/(1.0+(d)));
	c3 = (1.0/(1.0+(d)))*q;
	c4 = (1.0/(1.0+(d)))*r;
	c5 = (1-(d));		
}

void LangevinNVT::setmom(matrix<double> &ps) {
	if(ps.getNsafe() != this->getN() || ps.getncols() != this->getdimension() ) error("attempting to set momenta matrix which is not the same size as the system");
	mom = new matrix<double>(ps);
}

void LangevinNVT::printparams() {
	cout << "params are:" << endl;
	cout << "gamma" << gamma << endl;
	cout << "dt" << dt << endl;
	cout << "kT" << kT << endl;
	cout << "m" << m << endl;
	cout << "d" << d << endl;//((gamma*dt/2.));
	cout << "q" << q << endl;//((dt)/2.);
	cout << "r" << r << endl;//(sqrt(0.5*kT*(gamma)*(m)*(dt)));
	cout << endl;
}

matrix<double>& LangevinNVT::getmom() {
	return (*this->mom);
}


vector1<double> LangevinNVT::avmom() {
	int ds =  this->getdimension();
	vector1<double> momav(ds);
	for(int i = 0 ; i < (*mom).getNsafe() ; i++) {
		for(int i1 = 0 ; i1 < ds ; i1++ ) {
			momav[i1]+=((*mom)(i,i1)*(*mom)(i,i1))/(2*m);
		}
	}
	momav/=((double)(mom->getNsafe()));
	return momav;
}

void LangevinNVT::advance_pos() { 
	int ds = this->getdimension();
	//int locald = this->getdimension();
	for(int i = 0 ; i < (*dat).getNsafe() ;  i++ ) {
	for(int i1 = 0 ; i1 < ds ; i1++ ) {
	(*dat)(i,i1) = (*dat)(i,i1)+ c1*(*mom)(i,i1);
	//temp1->operator[](i1) = temp1->operator[](i1) + c1*(temp2->operator[](i1));
	}
	
	
	}
	(*this->geo).correct_position_and_momentum(*dat,*mom);

}

// double c1 = (dt/m);
// double c2 = (1.0/(1.0+(d)));
// double c3 = (1.0/(1.0+(d)))*q;
// double c4 = (1.0/(1.0+(d)))*r;
// double c5 = (1-(d));

void LangevinNVT::advance_mom(matrix<double> &F,matrix<double> &R){ 
	int ds = this->getdimension();
	// for(int i = 0 ; i < (*mom).getNsafe() ; i++) {
	// 	for(int i1 = 0; i1 < ds ; i1++ ) {
	// 		(mom)->operator()(i,i1) = c2*((mom)->operator()(i,i1)) + (c3)*F(i,i1) + (c4)*R(i,i1);
	// 	}
	// }
	for(int i = 0 ; i < (*mom).getNsafe() ; i++) {	
		for(int i1 = 0 ;  i1 < ds ; i1++ ) {
			//(mom)->operator()(i,i1) = c5*((mom)->operator()(i,i1))+(q)*F(i,i1)+(r)*R(i,i1);
			(mom)->operator()(i,i1) = c5*c2*((mom)->operator()(i,i1)) + (c5*(c3)+q)*F(i,i1) + (c5*(c4)+r)*R(i,i1);
		}
	}


}

void LangevinNVT::initialize(matrix<int> &pairs) {
	matrix<double> F(this->calculateforces(pairs));
	matrix<double> R(this->getN(),this->getdimension());
	int tempd = this->getdimension();
	for(int i = 0 ; i < this->getN() ; i++) {
		for(int j = 0 ; j < this->getdimension() ; j++) {
			R(i,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
		}
	}
	for(int i = 0 ; i < (*mom).getNsafe() ; i++) {	
		for(int i1 = 0 ;  i1 < tempd ; i1++ ) {
			(*mom)(i,i1) = (1-(d))*(*mom)(i,i1)+(q)*F(i,i1)+(r)*R(i,i1);
		}
	}
	for(int i = 0 ; i < (*dat).getNsafe() ;  i++ ) {
	for(int i1 = 0 ; i1 < tempd ; i1++ ) {
	(*dat)(i,i1) = (*dat)(i,i1) + (dt/m)*(*mom)(i,i1);
	}
	
	}
	(*this->geo).correct_position_and_momentum((*dat),(*mom));	
}

bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}



void LangevinNVT::adv(matrix<int> &pairs) {



matrix<double> F(this->calculateforces(pairs));


// int kj = 0;
// string filename;
// string filenamefirst = "F";
// stringstream ss;
// ss << kj;

// string extension = ".csv";
// filename = filenamefirst+ss.str();


// string filenamemomfirst = "mom";
// string posnamefirst = "pos";
// string rannamefirst = "Ran";

// string momfile = filenamemomfirst + ss.str();

// string posfile = posnamefirst + ss.str();

// string ranfile = rannamefirst + ss.str();


// while(fileExists(filename+extension)) {
// kj++;
// stringstream sstemp;
// sstemp << kj;
// filename = filenamefirst+sstemp.str();

// momfile = filenamemomfirst + sstemp.str();
// posfile = posnamefirst + sstemp.str();
// ranfile = rannamefirst + sstemp.str();
// }




// outfunc(F,filename);
//cout << F.max() << endl;
int dd = this->getdimension();
matrix<double> R(this->getN(),dd);
for(int i = 0 ; i < this->getN() ; i++) {
	for(int j = 0 ; j < dd ; j++) {
		R(i,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
	}
}
//outfunc(R,ranfile);
//cout << "advance mom commence" << endl;

this->advance_mom(F,R);

// cout << "advance mom complete" << endl;
// vector1<double> x(this->getN()),y(this->getN()),z(this->getN());
// for(int i = 0 ; i < this->getN() ; i++) {
// 	x[i]=(*mom)(i,0);
// 	y[i]=(*mom)(i,1);
// 	z[i]=(*mom)(i,2);
// }

// cout << maxval(x) << endl;
// cout << maxval(y) << endl;
// cout << maxval(z) << endl;
//outfunc((*mom),momfile);
this->advance_pos();
//outfunc((*dat),posfile);
//cout << "advance pos complete" << endl;
}