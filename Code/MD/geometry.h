struct geometry {
	int dimension;

	virtual double distance(const vector1<double>&,const vector1<double>&)=0; //bare distance

//	virtual bool distance_less_than(const vector1<double>&,const vector1<double>&,double)=0; //is the distance less than a certain value
	
//	virtual vector1<double> unit_vector(vector1<double>,vector1<double>,double&)=0; // return the unit vector between the two positions as well as the total distance

	virtual double distance(const matrix<double>&,int,int)=0; //bare distance between elements of a matrix

	virtual bool distance_less_than(matrix<double>&,int,int,double)=0; //is the distance less than a certain value between elements of matrix
	
//	virtual vector1<double> unit_vector(const matrix<double>&,int,int,double&)=0; // return the unit vector between the two positions as well as the total distance

//	virtual void unit_vector(const matrix<double>&,int,int,vector1<double>&,double&)=0; // return the unit vector between the two positions as well as the total distance
	virtual void distance3v(double&,double&,double&,double&,double&,double&,vector1<double> &un, double&)=0;

	virtual void distance_vector(matrix<double>&,int&,int&,vector1<double>&,double&)=0; // return the unit vector between the two positions as well as the total distance

	//virtual void correct_position(vector1<double>&)=0;
	virtual bool exceed_dis(matrix<double>&,matrix<double>&,double)=0;

	virtual void correct_position(matrix<double>&)=0; //using the periodic boundary conditions, correct the position and momentum 

	virtual void correct_position_and_momentum(matrix<double>&,matrix<double>&)=0;

	virtual int assign_box(matrix<double>&,int,vector1<int> &, int )=0;

	virtual matrix<int> generate_boxes_relationships(int,int&)=0;

	//virtual bool geometrytest(vector1<double>&) = 0; // is the vector1 within the simulation box ?

	virtual geometry* clone() const = 0;
};

struct cube : geometry {
	double l;

	vector1<bool> pb;

	bool periodic; //Kluge
	double l2s;

	cube(double ll, vector1<bool> pbb, int dim) : l(ll),pb(pbb) {
		dimension = dim;
		periodic = pbb[0];
		l2s = SQR(ll/2.);
		if(pbb.getsize() != dim) error("warning: specification of boundary conditions must be the same dimension as the system");
	
	}

	double distance(const vector1<double> &x1,const vector1<double> &x2) {
		double temp = 0.0;
		//vector1<double> dx(dimension);
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = x1.gpcons(i1)-x2.gpcons(i1);
			if(pb[i1]) { 
				if(abs(dx) > l*0.5 ) 
					dx = dx - SIGN(l,dx); 
			}
			
				temp += SQR(dx);
			
		}

		return sqrt(temp);	

	}

	double distance(const matrix<double> &r,int i, int j) {
		double temp = 0.0;
		//vector1<double> dx(dimension);
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = r.gpcons(i,i1)-r.gpcons(j,i1);
			if(pb[i1]) { 
				if(abs(dx) > l*0.5 ) 
					dx = dx - SIGN(l,dx); 
			}
			
				temp += SQR(dx);
			
		}

		return sqrt(temp);	

	}


	// bool distance_less_than(const vector1<double> &x1,const vector1<double> &x2, double val) {
	// 	double temp = 0.0;
	// 	//vector1<double> dx(dimension);
	// 	for(int i1  = 0; i1 < dimension ; i1++ ) {
	// 		double dx = x1.gpcons(i1)-x2.gpcons(i1);
	// 		if(pb[i1]) { 
	// 			if(abs(dx) > l*0.5 ) 
	// 				dx = dx - SIGN(l,dx); 
	// 		}
	// 		if(abs(dx)>val) {return false; }
	// 		temp += SQR(dx);
			
	// 	}

	// 	if(temp < SQR(val)) return true;
	// 	else return false;	

	// }

	bool distance_less_than(matrix<double> &r,int i, int j, double val) {
		double temp = 0.0;
		//vector1<double> dx(dimension);
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = r.mat[i*dimension+i1]-r.mat[j*dimension+i1];

			if(periodic) { 
				if(fabs(dx) > l*0.5 ) 
					dx = dx - SIGN(l,dx); 
			}
			dx = SQR(dx);
			if(dx>SQR(val)) {return false; }
			temp += dx;
			
		}

		if(temp < SQR(val)) return true;
		else return false;	

	}

	bool exceed_dis(matrix<double> &r1, matrix<double> &r2, double d) { //returns the square distance to d and the distance between to uv
		//vector1<double> dx(dimension);
		
		//vector1<double> dx(dimension);
		for(int i = 0  ; i < r1.getNsafe() ; i++) {
			double temp = 0.0;
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = r1.mat[i*dimension+i1]-r2.mat[i*dimension+i1];

			if(periodic) { 
				if(fabs(dx) > l*0.5 ) 
					dx = dx - SIGN(l,dx); 
			}
			dx = SQR(dx);
			if(dx>d) {return true; }
			temp += dx;
			
		}

		if(temp > d) return true;	
		}
		return false;
		//uv*=(1./d);
		// for(int i1 = 0 ; i1 < dimension ; i1++) {
		// 	dx[i1]=d2*dx[i1];
		// }
		//return dx;	
	}		

	// vector1<double> unit_vector(vector1<double> x1, vector1<double> x2, double &d) {
	// 	vector1<double> dx(dimension);
	// 	double temp = 0.0;
	// 	for(int i1  = 0; i1 < dimension ; i1++ ) {
	// 		double as = x1[i1]-x2[i1];
	// 		if(pb[i1]) { if(fabs(as) > l*0.5 ) 
	// 			as = as - SIGN(l,as); }
	// 		temp += SQR(as);
	// 		dx[i1] = as;
	// 	}
	// 	d = sqrt(temp);
	// 	dx*=(1./d);
	// 	// for(int i1 = 0 ; i1 < dimension ; i1++) {
	// 	// 	dx[i1]=d2*dx[i1];
	// 	// }
	// 	return dx;		
	// }

	// vector1<double> unit_vector(const matrix<double> &r,int i, int j, double &d) {
	// 	vector1<double> dx(dimension);
	// 	double temp = 0.0;
	// 	for(int i1  = 0; i1 < dimension ; i1++ ) {
	// 		double as = r.gpcons(i,i1)-r.gpcons(j,i1);
	// 		if(pb[i1]) { if(fabs(as) > l*0.5 ) 
	// 			as = as - SIGN(l,as); }
	// 		temp += SQR(as);
	// 		dx[i1] = as;
	// 	}
	// 	d = sqrt(temp);
	// 	dx*=(1./d);
	// 	// for(int i1 = 0 ; i1 < dimension ; i1++) {
	// 	// 	dx[i1]=d2*dx[i1];
	// 	// }
	// 	return dx;	
	// }


	// void unit_vector(const matrix<double> &r,int i, int j, vector1<double> &uv, double &d) {
	// 	//vector1<double> dx(dimension);
	// 	double temp = 0.0;
	// 	for(int i1  = 0; i1 < dimension ; i1++ ) {
	// 		double as = r.gpcons(i,i1)-r.gpcons(j,i1);
	// 		if(pb[i1]) { if(fabs(as) > l*0.5 ) 
	// 			as = as - SIGN(l,as); }
	// 		temp += SQR(as);
	// 		uv[i1] = as;
	// 	}
	// 	d = sqrt(temp);
	// 	uv*=(1./d);
	// 	// for(int i1 = 0 ; i1 < dimension ; i1++) {
	// 	// 	dx[i1]=d2*dx[i1];
	// 	// }
	// 	//return dx;	
	// }

	void distance3v(double &x1, double &x2, double &y1, double &y2, double &z1, double &z2, vector1<double> &uv, double &d) { //returns the square distance to d and the distance between to uv
		//vector1<double> dx(dimension);
		
		
		double dx = x1-x2;
		double dy = y1-y2;
		double dz = z1-z2;
		if(periodic) { if(SQR(dx) > l2s ) 
			dx = dx - SIGN(l,dx);
			if(SQR(dy) > l2s ) 
			dy = dy - SIGN(l,dy);
			if(SQR(dz) > l2s ) 
			dz = dz - SIGN(l,dz);
			}
		d = SQR(dx)+SQR(dy)+SQR(dz);
		uv[0] = dx;
		uv[1] = dy;
		uv[2] = dz;
		
		
		//uv*=(1./d);
		// for(int i1 = 0 ; i1 < dimension ; i1++) {
		// 	dx[i1]=d2*dx[i1];
		// }
		//return dx;	
	}	

	void distance_vector(matrix<double> &r, int &i, int &j, vector1<double> &uv, double &d) { //returns the square distance to d and the distance between to uv
		//vector1<double> dx(dimension);
		double temp = 0.0;
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double as = r.mat[i*dimension+i1]-r.mat[j*dimension+i1];
			if(periodic) { if(SQR(as) > l2s ) 
				as = as - SIGN(l,as); }
			temp += SQR(as);
			uv[i1] = as;
		}
		d = temp;
		//uv*=(1./d);
		// for(int i1 = 0 ; i1 < dimension ; i1++) {
		// 	dx[i1]=d2*dx[i1];
		// }
		//return dx;	
	}	



	// void correct_position(vector1<double> &r) {
	// 		for(int i = 0 ; i < dimension ; i++) {
	// 			if(pb[i]) {
	// 				if(r[i] < 0) {	
	// 					while(r[i] < 0) {
	// 						r[i] += l;
	// 					}
	// 				}
	// 				else if(r[i] > l) {
	// 					while(r[i] > l) {
	// 						r[i] -= l;
	// 					}
	// 				}
	// 				else {
	// 			//do nothing
	// 				}
	// 			}
	// 			else {
	// 				if(r[i]<0) {
	// 					r[i] = -r[i];
	// 				}
	// 				else if(r[i]>l) {
	// 					r[i] = l-(r[i]-l);
	// 				}
	// 				else {
	// 				//do nothing
	// 				}
	// 			}
	// 		}
		
	// }

	void correct_position(matrix<double> &r) {
		for(int j = 0 ; j < r.getNsafe() ; j++) {
			for(int i = 0 ; i < dimension ; i++) {
				if(pb[i]) {
					if(r(j,i) < 0) {	
						while(r(j,i) < 0) {
							r(j,i) += l;
						}
					}
					else if(r(j,i) > l) {
						while(r(j,i) > l) {
							r(j,i) -= l;
						}
					}
					else {
				//do nothing
					}
				}
				else {
					if(r(j,i)<0) {
						r(j,i) = -r(j,i);
					}
					else if(r(j,i)>l) {
						r(j,i) = l-(r(j,i)-l);
					}
					else {
					//do nothing
					}
				}
			}
		}
	}

	void correct_position_and_momentum(matrix<double> &r, matrix<double> &p) {
		//for(int j = 0 ; j < r.getNsafe() ; j++) {
			for(int i = 0 ; i < dimension ; i++) {
				if(pb[i]) {
					for(int j = 0 ; j < r.getNsafe() ; j++) {
						double temp = r(j,i);
						if(temp < 0) {	
							while(temp < 0) {
								temp += l;
							}
							r(j,i) = temp;
						}
						else if(temp > l) {
							while(temp > l) {
								temp -= l;
							}
							r(j,i) = temp;
						}
						else {
					//do nothing
						}
					}
				}
				else {
					for(int j = 0 ; j < r.getNsafe() ; j++) {
						if(r(j,i)<0) {
							r(j,i) = -r(j,i);
							p(j,i) = -p(j,i);
						}
						else if(r(j,i)>l) {
							r(j,i) = l-(r(j,i)-l);
							p(j,i) = -p(j,i);
						}
						else {
						//do nothing
						}
					}
				}
			}
		}
	

	int assign_box(matrix<double> &v1,int i,vector1<int> &dim, int n) {
		vector1<int> cs(dimension);
		for(int i1 = 0; i1 < dimension ; i1++ ) {
			int m = fasterfloor(n*(v1(i,i1))/l);
			if(m<0) m=0;
			else if(m>n-1) m=n-1;
			cs[i1]=m;
		}

		int c = scalar(dim,cs);
		return c;
	}

	matrix<int> generate_boxes_relationships(int n,int &tb) {
		tb=pow(n,dimension);
		double box_size = l/(double)n;
		//for(int i = 0 ; i < n ; i++)
		matrix<int> boxcheckk(tb,pow(3,dimension));
		vector1<int> cvec(dimension);
		for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= n;
		}
		cvec[i] = ij;
		}

		int iter  = 0;
		vector1<int> iterator(dimension,0);

		for(;;) {
		//cout << iterator << endl;
		for(int j = 0 ; j < dimension-1 ; j++ ) {
		iterator[j+1]=floor(iterator[j]/n);
		}

		vector1<int> iterator2(dimension);
		for(int j = 0 ; j < dimension ; j++ ) {
		iterator2[j] = (iterator[j])%(n);
		}


		int c = scalar(iterator2,cvec);

		int j= 0;

		vector1<int> temp(dimension,0);
		for(;;) {
			for(int j1 = 0 ; j1 < dimension-1 ; j1++ ) {
			temp[j1+1]=floor(temp[j1]/3);
			}
		vector1<int> temp2(dimension);
		for(int j1 = 0 ; j1 < dimension ; j1++ ) {
		temp2[j1]=(temp[j1]%3)-1;
		}

		vector1<int> cint(dimension);
		for(int j1 = 0 ; j1 < dimension ; j1++ ) {
		cint[j1] = iterator2[j1]+temp2[j1];
		if(cint[j1]<0) { cint[j1] = cint[j1] + n; }
		else if(cint[j1] > n-1) {cint[j1] = cint[j1] -n; }
		else { }

		}
		int c2  = scalar(cint,cvec);
		boxcheckk(c,temp[0])=c2;

		temp[0] = temp[0]+1;
		if(temp[0] > pow(3,dimension)-1 ) break;
		}
		iterator[0]=iterator[0]+1;
		if(c>=tb-1) break;
		}

		vector<int> chckvec;
		for(int i = 0 ; i < boxcheckk.getncols() ; i++) {
			chckvec.push_back(boxcheckk(0,i));
		}

		sort(chckvec.begin(),chckvec.end());
		chckvec.erase(unique(chckvec.begin(), chckvec.end()), chckvec.end());

		matrix<int> boxcheckk2(boxcheckk.getNsafe(),chckvec.size());

		for(int i =0 ; i < boxcheckk2.getNsafe() ; i++) {
			vector<int> chckvec2;
			for(int i1 = 0 ; i1 < boxcheckk.getncols() ; i1++) {
				chckvec2.push_back(boxcheckk(i,i1));
			}
			sort(chckvec2.begin(),chckvec2.end());
			chckvec2.erase(unique(chckvec2.begin(), chckvec2.end()), chckvec2.end());	

			for(int j = 0  ; j < chckvec2.size() ; j++ ) {
				boxcheckk2(i,j)=chckvec2[j];
			} 		
		}

		matrix<int> res(boxcheckk2);
		return res;
	}


	// bool geometrytest(vector1<double> &r) {
	// 	bool geo = true;
	// 	for(int i = 0 ; i < dimension ; i++) {
	// 		if(!pb[i] && (r[i] > l || r[i] < 0)) {
	// 			geo = false;
	// 			break;
	// 		}
	// 	}
	// 	return geo;
	// }

	cube* clone() const {
		return new cube(*this);
	}

};

/*
struct cuboid : geometry {
	vector1<double> l;

	vector1<bool> pb;


	cuboid(vector1<double> ll, vector1<bool> pbb, int dim) : l(ll),pb(pbb) {
		dimension = dim;
		if(pbb.getsize() != dim) error("warning: specification of boundary conditions must be the same dimension as the system");
	
	}

	double distance(const vector1<double> &x1,const vector1<double> &x2) {
		double temp = 0.0;
		//vector1<double> dx(dimension);
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = x1.gpcons(i1)-x2.gpcons(i1);
			if(pb[i1]) { 
				if(abs(dx) > l[i1]*0.5 ) 
					dx = dx - SIGN(l[i1],dx); 
			}
			
				temp += SQR(dx);
			
		}

		return sqrt(temp);	

	}

	bool distance_less_than(const vector1<double> &x1,const vector1<double> &x2, double val) {
		double temp = 0.0;
		//vector1<double> dx(dimension);
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			double dx = x1.gpcons(i1)-x2.gpcons(i1);
			if(pb[i1]) { 
				if(abs(dx) > l[i1]*0.5 ) 
					dx = dx - SIGN(l[i1],dx); 
			}
			if(abs(dx)>val) {return false; }
			temp += SQR(dx);
			
		}

		if(temp < SQR(val)) return true;
		else return false;	

	}

	vector1<double> unit_vector(vector1<double> x1, vector1<double> x2, double &d) {
		vector1<double> dx(dimension);
		double temp = 0.0;
		for(int i1  = 0; i1 < dimension ; i1++ ) {
			dx[i1] = x1[i1]-x2[i1];
			if(pb[i1]) { if(abs(dx[i1]) > l[i1]*0.5 ) 
				dx[i1] = dx[i1] - SIGN(l[i1],dx[i1]); }
			temp += SQR(dx[i1]);
		}
		d = sqrt(temp);
		dx /= d;
		return dx;		
	}
	void correct_position(vector1<double> &r) {
			for(int i = 0 ; i < dimension ; i++) {
				if(pb[i]) {
					if(r[i] < 0) {	
						while(r[i] < 0) {
							r[i] += l[i];
						}
					}
					else if(r[i] > l[i]) {
						while(r[i]> l[i]) {
							r[i] -= l[i];
						}
					}
					else {
				//do nothing
					}
				}
				else {
					if(r[i]<0) {
						r[i] = -r[i];
					}
					else if(r[i]>l[i]) {
						r[i] = l[i]-(r[i]-l[i]);
					}
					else {
					//do nothing
					}
				}
			}
	}

	void correct_position(matrix<double> &r) {
		for(int j = 0 ; j <  r.getNsafe() ; j++ ) {
			for(int i = 0 ; i < dimension ; i++) {
				if(pb[i]) {
					if(r(j,i) < 0) {	
						while(r(j,i) < 0) {
							r(j,i) += l[i];
						}
					}
					else if(r(j,i) > l[i]) {
						while(r(j,i)> l[i]) {
							r(j,i) -= l[i];
						}
					}
					else {
				//do nothing
					}
				}
				else {
					if(r(j,i)<0) {
						r(j,i) = -r(j,i);
					}
					else if(r(j,i)>l[i]) {
						r(j,i) = l[i]-(r(j,i)-l[i]);
					}
					else {
					//do nothing
					}
				}
			}
		}
	}

	void correct_position_and_momentum(matrix<double> &r, matrix<double> &p) {
		for(int j = 0  ;  j < r.getNsafe() ; j++ ) {
			for(int i = 0 ; i < dimension ; i++) {
				if(pb[i]) {
					if(r(j,i) < 0) {	
						while(r(j,i) < 0) {
							r(j,i) += l[i];
						}
					}
					else if(r(j,i) > l[i]) {
						while(r(j,i) > l[i]) {
							r(j,i) -= l[i];
						}
					}
					else {
				//do nothing
					}
				}
				else {
					if(r(j,i)<0) {
						r(j,i) = -r(j,i);
						p(j,i) = -p(j,i);
					}
					else if(r(j,i)>l[i]) {
						r(j,i) = l[i]-(r(j,i)-l[i]);
						p(j,i) = -p(j,i);
					}
					else {
					//do nothing
					}
				}
			}
		}
	}

	int assign_box(vector1<double> v1,vector1<int> &dim, int n) {
		return 0;
	}

	matrix<int> generate_boxes_relationships(int c, int &a) {
		return matrix<int>(1,3);
	}

	bool geometrytest(vector1<double> &r) {
		bool geo = true;
		for(int i = 0 ; i < dimension ; i++) {
			if(!pb[i] && (r[i] > l[i] || r[i] < 0)) {
				geo = false;
				break;
			}
		}
		return geo;
	}

	cuboid* clone() const {
		return new cuboid(*this);
	}
};

// */

// struct freespace : geometry {

// 	freespace(int dim) { dimension = dim; }

// 	double distance(const vector1<double> &x1,const vector1<double> &x2) {
// 		double temp = 0.0;
// 		//vector1<double> dx(dimension);
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			double dx = x1.gpcons(i1)-x2.gpcons(i1);
// 			temp += SQR(dx);
// 		}
// 		return sqrt(temp);	
// 	}

// 	double distance(const matrix<double> &r,int i, int j) {
// 		double temp = 0.0;
// 		//vector1<double> dx(dimension);
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			double dx = r.gpcons(i,i1)-r.gpcons(j,i1);
// 			temp += SQR(dx);
// 		}
// 		return sqrt(temp);	
// 	}

// 	bool distance_less_than(const vector1<double> &x1,const vector1<double> &x2, double val) {
// 		double temp = 0.0;
// 		//vector1<double> dx(dimension);
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			double dx = x1.gpcons(i1)-x2.gpcons(i1);
// 			if(fabs(dx)>val) return false;
// 			temp += SQR(dx);
// 		}
// 		if(temp < SQR(val)) return true;
// 		else return false;			

// 	}
// 	bool distance_less_than(matrix<double> &r,int i, int j, double val) {
// 		double temp = 0.0;
// 		//vector1<double> dx(dimension);
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			double dx = r(i,i1)-r(j,i1);
// 			dx = SQR(dx);
// 			if(dx>SQR(val)) return false;
// 			temp += dx;
// 		}
// 		if(temp < SQR(val)) return true;
// 		else return false;			

// 	}

// 	vector1<double> unit_vector(vector1<double> x1, vector1<double> x2, double &d) {
// 		vector1<double> dx(dimension);
// 		double temp = 0.0;
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			dx[i1] = x1[i1]-x2[i1];
// 			temp += SQR(dx[i1]);
// 		}
// 		d = sqrt(temp);
// 		dx /= d;
// 		return dx;		
// 	}

// 	vector1<double> unit_vector(const matrix<double> &r,int i, int j, double &d) {
// 		vector1<double> dx(dimension);
// 		double temp = 0.0;
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			dx[i1] = r.gpcons(i,i1)-r.gpcons(j,i1);
// 			temp += SQR(dx[i1]);
// 		}
// 		d = sqrt(temp);
// 		dx /= d;
// 		return dx;		
// 	}

// 	void unit_vector(const matrix<double> &r,int i, int j, vector1<double> &dx, double &d) {
// 		//vector1<double> dx(dimension);
// 		double temp = 0.0;
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			dx[i1] = r.gpcons(i,i1)-r.gpcons(j,i1);
// 			temp += SQR(dx[i1]);
// 		}
// 		d = sqrt(temp);
// 		dx /= d;
// 		//return dx;		
// 	}

// 	void distance_vector(matrix<double> &r,int i, int j, vector1<double> &dx, double &d) {
// 		//vector1<double> dx(dimension);
// 		double temp = 0.0;
// 		for(int i1  = 0; i1 < dimension ; i1++ ) {
// 			dx[i1] = r(i,i1)-r(j,i1);
// 			temp += SQR(dx[i1]);
// 		}
// 		d = sqrt(temp);
// 		//dx /= d;
// 		//return dx;		
// 	}


// 	void correct_position(vector1<double> &r) {
// 		//do nothing
// 	}

// 	void correct_position(matrix<double> &r) {
// 		// free space
// 	}	
// 	void correct_position_and_momentum(matrix<double> &r, matrix<double> &p) {
// 		// free space
// 	}	

// 	int assign_box(vector1<double> v1,vector1<int> &dim, int n) {
// 		return 0;
// 	}

// 	matrix<int> generate_boxes_relationships(int c, int &a) {
// 		return matrix<int>(1,3);
// 	}


// 	bool geometrytest(vector1<double> &r) {
// 		return true;
// 	}

// 	freespace* clone() const {
// 		return new freespace(*this);
// 	}	
// };
