

// vector1<double> NCGas::forcecalc1(double x1, double x2, double y1, double y2, double z1, double z2) {


// 	double dx = x2 - x1;
// 	double dy = y2 - y1;
// 	double dz = z2 - z1;
// 	double d = SQR(dx) + SQR(dy) + SQR(dz);
// 	vector1<double> force(3);
// 	force[0]=(-2*eps*(dy))/(exp((d)/lambda)*lambda) + (4*eps*(dz))/(exp((2*(SQR(dx)+SQR(dy) + SQR(dz)))/lambda)*lambda);
// 	force[1]=(2*eps*(dx))/(exp((d)/lambda)*lambda) - (2*eps*(dz))/(exp((d)/lambda)*lambda);
// 	force[2]=(-4*eps*(dx))/(exp((2*(d))/lambda)*lambda) + (2*eps*(dy))/(exp((d)/lambda)*lambda);
// 	return force;
// }

// vector1<double> NCGas::forcecalc2(int i, int j) {
// 	double dis = 0.0;
// 	vector1<double> uv = obj->unit_vector(i,j,dis);
// 	double dx = uv.gpcons(0)*dis;
// 	double dy = uv.gpcons(1)*dis;
// 	double dz = uv.gpcons(2)*dis;
// 	//double d = SQR(dis);
// 	vector1<double> force(3);
// 	force[0]=((-dy + dz)*eps)/(dis*exp(dis/lambda)*lambda);
// 	force[1]=((dx - dz)*eps)/(dis*exp(dis/lambda)*lambda);
// 	force[2]=((-dx + dy)*eps)/(dis*exp(dis/lambda)*lambda);
// 	return force;
// }

// vector1<double> NCGas::forcecalc2(int i, int j) {
// 	double dis = 0.0;
// 	vector1<double> uv = obj->unit_vector(i,j,dis);
// 	double dx = uv.gpcons(0)*dis;
// 	double dy = uv.gpcons(1)*dis;
// 	double dz = uv.gpcons(2)*dis;
// 	//double d = SQR(dis);
// 	vector1<double> force(3);
// 	force[0]=-((dx*SQR(dy)*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda)) + (dx*SQR(dz)*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda);
// 	force[1]=(SQR(dx)*dy*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda) - (dy*SQR(dz)*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda);
// 	force[2]=-((SQR(dx)*dz*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda)) + (SQR(dy)*dz*exp(-SQR(dis)/(2.*SQR(lambda))))/SQR(lambda);
// 	return (eps/SQR(lambda))*force;
// }

// vector1<double> NCGas::forcecalc2(int i, int j, vector1<double> &force) {
// 	double dis = 0.0;
// 	//vector1<double> uv = obj->unit_vector(i,j,dis);
// 	//(obj->getcoordinate(i,0))-(obj->getcoordinate(j,0))
// 	double dx = (obj->getcoordinate(i,0))-(obj->getcoordinate(j,0));

// 	if(fabs(dx) > chm) dx = dx - SIGN(l,dx);

// 	double dy = (obj->getcoordinate(i,1))-(obj->getcoordinate(j,1));
// 	if(fabs(dy) > chm) dy = dy - SIGN(l,dy);

// 	double dz = (obj->getcoordinate(i,2))-(obj->getcoordinate(j,2));
// 	if(fabs(dz) > chm) dz = dz - SIGN(l,dz);

// 	double d = SQR(dx)+SQR(dy)+SQR(dz);

// 	force[0]=((-dy + dz)*eps)/(d*exp(d/lambda)*lambda);
// 	force[1]=((dx - dz)*eps)/(d*exp(d/lambda)*lambda);
// 	force[2]=((-dx + dy)*eps)/(d*exp(d/lambda)*lambda);
	
// }

// void NCGas::forcecalc2(int i, int j, vector1<double> &force) {
// 	double dis = 0.0;
// 	//vector1<double> uv = obj->unit_vector(i,j,dis);
// 	//(obj->getcoordinate(i,0))-(obj->getcoordinate(j,0))
// 	double dx = (obj->getcoordinate(i,0))-(obj->getcoordinate(j,0));

// 	if(fabs(dx) > chm) dx = dx - SIGN(l,dx);

// 	double dy = (obj->getcoordinate(i,1))-(obj->getcoordinate(j,1));
// 	if(fabs(dy) > chm) dy = dy - SIGN(l,dy);

// 	double dz = (obj->getcoordinate(i,2))-(obj->getcoordinate(j,2));
// 	if(fabs(dz) > chm) dz = dz - SIGN(l,dz);

// 	double d = SQR(dx)+SQR(dy)+SQR(dz);



// 	force[0]=((-dy + dz)*eps)/(exp(d/lambda)*lambda);
// 	force[1]=((dx - dz)*eps)/(exp(d/lambda)*lambda);
// 	force[2]=((-dx + dy)*eps)/(exp(d/lambda)*lambda);


// }

void NCGas::forcecalc2(int i, int j, vector1<double> &force) {
	double dis = 0.0;
	//vector1<double> uv = obj->unit_vector(i,j,dis);
	//(obj->getcoordinate(i,0))-(obj->getcoordinate(j,0))
	double dx = (obj->getcoordinate(i,0))-(obj->getcoordinate(j,0));

	if(fabs(dx) > chm) dx = dx - SIGN(l,dx);

	double dy = (obj->getcoordinate(i,1))-(obj->getcoordinate(j,1));
	if(fabs(dy) > chm) dy = dy - SIGN(l,dy);

	double dz = (obj->getcoordinate(i,2))-(obj->getcoordinate(j,2));
	if(fabs(dz) > chm) dz = dz - SIGN(l,dz);

	double d = SQR(dx)+SQR(dy)+SQR(dz);

	double fac =  2.*eps*(1./exp(d/SQR(lambda)));


	force[0]=fac*-((SQR(dx)*dy + dz*(-dz + lambda)*(dz + lambda))/lam4);
	force[1]=fac*(SQR(dx)*dx - SQR(dy)*dz - dx*lam2)/lam4;
	force[2]=fac*-((-SQR(dy)*dy + dx*SQR(dz) + dy*lam2)/lam4);

}

// void NCGas::forcecalc2(int i, int j, vector1<double> &force) { //periodic assumption
// 	double dis = 0.0;
// 	//vector1<double> uv = obj->unit_vector(i,j,dis);
// 	//(obj->getcoordinate(i,0))-(obj->getcoordinate(j,0))
// 	double dx = (obj->getcoordinate(i,0))-(obj->getcoordinate(j,0));

// 	if(fabs(dx) > chm) dx = dx - SIGN(l,dx);

// 	double dy = (obj->getcoordinate(i,1))-(obj->getcoordinate(j,1));
// 	if(fabs(dy) > chm) dy = dy - SIGN(l,dy);

// 	double dz = (obj->getcoordinate(i,2))-(obj->getcoordinate(j,2));
// 	if(fabs(dz) > chm) dz = dz - SIGN(l,dz);

// 	double dis2 = SQR(dx)+SQR(dy)+SQR(dz);
// 	//double d = SQR(dis);
// 	//vector1<double> force(3);
// 	force[0]=eps*(-((dx*SQR(dy)*exp(-(dis2)/(2.*(lam2))))/(lam4)) + (dx*SQR(dz)*exp(-(dis2)/(2.*(lam2))))/(lam4));
// 	force[1]=eps*((SQR(dx)*dy*exp(-(dis2)/(2.*(lam2))))/(lam4) - (dy*SQR(dz)*exp(-(dis2)/(2.*(lam2))))/(lam4));
// 	force[2]=eps*(-((SQR(dx)*dz*exp(-(dis2)/(2.*(lam2))))/(lam4)) + (SQR(dy)*dz*exp(-(dis2)/(2.*(lam2))))/(lam4));
// 	//return force;
// }

// vector1<double> NCGas::forcecalc3(int i) {
// 	// double dis = 0.0;
// 	// vector1<double> uv = obj->unit_vector(i,j,dis);	
// 	double x = obj->getcoordinate(i,0);
// 	double y = obj->getcoordinate(i,1);
// 	double z = obj->getcoordinate(i,2);
// 	// double d = SQR(dis);
// 	vector1<double> force(3);
// 	force[0]=-((exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + y))/SQR(lambda)) + (exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + z))/SQR(lambda);
// 	force[1]=(exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + x))/SQR(lambda) - (exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + z))/SQR(lambda);
// 	force[2]=-((exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + x))/SQR(lambda)) + (exp((-SQR(-6.02 + x) - SQR(-6.02 + y) - SQR(-6.02 + z))/(2.*SQR(lambda)))*eps*(-6.02 + y))/SQR(lambda);
// 	return force;
// }


NCGas::NCGas(double ll, int totalN) {
	obj = new LangevinNVT;

	l = ll;
	num = floor(ll/4.0);
	double pi = 3.1415;
	eps = 1.0;
	eqeps = 1.0;
	lambda = 1.5;

	chm = 0.5*l;
	lam2 = SQR(lambda);
	lam4 = SQR(lam2);

	vector1<bool> pb(3,true);
	cube bc(ll,pb,3);

	
	matrix<double> store(totalN,3);



//	intmatrix f(totalN);
	double sigma  = 1.0;
	WCAPotential nswca(1.0, sigma, eqeps);

	//double specific_volume  = sigma*sigma*sigma;
	//double dls = sigma;
/*	int pp = floor(l/sigma);
	vector<vector<double> > possible_pos;
	possible_pos.reserve(pp*pp*pp);
	for(int i = 0 ; i < pp ; i++)
		for(int j = 0 ; j < pp ; j++)
			for(int k = 0  ; k < pp ; k++) {
				double x =  0.5*sigma+i*sigma;
				double y =  0.5*sigma+j*sigma;
				double z =  0.5*sigma+k*sigma;
				vector<double> b;
				b.push_back(x);
				b.push_back(y);
				b.push_back(z);

				possible_pos.push_back(b);
			}

			cout << "system possible set up complete" << endl;

	for(int i = 0 ; i < totalN ; i++) {
		cout << i << endl;
		int randint = rand() % (possible_pos.size());
		store(i,0) =  possible_pos[randint][0];
		store(i,1) =  possible_pos[randint][1];
		store(i,2) =  possible_pos[randint][2];
		possible_pos.erase(possible_pos.begin()+randint);

	}	

	cout << "system actual set up complete" << endl;*/

	//cout << store << endl;


	//LJPotential nswca(1.0,1.0);
	//ExponentialPotential nswca(eqeps,1.0);
	
//	f.add_interaction_all(nswca);

	LangevinNVT b(bc);


	double kT = 1.0;
	double dt = 0.005;
	double eta = 50.;

	// int sweeps = 5;
	// double max_step = 0.2;

	// s_matrix<double> tempdat(store);

	// MC_Metropolis a(bc);
	// a.setdat(tempdat);
	// a.set_sweeps(sweeps);
	// a.setkT(kT);
	// a.set_max_stepsize(max_step);
	// a.setinteractions(f);

	// matrix<int> *yolo = a.calculatepairs();
	// cout << "pairs calculated" << endl;

	// for(int i = 0 ; i < 1000*a.getN() ; i++) {	
		
	// 	a.adv(*yolo);

	// 	if((*yolo)(yolo->getNsafe()-1,0)==1) {
	// 		delete yolo;
	// 		cout << 100.0*(double)i/(1000.*a.getN()) << "%\n"; 
	// 		yolo = a.calculatepairs();
	// 		}
	// }


	// delete yolo;
	cout << "STARTING IMPORT" << endl;
	double T;
	bool err1;
	matrix<double> stat = importcsv("/home/dino/Documents/Demons/Data/exampledata.csv",T,err1);
	cout << "imported" << endl;
	if(err1) error("failed to impport");

	// double alpha = 2*pi*rand();
	// double beta = 2*pi*rand();
	// double gamma = 2*pi*rand();

	// matrix<double> rot(3,3);

	// rot(0,0)=cos(alpha)*cos(beta);
	// rot(0,1)=cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma);
	// rot(0,2)=cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
	// rot(1,0)=sin(alpha)*cos(beta);
	// rot(1,1)=sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma);
	// rot(1,2)=sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
	// rot(2,0)=-sin(beta);
	// rot(2,1)=cos(beta)*sin(gamma);
	// rot(2,2)=cos(beta)*cos(gamma);

	// for(int i = 0 ; i < stat.getNsafe() ; i++ ) {
	// 	vector1<double> pos(3);
	// 	pos[0]=stat(i,0)-80.572/2.;
	// 	pos[1]=stat(i,1)-80.572/2.;
	// 	pos[2]=stat(i,2)-80.572/2.;

	// 	vector1<double> pos2 = rot*pos;
	// 	stat(i,0)=pos2[0];
	// 	stat(i,1)=pos2[1];
	// 	stat(i,2)=pos2[2];
	// }

	// bc.correct_position(stat);

	//apply a random rotation

	b.setdat(stat);
	
	//b.setdat(store);
	b.setinteractions(nswca);
	b.setkT(kT);
	b.setdt(dt);
	b.setgamma(eta);
	b.setm(1.0);
	matrix<double> moms(b.getN(),3);
	b.setmom(moms);


	*obj  = b;



};



void NCGas::seteps(double epss) {
	eps = epss;
}

void NCGas::setkT(double kT) {
	obj->setkT(kT);
}

void NCGas::seteqeps(double eqepss) {
	int totalN = obj->getN();
	eqeps = eqepss;
	double sigma = 1.0;
	WCAPotential nswca(1.0, sigma, eqeps);
	obj->setinteractions(nswca);

	// vector1<double> inta(3);
	// inta[0]=1.0;
	// inta[1]=1.0;
	// inta[2]=eqeps;
	// for(int i  = 0 ; i < totalN ; i++) {
	// 	for(int j = i+1 ; j < totalN ; j++) {
	// 		obj->change_interaction_parameters(i,j,0,inta);
	// 	}
	// }
	
}


void NCGas::setgeo(cube &bc) {
	(*(this->obj)).setgeometry(bc);
	l=bc.l;
	num = floor(bc.l/4.0);
	// vector1<double> inta(3);
	// inta[0]=1.0;
	// inta[1]=1.0;
	// inta[2]=eqeps;
	// for(int i  = 0 ; i < totalN ; i++) {
	// 	for(int j = i+1 ; j < totalN ; j++) {
	// 		obj->change_interaction_parameters(i,j,0,inta);
	// 	}
	// }
	
}

// void NCGas::seteqeps(double eqepss) {
// 	int totalN = obj->getN();
// 	eqeps = eqepss;
// 	vector1<double> inta(2);
// 	inta[0]=eqeps;
// 	inta[1]=1.0;
// 	for(int i  = 0 ; i < totalN ; i++) {
// 		for(int j = i+1 ; j < totalN ; j++) {
// 			obj->change_interaction_parameters(i,j,0,inta);
// 		}
// 	}
	
// }

matrix<double> NCGas::calculate_virial() {
	int ccc;
	int totalN = obj->getN();
	//int sibdiv = floor(ll/4.0);

	matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);


	matrix<int> *froyo = obj->calculatepairs(boxes,3.5);

	//cout << boxes << endl;

	//matrix<double> F((*obj).calculateforces(*froyo));
	//s_matrix<double> F((*obj).calculateforces(*froyo));
	//matrix<double> v(3,3);
	matrix<double> state = obj->getdat();
	matrix<double> v = obj->calculatestress(*froyo);

	//cout << v << endl;
	//cout << "calculated eq stress" << endl;
	matrix<double> v2(3,3);

		#pragma omp parallel for
		for(int p = 0 ; p < (*froyo).getNsafe() ; p++) {
			
					int i = froyo->operator()(p,0);
					int j = froyo->operator()(p,1);
			//for(int j = i+1 ; j < totalN ; j++) {
					vector1<double> uv(3);
					double dis = 0;
					
					((obj)->getgeo()).distance_vector(state, i, j, uv, dis); //returns the square distance to d and the distance between to uv
		//vector1<double> dx(dimension);
					
					// //vector1<double> as = forcecalc1((*obj).getcoordinate(i,0), (*obj).getcoordinate(j,0), (*obj).getcoordinate(i,1), (*obj).getcoordinate(j,1), (*obj).getcoordinate(i,2), (*obj).getcoordinate(j,2));
					// cout << "p1: " << (*obj).getcoordinate(i,0) << " " << (*obj).getcoordinate(i,1) << " " << (*obj).getcoordinate(i,2) << endl;
					// cout << "p2: " << (*obj).getcoordinate(j,0) << " " << (*obj).getcoordinate(j,1) << " " << (*obj).getcoordinate(j,2) << endl;
					vector1<double> as(3); 
					this->forcecalc2(i,j,as);

					
						for(int i1 = 0 ; i1 < 3 ; i1++) {
							for(int i2 = 0 ; i2 < 3 ; i2++) {
								v2(i1,i2)+=uv[i1]*as[i2];
							}
						}
					// cout << as << endl;
					// pausel();


					// F(i,0)+=as.gpcons(0);
					// F(i,1)+=as.gpcons(1);
					// F(i,2)+=as.gpcons(2);
					// F(j,0)+=-as.gpcons(0);
					// F(j,1)+=-as.gpcons(1);
					// F(j,2)+=-as.gpcons(2);
			
			//}
		}



		matrix<double> vx(3,3);
		
		for(int i = 0 ; i < 3 ; i++) {
			for(int j = 0 ; j < 3 ; j++) {
				vx(i,j)=v(i,j);
			}
		}

		for(int i = 0 ; i < 3 ; i++) {
			for(int j = 0 ; j < 3 ; j++) {
				vx(i,j)+=v2(i,j);
			}
		}

		delete froyo;

		return vx;

}

matrix<double> NCGas::calculate_virial_with_matrices(matrix<double>&x,matrix<double> &F) {
	// int ccc;
	// int totalN = obj->getN();
	// //int sibdiv = floor(ll/4.0);

	// matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);


	// matrix<int> *froyo = obj->calculatepairs(boxes,3.5);

	// //cout << boxes << endl;

	// matrix<double> F((*obj).calculateforces(*froyo));
	// //s_matrix<double> F((*obj).calculateforces(*froyo));



	// 	#pragma omp parallel for
	// 	for(int p = 0 ; p < (*froyo).getNsafe() ; p++) {
	// 				int i = froyo->operator()(p,0);
	// 				int j = froyo->operator()(p,1);
	// 		//for(int j = i+1 ; j < totalN ; j++) {
				
	// 				// //vector1<double> as = forcecalc1((*obj).getcoordinate(i,0), (*obj).getcoordinate(j,0), (*obj).getcoordinate(i,1), (*obj).getcoordinate(j,1), (*obj).getcoordinate(i,2), (*obj).getcoordinate(j,2));
	// 				// cout << "p1: " << (*obj).getcoordinate(i,0) << " " << (*obj).getcoordinate(i,1) << " " << (*obj).getcoordinate(i,2) << endl;
	// 				// cout << "p2: " << (*obj).getcoordinate(j,0) << " " << (*obj).getcoordinate(j,1) << " " << (*obj).getcoordinate(j,2) << endl;
	// 				vector1<double> as(3); 
	// 				this->forcecalc2(i,j,as);

	// 				// cout << as << endl;
	// 				// pausel();


	// 				F(i,0)+=as.gpcons(0);
	// 				F(i,1)+=as.gpcons(1);
	// 				F(i,2)+=as.gpcons(2);
	// 				F(j,0)+=-as.gpcons(0);
	// 				F(j,1)+=-as.gpcons(1);
	// 				F(j,2)+=-as.gpcons(2);
			
	// 		//}
	// 	}

	int totalN = x.getNsafe();

		matrix<double> v(3,3);
		for(int p = 0 ; p < totalN ; p++)
		for(int i = 0 ; i < 3 ; i++) {
			for(int j = 0 ; j < 3 ; j++) {
				v(i,j)+=x(p,i)*F(p,j);
			}
		}
		for(int i = 0 ; i < 3 ; i++) {
			for(int j = 0 ; j < 3 ; j++) {
				v(i,j)/=(double)(totalN);
			}
		}
		//delete froyo;
		return v;

}


void NCGas::run(int runtime)
{
	int ccc;
	int totalN = obj->getN();
	//int sibdiv = floor(ll/4.0);

	matrix<int> boxes = (obj)->getgeo().generate_boxes_relationships(num,ccc);

	//double ll = ((*obj).getgeo()).l;	
	
	matrix<int> *froyo = obj->calculatepairs(boxes,3.5);

	//matrix<int> *froyo2 = obj->calculatepairs_sorted(boxes,3.5);
	// cout << "pairs calculated" << endl;
	// ofstream myfilep;
	// myfilep.open("froyo.csv");
	// myfilep <<= *froyo;
	// myfilep.close();
	// pausel();
	// ofstream myfilep2;
	// myfilep2.open("froyo2.csv");
	// myfilep2 <<= *froyo2;
	// myfilep2.close();
	// pausel();

	matrix<double> state(obj->getdat()); //the state of the system

	int every = 1000;
	int i;
	bool up = false;

	for(i = 0 ; i < runtime ; i++) {
		cout << i << endl;
		cout << (*obj).avmom() << endl;
	if(i%15==0) {
		delete froyo;
		// cout << "updated after: " << i << endl;
		// state = obj->getdat();
		froyo = (*obj).calculatepairs(boxes,3.5); 
	}


		
		matrix<double> F((*obj).calculateforces(*froyo));


		

		int dd = obj->getdimension();
		matrix<double> R(obj->getN(),dd);
		for(int i = 0 ; i < obj->getN() ; i++) {
			for(int j = 0 ; j < dd ; j++) {
				R(i,j) = (3.464101615 * ((double) rand() / (RAND_MAX)) - 1.732050808);
			}
		}
		
		// #pragma omp parallel for
		// for(int i = 0 ; i < obj->getN() ; i++) {
		// 			// int i = froyo->operator()(p,0);
		// 			// int j = froyo->operator()(p,1);
		// 	//for(int j = i+1 ; j < totalN ; j++) {
				
		// 			//vector1<double> as = forcecalc1((*obj).getcoordinate(i,0), (*obj).getcoordinate(j,0), (*obj).getcoordinate(i,1), (*obj).getcoordinate(j,1), (*obj).getcoordinate(i,2), (*obj).getcoordinate(j,2));
					
		// 			vector1<double> as = forcecalc3(i);
		// 			F(i,0)+=as.gpcons(0);
		// 			F(i,1)+=as.gpcons(1);
		// 			F(i,2)+=as.gpcons(2);
		// 			// F(j,0)+=-as.gpcons(0);
		// 			// F(j,1)+=-as.gpcons(1);
		// 			// F(j,2)+=-as.gpcons(2);
			
		// 	//}
		// }
		matrix<double> F2(F.getNsafe(),3);

		#pragma omp parallel for
		for(int p = 0 ; p < (*froyo).getNsafe() ; p++) {
					int i = froyo->operator()(p,0);
					int j = froyo->operator()(p,1);
			//for(int j = i+1 ; j < totalN ; j++) {
				
					// //vector1<double> as = forcecalc1((*obj).getcoordinate(i,0), (*obj).getcoordinate(j,0), (*obj).getcoordinate(i,1), (*obj).getcoordinate(j,1), (*obj).getcoordinate(i,2), (*obj).getcoordinate(j,2));
					// cout << "p1: " << (*obj).getcoordinate(i,0) << " " << (*obj).getcoordinate(i,1) << " " << (*obj).getcoordinate(i,2) << endl;
					// cout << "p2: " << (*obj).getcoordinate(j,0) << " " << (*obj).getcoordinate(j,1) << " " << (*obj).getcoordinate(j,2) << endl;
					vector1<double> as(3);
					this->forcecalc2(i,j,as);

					F2(i,0)+=as.gpcons(0);
					F2(i,1)+=as.gpcons(1);
					F2(i,2)+=as.gpcons(2);
					F2(j,0)+=-as.gpcons(0);
					F2(j,1)+=-as.gpcons(1);
					F2(j,2)+=-as.gpcons(2);
			
			//}
		}	

	// matrix<double> F1((*obj).calculateforces_truncated((*froyo),1.12246));		
	// double totv = 0.0;
	// for(int i1 = 0  ; i1 < F.getNsafe() ; i1++ ) {
	// 	for(int i2 = 0  ; i2 < 3 ; i2++ ) {
	// 		totv += F1(i1,i2)*F2(i1,i2);
	// 	}
	// }
	// cout << totv << endl;
	
	// double totv2 = 0.0;
	// for(int i1 = 0  ; i1 < F.getNsafe() ; i1++ ) {
	// 	for(int i2 = 0  ; i2 < 3 ; i2++ ) {
	// 		totv2 += F(i1,i2)*F2(i1,i2);
	// 	}
	// }
	// cout << totv2 << endl;
	

	F+=F2;

	if(i>0&&i%every==0) { 
// vector1<double> x(this->getN()),y(this->getN()),z(this->getN());
// for(int i = 0 ; i < this->getN() ; i++) {
// 	x[i]=(*mom)(i,0);
// 	y[i]=(*mom)(i,1);
// 	z[i]=(*mom)(i,2);
// }

// cout << maxval(x) << endl;
// cout << maxval(y) << endl;
// cout << maxval(z) << endl;

		cout << i << endl;
		stringstream ss2;
		// ss2 <<i/every;
		// string pairlist = "list";

		 stringstream kts;
		 kts << (*obj).getkT();

		 stringstream epi;
		 epi << eps;

		 stringstream epieq;
		 epieq << eqeps;

		 stringstream len;
		 len << l;

		 string extension =  "_eqeps="+epieq.str()+"_eps="+epi.str()+"_kT="+kts.str()+"_l="+len.str()+".csv";
		// pairlist += ss2.str();
		// pairlist += extension;

		// ofstream myfilex;
		// myfilex.open(pairlist.c_str());
		// myfilex <<= froyo;
		// myfilex.close();		
	
		stringstream ss;
		ss <<i/every;
		string filename = "x";
		filename += ss.str();
		filename += extension;

		string momname = "F";
		momname += ss.str();
		momname += extension;

		string momname2 = "Fnc";
		momname2 += ss.str();
		momname2 += extension;

		string pairname = "pairs";
		pairname += ss.str();
		pairname += extension;



		ofstream myfile;
		myfile.open(filename.c_str());
		myfile <<= (*obj).getdat();
		myfile.close();

		matrix<double> F1((*obj).calculateforces_truncated((*froyo),1.12246));

		ofstream myfile2;
		myfile2.open(momname.c_str());
		myfile2 <<= F1;
		myfile2.close();

		ofstream myfile3;
		myfile3.open(momname2.c_str());
		myfile3 <<= F2;
		myfile3.close();		
		}	
		// for(int i  = 0 ; i < obj->getN() ; i++ ) {
		// 	double x1 = (*obj).getcoordinate(i,0);
		// 	double y1 = (*obj).getcoordinate(i,1);
		// 	double z1 = (*obj).getcoordinate(i,2);



		// 	double dx = x1 - ll/2.;
		// 	double dy = y1 - ll/2.;
		// 	double dz = z1 - ll/2.;
		// 	double d = SQR(dx) + SQR(dy) + SQR(dz);

		// 	F(i,0) += -(dy/(exp(d/(2.*SQR(lambda)))*SQR(lambda))) + dz/(exp(d/(2.*SQR(lambda)))*SQR(lambda));
		// 	F(i,1) +=  dx/(exp(d/(2.*SQR(lambda)))*SQR(lambda)) - dz/(exp(d/(2.*SQR(lambda)))*SQR(lambda))
		// 	F(i,2) += -(dx/(exp(d/(2.*SQR(lambda)))*SQR(lambda))) + dy/(exp(d/(2.*SQR(lambda)))*SQR(lambda))
		// }



		(*obj).advance_mom(F,R);


		(*obj).advance_pos();

		//up = (obj)->getgeo().exceed_dis(obj->getdat(), state,1.0);

	}
}

