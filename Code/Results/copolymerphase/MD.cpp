MD::MD()  {
	//dat = new matrix<double>();

}

MD::MD(const MD &old) {
	// geometry *geotemp = (old.geo);
	// //geo = &(old.geo);
	
	geometry *geotemp = (old.geo)->clone();
	geo = geotemp;
	//cout << "yo" << endl;
	//geo = old.geo;
	//dat = new matrix<double>;

	potential *potnew =  (old.ints)->clone();
	ints = potnew;

	matrix<double> *groleo = (old.dat)->clone();
	dat = groleo;

	dimension = dat->getncols();

	if (dimension != (*(this->geo)).dimension) error("dimensions must match in MD");
}


MD::~MD() {
	// delete dat;
	// delete geo;
}

void MD::setgeometry(geometry &a) {
	geometry* q = a.clone();
	geo = q;


}

void MD::setdat(matrix<double> &a) {
	matrix<double> *res =  a.clone();
	dat = res;
	dimension = dat->getncols();
	if (dimension != (*(this->geo)).dimension) error("dimensions must match in MD");
}

void MD::setinteractions(potential &a) {
	potential* q = a.clone();
	ints = q;
}



double MD::getcoordinate(int i, int j) {
	return (*dat).mat[i*dimension+j];
}

int MD::getdimension() {
	return dimension;
}

int MD::getN() {
	return (*this->dat).getNsafe();
}

matrix<double>& MD::getdat() {
	return *(this->dat);
}

geometry& MD::getgeo() {
	return *(this->geo);
}

potential& MD::getints() {
	return *ints;
}

// double MD::distance(int &i, int &j) {
// return (*this->geo).distance((*dat)[i],(*dat)[j]);
// }

// bool MD::distance_less_than(const int &i,const int &j, double R) {
// return (*this->geo).distance_less_than((*dat)[i],(*dat)[j],R);
// }




// matrix<int>* MD::calculatepairs(matrix<int> &boxlist) {
// 	vector<int> index1;
// 	vector<int> index2;
// 	vector<int> index3;
// 	//#pragma omp parallel for
// 	int Ns = (*dat).getNsafe();
// 	int np = ints.number_of_potentials();
// 	//vector1<double> int_dis(np);
// 	vector1<double> int_dis(np);
// 	vector1<bool> int_dl(np);
	
// 	for(int i = 0 ; i < np ; i++) {
// 		int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
// 	}
	
// 	for(int i = 0 ; i < np ; i++) {
// 		int_dl[i] = (ints.access_potential(i).dl);
// 	}


// 	int total_cubes = boxlist.getNsafe();

// 	double dims = (double)this->getdimension();

// 	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

// 	vector<vector<int> > b;
// 	b.reserve(total_cubes);
//     for(int j = 0 ; j < total_cubes ; j++) {
//         vector<int> temp;
//         b.push_back(temp);
//     }
//     int dimension = this->getdimension();


// 	vector1<int> dim(dimension);
// 	for(int i = 0 ; i < dimension ; i++ ) {
// 		int ij = 1;
// 		for(int j = 0 ; j < i ; j++ ) {
// 		ij*= cubes_per_length;
// 		}
// 		dim[i] = ij;
// 	}



// 	for(int i = 0 ; i < this->getN() ; i++) {

// 		int c = geo->assign_box((*dat)[i],dim,cubes_per_length);

// 		b[c].push_back(i);
// 	}



// 	int ss = boxlist[0].getsize();

// 	#pragma omp parallel
// 	{
// 	vector<int> index1_private;
// 	vector<int> index2_private;
// 	vector<int> index3_private;
// 	#pragma omp for nowait schedule(static)
// 	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
// 		for(int c2 = 0 ; c2 < ss ; c2++) {
// 			int box1 = c1;
// 			int box2 = boxlist(c1,c2);
// 			if(box1==box2) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = i+1 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							int q = ints.get_potential_number(iterator1,iterator2,l);

// 							if(int_dl[q]) {
// 								bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 								if(cond) {
// 									index1_private.push_back(iterator1);
// 									index2_private.push_back(iterator2);
// 									index3_private.push_back(l);					
// 								}

// 							}
// 							else {
// 								index1_private.push_back(iterator1);
// 								index2_private.push_back(iterator2);
// 								index3_private.push_back(l);
// 							}							
// 						}
// 					}
// 				}

// 			}
// 			else if(box2>box1) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = 0 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							int q = ints.get_potential_number(iterator1,iterator2,l);

// 							if(int_dl[q]) {
// 								bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 								if(cond) {
// 									index1_private.push_back(iterator1);
// 									index2_private.push_back(iterator2);
// 									index3_private.push_back(l);					
// 								}

// 							}
// 							else {
// 								index1_private.push_back(iterator1);
// 								index2_private.push_back(iterator2);
// 								index3_private.push_back(l);
// 							}							
// 						}						
// 					}
// 				}
// 			}
// 			else {

// 			}
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
// 	}

// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index3.insert(index3.end(),index3_private.begin(),index3_private.end());
// 	}
// 	}


// 	matrix<int> * a = new s_matrix<int>(index1.size(),3);
// 	//s_matrix<int> pairs(index1.size(),3);
// 	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
// 		(*a)(i,0) = index1[i];
// 		(*a)(i,1) = index2[i];
// 		(*a)(i,2) = index3[i];
// 	}

	
	
// 	return a;


// }


// matrix<int>* MD::calculatepairs_global(matrix<int> &boxlist, double cut_off) {
// 	vector<int> index1;
// 	vector<int> index2;
// 	vector<int> index3;
// 	//#pragma omp parallel for
// 	// int Ns = (*dat).getNsafe();
// 	// int np = ints.number_of_potentials();
// 	// //vector1<double> int_dis(np);
// 	// vector1<double> int_dis(np);
// 	// vector1<bool> int_dl(np);
	
// 	// for(int i = 0 ; i < np ; i++) {
// 	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
// 	// }
	
// 	// for(int i = 0 ; i < np ; i++) {
// 	// 	int_dl[i] = (ints.access_potential(i).dl);
// 	// }


// 	int total_cubes = boxlist.getNsafe();

// 	double dims = (double)this->getdimension();

// 	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

// 	vector<vector<int> > b;
// 	b.reserve(total_cubes);
//     for(int j = 0 ; j < total_cubes ; j++) {
//         vector<int> temp;
//         b.push_back(temp);
//     }
//     int dimension = this->getdimension();


// 	vector1<int> dim(dimension);
// 	for(int i = 0 ; i < dimension ; i++ ) {
// 		int ij = 1;
// 		for(int j = 0 ; j < i ; j++ ) {
// 		ij*= cubes_per_length;
// 		}
// 		dim[i] = ij;
// 	}



// 	for(int i = 0 ; i < this->getN() ; i++) {

// 		int c = geo->assign_box((*dat)[i],dim,cubes_per_length);

// 		b[c].push_back(i);
// 	}



// 	int ss = boxlist[0].getsize();

// 	#pragma omp parallel
// 	{
// 	vector<int> index1_private;
// 	vector<int> index2_private;
// 	vector<int> index3_private;
// 	#pragma omp for nowait schedule(static)
// 	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
// 		for(int c2 = 0 ; c2 < ss ; c2++) {
// 			int box1 = c1;
// 			int box2 = boxlist(c1,c2);
// 			if(box1==box2) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = i+1 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];
// 						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);

// 							//if(int_dl[q]) {
// 						bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
// 						if(cond) {
// 							for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);
// 							index1_private.push_back(iterator1);
// 							index2_private.push_back(iterator2);
// 							index3_private.push_back(l);	
// 							}				
// 						}							
// 					}
// 				}
// 			}
// 			else if(box2>box1) {
// 				for(int i = 0 ; i < b[box1].size() ; i++) {
// 					for(int j = 0 ;  j < b[box2].size() ; j++) {
// 						int iterator1 = (b[box1])[i];
// 						int iterator2 = (b[box2])[j];

// 						bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
// 						if(cond) {
// 							for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 							//int q = ints.get_potential_number(iterator1,iterator2,l);
// 							index1_private.push_back(iterator1);
// 							index2_private.push_back(iterator2);
// 							index3_private.push_back(l);	
// 							}				
// 						}						
// 						// for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
// 						// 	int q = ints.get_potential_number(iterator1,iterator2,l);

// 						// 	if(int_dl[q]) {
// 						// 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
// 						// 		if(cond) {
// 						// 			index1_private.push_back(iterator1);
// 						// 			index2_private.push_back(iterator2);
// 						// 			index3_private.push_back(l);					
// 						// 		}

// 						// 	}
// 						// 	else {
// 						// 		index1_private.push_back(iterator1);
// 						// 		index2_private.push_back(iterator2);
// 						// 		index3_private.push_back(l);
// 						// 	}							
// 						// }						
// 					}
// 				}
// 			}
// 			else {

// 			}
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
// 	}
// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
// 	}

// 	#pragma omp for schedule(static) ordered
// 	for(int i = 0 ; i < omp_get_num_threads(); i++) {
// 		#pragma omp ordered
// 		index3.insert(index3.end(),index3_private.begin(),index3_private.end());
// 	}
// 	}
// 	matrix<int> * a = new s_matrix<int>(index1.size(),3);
// 	//s_matrix<int> pairs(index1.size(),3);
// 	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
// 		(*a)(i,0) = index1[i];
// 		(*a)(i,1) = index2[i];
// 		(*a)(i,2) = index3[i];
// 	}

	
	
// 	return a;


// }

matrix<int>* MD::calculatepairs(matrix<int> &boxlist, double cut_off) {
	vector<int> index1;
	vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b;
	b.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
        b.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}



	for(int i = 0 ; i < this->getN() ; i++) {

		int c = geo->assign_box((*dat),i,dim,cubes_per_length);

		b[c].push_back(i);
	}



	int ss = boxlist.getncols();

	#pragma omp parallel
	{
	vector<int> index1_private;
	vector<int> index2_private;
	//vector<int> index3_private;
	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
		for(int c2 = 0 ; c2 < ss ; c2++) {
			int box1 = c1;
			int box2 = boxlist(c1,c2);
			if(box1==box2) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = i+1 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];
						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
						bool cond = geo->distance_less_than(*dat,iterator1,iterator2,cut_off);
						//bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(iterator1);
							index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
							//}				
						}							
					}
				}
			}
			else if(box2>box1) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = 0 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];

						bool cond = geo->distance_less_than(*dat,iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(iterator1);
							index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
						//	}				
						}						
						// for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
						// 	int q = ints.get_potential_number(iterator1,iterator2,l);

						// 	if(int_dl[q]) {
						// 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
						// 		if(cond) {
						// 			index1_private.push_back(iterator1);
						// 			index2_private.push_back(iterator2);
						// 			index3_private.push_back(l);					
						// 		}

						// 	}
						// 	else {
						// 		index1_private.push_back(iterator1);
						// 		index2_private.push_back(iterator2);
						// 		index3_private.push_back(l);
						// 	}							
						// }						
					}
				}
			}
			else {

			}
	}
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
	}
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index2.insert(index2.end(),index2_private.begin(),index2_private.end());
	}

	}
	matrix<int> * a = new matrix<int>(index1.size(),2);
	//s_matrix<int> pairs(index1.size(),3);
	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
		(*a)(i,0) = index1[i];
		(*a)(i,1) = index2[i];
	}

	
	
	return a;


}



matrix<int>* MD::calculatepairs_sorted(matrix<int> &boxlist, double cut_off) {
	vector<mdpair> index1;
	//vector<int> index2;
	//#pragma omp parallel for
	// int Ns = (*dat).getNsafe();
	// int np = ints.number_of_potentials();
	// //vector1<double> int_dis(np);
	// vector1<double> int_dis(np);
	// vector1<bool> int_dl(np);
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dis[i] = 1.4*(ints.access_potential(i).interaction_distance);
	// }
	
	// for(int i = 0 ; i < np ; i++) {
	// 	int_dl[i] = (ints.access_potential(i).dl);
	// }


	int total_cubes = boxlist.getNsafe();

	double dims = (double)this->getdimension();

	int cubes_per_length = (int)round(exp(log(total_cubes)/dims));

	vector<vector<int> > b;
	b.reserve(total_cubes);
    for(int j = 0 ; j < total_cubes ; j++) {
        vector<int> temp;
        b.push_back(temp);
    }
    //int dimension = this->getdimension();


	vector1<int> dim(dimension);
	for(int i = 0 ; i < dimension ; i++ ) {
		int ij = 1;
		for(int j = 0 ; j < i ; j++ ) {
		ij*= cubes_per_length;
		}
		dim[i] = ij;
	}



	for(int i = 0 ; i < this->getN() ; i++) {

		int c = geo->assign_box((*dat),i,dim,cubes_per_length);

		b[c].push_back(i);
	}



	int ss = boxlist.getncols();


	#pragma omp parallel
	{
	vector<mdpair> index1_private;
	//vector<int> index2_private;
	//vector<int> index3_private;
	#pragma omp for nowait schedule(static)
	for(int c1 = 0 ; c1 < boxlist.getNsafe() ; c1++)
		for(int c2 = 0 ; c2 < ss ; c2++) {
			int box1 = c1;
			int box2 = boxlist(c1,c2);

			if(box1==box2) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = i+1 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];
						//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);

							//if(int_dl[q]) {
						bool cond = geo->distance_less_than(*dat,iterator1,iterator2,cut_off);
						//bool cond =  distance_less_than(iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(mdpair(iterator1,iterator2));
							//index2_private.push_back(iterator2);
							//index3_private.push_back(0);	
							//}				
						}							
					}
				}
			}
			else if(box2>box1) {
				for(int i = 0 ; i < b[box1].size() ; i++) {
					for(int j = 0 ;  j < b[box2].size() ; j++) {
						int iterator1 = (b[box1])[i];
						int iterator2 = (b[box2])[j];

						bool cond = geo->distance_less_than(*dat,iterator1,iterator2,cut_off);
		
						if(cond) {
							//for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
							//int q = ints.get_potential_number(iterator1,iterator2,l);
							index1_private.push_back(mdpair(iterator1,iterator2));
							//index3_private.push_back(0);	
						//	}				
						}						
						// for(int l = 0 ; l < ints.getnints(iterator1,iterator2) ; l++) {
						// 	int q = ints.get_potential_number(iterator1,iterator2,l);

						// 	if(int_dl[q]) {
						// 		bool cond =  distance_less_than(iterator1,iterator2,int_dis[q]);
				
						// 		if(cond) {
						// 			index1_private.push_back(iterator1);
						// 			index2_private.push_back(iterator2);
						// 			index3_private.push_back(l);					
						// 		}

						// 	}
						// 	else {
						// 		index1_private.push_back(iterator1);
						// 		index2_private.push_back(iterator2);
						// 		index3_private.push_back(l);
						// 	}							
						// }						
					}
				}
			}
			else {

			}
	}
	#pragma omp for schedule(static) ordered
	for(int i = 0 ; i < omp_get_num_threads(); i++) {
		#pragma omp ordered
		index1.insert(index1.end(),index1_private.begin(),index1_private.end());
	}
	// #pragma omp for schedule(static) ordered
	// for(int i = 0 ; i < omp_get_num_threads(); i++) {
	// 	#pragma omp ordered
	// 	index2.insert(index2.end(),index2_private.begin(),index2_private.end());
	// }
	}
	sort(index1.begin(),index1.end());
	matrix<int> * a = new matrix<int>(index1.size(),2);
	//s_matrix<int> pairs(index1.size(),3);
	for(int i = 0 ; i < (*a).getNsafe() ; i++) {
		(*a)(i,0) = index1[i].a;
		(*a)(i,1) = index1[i].b;
	}

	
	
	return a;


}

matrix<double> MD::calculateforces_old(matrix<int> &pairs) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();

	#pragma omp parallel for
	for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {
		int p1 = pairs.mat[i*2+0];
		int p2 = pairs.mat[i*2+1];
		//int i1 = pairs(i,2);

		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		geo->distance_vector(*dat,p1,p2,un,dis);
		double f1  = (*ints).force(sqrt(dis));

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j]/sqrt(dis);
			(forces).mat[p1*dimension+j] += fac;
			(forces).mat[p2*dimension+j] += -fac;
		}
	}



	return forces;
}

matrix<double> MD::calculateforces(matrix<int> &pairs) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();

	#pragma omp parallel for
	for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo->distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo->distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j];
			(forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			(forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
	}



	return forces;
}

matrix<double> MD::calculatestress(matrix<int> &pairs) {
	
	//matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();
	matrix<double> stress(dimension,dimension);
	#pragma omp parallel for
	for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo->distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo->distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);

		for(int k = 0 ; k < dimension ; k++) {
		for(int j = 0 ; j < dimension ; j++) {
			
			double fac = f1*sqrt(dis)*un[k]*un[j];
			stress(k,j) += fac;
			// double fac2 = 
			// (forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			// (forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
		}
	}



	return stress;
}

matrix<double> MD::calculateforces_truncated(matrix<int> &pairs, double above) {
	
	matrix<double> forces((*dat).getNsafe(),dimension);
	//vec_vec<double> outputs(pairs.getn());
	// ofstream myfile;
	// myfile.open("forces.csv", ios::out | ios::app);
	//cout << pairs.getNsafe() << endl;
	//potential * pot = ints(0,1,0).clone();

	#pragma omp parallel for
	for(int i = 0 ; i < pairs.getNsafe() ; i++ ) {


		double dis;
		//vector1<double> un = unitvector((*dat)[p1],(*dat)[p2],dis);
		vector1<double> un(dimension);
		//int p1 = pairs.mat[i*2+0];
		//int p2 = pairs.mat[i*2+1];
		//geo->distance_vector(*dat,p1,p2,un,dis);

		int p1= pairs.mat[i*2+0]*dimension;
		int p2 =pairs.mat[i*2+1]*dimension;

		geo->distance3v((*dat).mat[p1+0],(*dat).mat[p2+0],(*dat).mat[p1+1],(*dat).mat[p2+1],(*dat).mat[p1+2],(*dat).mat[p2+2],un,dis);

		double f1  = (*ints).force2mdx(dis);
		if(dis < above) f1 = 0;

		for(int j = 0 ; j < dimension ; j++) {
			double fac = f1*un[j];
			(forces).mat[pairs.mat[i*2+0]*dimension+j] += fac;
			(forces).mat[pairs.mat[i*2+1]*dimension+j] += -fac;
		}
	}



	return forces;
}





//void MD::adv(matrix<int> &pairs) { error("adv function not overriden matrix");}