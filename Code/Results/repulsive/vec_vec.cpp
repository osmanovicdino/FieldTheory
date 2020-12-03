template <class T>
vec_vec<T>::vec_vec() {
	n = 0;
	dat = new vector<vector1<T> >;
	// (*dat).reserve(n);
	// (*dat).push_back(vector1<T>(1));
}

template <class T>
vec_vec<T>::vec_vec(int  nn) {
	n = nn;
	dat = new vector<vector1<T> >;
	(*dat).reserve(n);
	for(int i = 0 ; i < n ; i++) {
		vector1<T> bg;
		(*dat).push_back(bg);
	}
}

template <class T>
vec_vec<T>::vec_vec(int nn, int mm) {
	n = nn;
	dat = new vector<vector1<T> >;
	(*dat).reserve(n);
	for(int i = 0 ; i < n ; i++) {
		//vector1<T> lolz(mm);
		(*this->dat).push_back(vector1<T>(mm));
	}
}

template <class T>
vec_vec<T>::vec_vec(const vec_vec<T> &a) {
	
	n = a.n;
	dat = new vector<vector1<T> >;
	(*dat).reserve(n);
	for(int i = 0 ; i < n ; i++) {
		vector1<T> v1 = (*a.dat)[i];
		(*this->dat).push_back(v1);
	}
}

template <class T>
vec_vec<T>::vec_vec(const matrix<T> &a) {
	n = a.getnrows();
	dat = new vector<vector1<T> >;
	(*dat).reserve(n);
	
	for(int i = 0 ; i < n ; i++) {
		vector1<T> v1(a.getncols());
		for(int j = 0 ; j < v1.getsize() ; j++) {
			v1[j] = a.gpcons(i,j);
		}
		(*this->dat).push_back(v1);
	}
	
}

template <class T>
vec_vec<T>::~vec_vec() {
	(*dat).clear();
	delete dat;
}

template <class T>
inline vector1<T>& vec_vec<T>::operator[](int i) {
	//if(i < 0 || i > n-1) error("choice is out of bounds in vec_vec");
	return (*this->dat)[i];
}

template <class T>
inline T& vec_vec<T>::operator()(int i, int j) {
	//if(i < 0 || i > n-1) error("choice is out of bounds in vec_vec");
	//return ((*this->dat)[i])[j];
	return this->dat->operator[](i).data[j];
}

template<class T>
vec_vec<T>& vec_vec<T>::operator=(const vec_vec<T> &a) {
	
	(*dat).clear();
	delete dat;
	n = a.n;
	dat = new vector<vector1<T> >;
	(*dat).reserve(n);
	for(int i = 0 ; i < n ; i++) {
		vector1<T> v1 = (*a.dat)[i];
		(*this->dat).push_back(v1);
	}
	return *this;	
}

template <class T>
void vec_vec<T>::remove(const int &i) { //remove the i element of the vec_vec
	(*this->dat).erase((*this->dat).begin()+i);
	n = n-1;
}

template <class T>
void vec_vec<T>::setvalue(const int &i, const vector1<T> &v1) {
//	if( (*dat)[i].getsize() != v1.getsize() ) error("different vector sizes in assignment of vector");
//	for(int j = 0 ; j < v1.getsize() ; j++) {
		(*this->dat)[i] = v1;
	
}

template <class T>
void vec_vec<T>::append(const vector1<T> &v1) {
(*this->dat).push_back(v1);
n = n+1 ;
}

template <class T>
T vec_vec<T>::max() {
	vector1<double> x(n);
	for(int i = 0 ; i < n ; i++) {
		T y = maxvalabs((*dat)[i]);
		x[i] = y;
	}
	return maxvalabs(x);
}

template<class T>
ostream& operator<<(ostream &s, const vec_vec<T> &a) {
int nr = a.n;
s<<"{";

for(int i =0 ; i < nr-1 ; i++)
	s << (*a.dat)[i] << "\n";

if(nr>=1)
s << (*a.dat)[nr-1] << "};";
return s;
}

template <class T>
ostream& operator<<=(ostream &s, const vec_vec<T> &a) {
int nr = a.n;


for(int i =0 ; i < nr-1 ; i++) {
	s <<= (*a.dat)[i];
	s << "\n";
}

if(nr>=1)
s <<= (*a.dat)[nr-1];

return s;
}

template <class T>
void outfunc(const vec_vec<T> &a, string s) {
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if(!myfileg.is_open()) error("failed to open file");
    myfileg <<= a;
    myfileg.close();
}