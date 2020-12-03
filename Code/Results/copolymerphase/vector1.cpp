#include "basic.h"
#include "vector1.h"

#ifndef VECTOR1_CPP
#define VECTOR1_CPP


template <class T>
inline vector1<T>::vector1() { //null constructor
	size = 1;
	data = new T [size];
	data[0] = T(0);
}

template <class T>
vector1<T>::vector1(int n) { //create a vector of n elements
	size = n;
	data = new T [size];
	for ( int i = 0 ; i < n ; i++ )
		data[i] = T(0);
}
template <class T>
vector1<T>::vector1(int n, T a) { //create a vector of n elements initialized to a.
	size = n;
	data = new T [size];
	for ( int i = 0 ; i < n ; i++ )
		data[i] = a;
}

template <class T>
vector1<T>::vector1(const vector1<T> &v) { //copy constructor
	size = v.size;
	data = new T [size] ;
	//if (!vec)
	//	error("allocation failure in vector1::vector1(vector1&)!");
	for (int i = 0 ; i < size ; ++i) data[i] = v.data[i];
	}

template <class T>
vector1<T>::vector1(const vector1<T> &v, int n) {

    int p = v.size;
    size = p+n;
    data =  new T [size] ;
    for(int i = 0 ; i < p ; i++) {data[i] = v.data[i]; }
    for(int i = p ; i<size ;i++) {data[i] = T(0); } 
}

template <class T>
vector1<T>::vector1(int n, const T *d) { //vector from pointer to array d;
	size = n;
	data = new T [size];
	for ( int i = 0 ; i < size ; i++ ) data[i]=d[i];
}
template <class T>
vector1<T>::vector1(int n, T d[]) { // vector of elements from static array T d[], n is the size of static array
	size = n;
	data = new T [size];
	for ( int i = 0 ; i < size ; i++ ) data[i] = d[i];
}



template <class T> //destructor
inline vector1<T>::~vector1() { delete[] data; }


template <class T>
void vector1<T>::delete_element(int el) {
 if (el<0 || el > size-1) {
     cout << "index out of range in vector1 operator[]" << endl;
     cout << "current choice = " << el << endl;
     cout << "maxvalue = " << size << endl;
     error("i not in range");
 }
if (size < 1 ) { error("cannot delete elements in a vector of size 1"); }

  T* data2 =  new T [size-1];
  int iter = 0;
  for(int i = 0 ; i < size ; i++ ) {
  if(i == el ) { continue;}
  else { data2[iter] = data[i]; iter++; }
  }
  delete data;
  size = size - 1;
  data = new T[size];
  for(int i = 0 ; i < size ; i++ ) {
	data[i] = data2[i];
  }
  

}

// template <class T>
// inline T& vector1<T>::operator[](const int &i) { //safe method for accessing element
//  // if (i<0 || i > size-1) {
//  //     cout << "index out of range in vector1 operator[]" << endl;
//  //     cout << "current choice = " << i << endl;
//  //     cout << "maxvalue = " << size << endl;
//  //     error("i not in range");
//  // }
//  return data[i];
// }

// template <class T>
// inline T& vector1<T>::operator()(int i) { //unsafe method for accessing element
// 	//if ( i < 0 || i >= size ) error("index out of range in vector");
// 	return data[i];
// }




template <class T>
inline T& vector1<T>::gpcons(int i) const { //const
  	//if ( i < 0 || i >= size ) error("index out of range in vector");
	return data[i];  
}

template <class T>
inline T& vector1<T>::periodicselect(int i) {
    if( i < 0 ) {
    while ( i < 0 ) {
        i += size;
        if(i > 0) break;
    }
    }
    if ( i > size-1 ) {
    while ( i > size-1) {
        i -= size;
        if( i < size -1) break;
    }
    }
    return data[i];
    
}

template <class T> //unray plus, return the vector itself.
inline vector1<T>& vector1<T>::operator+() { return *this; }

template <class T>
vector1<T>& vector1<T>::operator=(const vector1<T> &v) { //assigment operator

    delete[] data;
    size = v.size;
    data = new T [size];
    for(int i = 0 ; i < size ; i++) {
        data[i] =  v.data[i];
    }
    return *this;
}

template <class T>
vector1<T>& vector1<T>::operator+=(const vector1<T>&v) //add a vector to another
{ if (size != v.size)
        error("diff size in vector1<T>& vector1<T>::op=(const vector1<T>&)!");
for (int i = 0 ; i < size; ++i) data[i] += v.data[i];
return *this;
}

template <class T>
vector1<T>& vector1<T>::operator-=(const vector1<T>&v) //minus a vector from another
{ if (size != v.size)
    error("diff size in vector1<T>& vector1<T>::op=(const vector1<T>&)!");
        for (int i = 0 ; i < size; ++i) data[i] -= v.data[i];
return *this;
}

template <class T>
vector1<T>& vector1<T>::operator+=(T x) { //add element x onto all elements of vector
for ( int i = 0 ; i < size ; ++i) data[i] += x;
return *this;

}

template <class T>
vector1<T>& vector1<T>::operator-=(T x) { //subtract element x from all elements of vector
for ( int i = 0 ; i < size ; ++i) data[i] -= x;
return *this;
}

template <class T>
vector1<T>& vector1<T>::operator*=(T x) //multiply vector by scalar
{ for (int i = 0 ; i < size; ++i )
data[i] *= x;
return *this;
}

template <class T>
vector1<T>& vector1<T>::operator/=(T x) //divide vector by scalar
{
//if ( x == 0 ) error("cannot divide by zero");
for (int i = 0 ; i < size; ++i )
data[i] /= x;
return *this;
}

template<class T>
bool& vector1<T>::operator==(const vector1<T> &x) {
    if( size != x.size) return false;
    else {
        for(int i = 0 ; i < size ; i++) {
            if( abs(data[i]-x.data[i])>1E-10) return false;
        }
        return true;
    }
}


template <class T>
void vector1<T>::swap(int i, int j){ // swap two elements
T tmp = data[i]; //temporary storage
data[i] = data[j];
data[j] = tmp;
}

template <class T>
void vector1<T>::resize(int n) { /*resize vector whilst destroying all data*/
    size = n; //set new size
    delete [] data; // delete data
    data = new T [n]; // initialize new data
    for(int i=0;i<n;i++)
        data[i]=0.0; //data is destroyed
}

template <class T>
vector1<T> operator-(const vector1<T> &v){ //unaray minus operator
int n = v.size;
vector1<T> u(n);
for ( int i = 0 ; i < n; ++i)
         u.data[i] = -v.data[i];
return u;
}

template <class T>
vector1<T> operator+(const vector1<T> &v1, const vector1<T> &v2) { //one vector add another
int n = v1.size;
 if ( v1.size != v2.size) error("vector1s are not of same dimension");
 vector1<T> v(n);
 for ( int i = 0 ; i < n; ++i)
         v.data[i] = v1.data[i]+v2.data[i]; //add elements
 return v;
}

template <class T>
bool operator==(const vector1<T> &v1, const vector1<T> &v2) { //check whether vectors are the same
if (v1.size != v2.size) return false;
for ( int i = 0 ; i < v1.size ; i++ ) {
if ( v1.data[i] != v2.data[i] ) return false; /* beware of this condition if the vector stores doubles*/
}
return true;

}


template <class T>
vector1<T> operator-(const vector1<T> &v1, const vector1<T> &v2) { /* one vector subtract another*/
int n = v1.size;
 if ( v1.size != v2.size) error("vector1<T>s are not of same dimension");
 vector1<T> v(n);
 for ( int i = 0 ; i < n; ++i)
         v.data[i] = v1.data[i]-v2.data[i];
 return v;
}

template <class T>
vector1<T> operator*(const vector1<T> &v, T d) { /* one vector multiplied by scalar d*/
vector1<T> vd = v;
vd *= d ;
return vd;
}

template <class T, class Q>
vector1<T> operator*(const vector1<T> &v, Q d) { /* multiply vector by scalar of different type*/
    vector1<T> vd = v;
    T dd(d); //convert Q to a T (if such a conversion exists)
    vd *= dd;
    return vd;
}

template <class T>
vector1<T> operator*(T d, const vector1<T> &v) { //opposite
return v * d ;
}

template <class T>
vector1<T> operator/(const vector1<T> &v, T d)	{ //divide vector by scalar
vector1<T> vd = v;
vd /= d;
return vd;
}

template <class Y>
vector1<Y> operator/(const vector1<Y> &a, const vector1<Y> &b) {
    if(a.getsize() != b.getsize() ) error("vector sizes not the same in operator/");
    
    vector1<Y> vd(a.getsize());
    for(int i = 0 ; i < a.getsize() ; i++)
        vd.data[i]=a.data[i]/b.data[i];
    
    return vd;
    
    
}

template <class T>
vector1<T> operator&(const vector1<T> &v, const vector1<T> &u) {
        /*  return vectors v and u with all elements multiplied by each other */
        int n = v.size;
        if ( n != u.size ) error("vectors not of same dimension in %");
        vector1<T> sol(n);
        for ( int i = 0 ; i < n ; i++ ) {
                sol.data[i] = v.data[i] * u.data[i] ; //return each element multiplied by the other corresponding index
        }
        return sol;
}

template <class T>
T minval(const vector1<T> &a) {
    T min = a.data[0];
    for(int i=1;i<a.size;i++) {
        if(  a.data[i] < min ) {
            min = a.data[i];
        }
    }
    return min;
}

template <class T>
T maxvalabs(const vector1<T> &a) {
    T max = abs(a.data[0]);
    for(int i=1;i<a.size;i++) {
        if(  abs(a.data[i]) > max ) {
            max = abs(a.data[i]);
        }
    }
    return max;
}

template <class T>
T minvalabs(const vector1<T> &a) {
    T min = abs(a.data[0]);
    for(int i=1;i<a.size;i++) {
        if(  abs(a.data[i]) < min ) {
            min = abs(a.data[i]);
        }
    }
    return min;
}

template <class T>
T maxval(const vector1<T> &a) {
    T max = a.data[0];
    for(int i=1;i<a.size;i++) {
        if(  a.data[i] > max ) {
            max = a.data[i];
        }
    }
    return max;
}

template <class T>
T elementmultiply(const vector1<T> &a) {
T res= 1;
for(int i = 0 ; i < a.size ; i++ ) {
res *= a.data[i];
}
return res;
}

template <class T>
int minindex(const vector1<T> &a) {
    T min = a.data[0];
    int j = 0;
    for(int i=1;i<a.size;i++) {
        if(  a.data[i] < min ) {
            min = a.data[i];
            j = i;
        }
    }
    return j;
}

template <class T>
int maxindex(const vector1<T> &a) {
    T max = a.data[0];
    int j  = 0;
    for(int i=1;i<a.size;i++) {
        if(  a.data[i] > max ) {
            max = a.data[i];
            j = i;
        }
    }
    return j;
}

template <class T>
double scalar( const vector1<T> &u, const vector1<T> &v ) { //scalar product
        double t = 0  ;
        int n = u.size ;
        if (u.size != v.size ) error("vector1s are not of same dimension");
        for (int i = 0 ; i < n ; ++i )
                t += double(u.data[i] * v.data[i]); //dot(scalar) product
        return t;
}

template <class T>
vector1<T> unitvector(const vector1<T> &u, const vector1<T> &v) {
    if (u.size != v.size || u.size == 0 ) error("vector1s are not of correct dimension in vector::unitvector");
    int n = u.size;
    vector1<T> res(n);
    double normfactor = 0.0;
    for(int i = 0 ; i < n ; i++) {
        res[i] = u.gpcons(i)-v.gpcons(i);
        normfactor += SQR(res[i]);
    }
    normfactor = sqrt(normfactor);
    if(normfactor < 1E-10) return vector1<T>(n);
    else return (1./normfactor)*res;

}

template <class T>
vector1<T> unitvector(const vector1<T> &u, const vector1<T> &v, double &norm) {
    if (u.size != v.size || u.size == 0 ) error("vector1s are not of correct dimension in vector::unitvector");
    int n = u.size;
    vector1<T> res(n);
    double normfactor = 0.0;
    for(int i = 0 ; i < n ; i++) {
        res[i] = u.gpcons(i)-v.gpcons(i);
        normfactor += SQR(res[i]);
    }
    normfactor = sqrt(normfactor);
    norm =  normfactor;
    if(normfactor < 1E-10) return vector1<T>(n);
    else return (1./normfactor)*res;

}

template <class T>
double norm( const vector1<T> &u, const vector1<T> &v) { //distance between two cartesian vectors
        double t = 0;
        int n = u.size;
        if ( n != v.size ) error("vector1s are not of same dimension");
        for ( int i = 0 ; i < n ; i++ ) {
                t += double((u.data[i]-v.data[i])*(u.data[i]-v.data[i])); //cartesian distance between two vectors
        }
        return sqrt(t);
}

template <class T>
double norm2( const vector1<T> &u) {
        double t = 0;
        int n = u.size;
        for ( int i = 0 ; i < n ; i++ ) {
                t += double(u.data[i]*u.data[i]); //dot product of vector with itself
        }
        return sqrt(t);
}

template <class T>
void vector1<T>::trimult(vector1<T> &d, vector1<T> &a, vector1<T> &c, vector1<T> &x) {
	/* matrix multiplication of x by a triadiagonal matrix with diagonal a, +1 diagonal = b and
	-1 diagonal = c */
	int n = x.size;
	data[0] = a.data[0]*x.data[0] + c.data[0] * x.data[1];
	for ( int i = 1 ; i < n-1 ;i++ ) {
	data[i] = a.data[i]*x.data[i] + c.data[i] * x.data[i+1] + d.data[i]*x.data[i-1];
	}
	data[n-1] = a.data[n-1]*x.data[n-1] + d.data[n-1]*x[n-2];
}


//vector solution to tridiagonal matrix equation, using thomas algorithm
template <class T>
void vector1<T>::trisols(vector1<T> &a, vector1<T> &b, vector1<T> c, vector1<T> &d) {
	// Using the thomas algorithm
        // a diagonal - 1
	// b diagonal
	// c diagonal + 1
	// d are the solutions
	int n = a.size;
	/* Modify the coefficients. */
	c.data[0] /= b.data[0];	/* Division by zero risk. */
	d.data[0] /= b.data[0];	/* Division by zero would imply a singular matrix. */
	for (int i = 1; i < n; i++){
		long double id = 1.0 / (b.data[i] - c.data[i-1] * a.data[i]);  /* Division by zero risk. */
		c.data[i] *= id;	                         /* Last value calculated is redundant. */
		d.data[i] = (d.data[i] - d.data[i-1] * a.data[i]) * id;
	}

	/* Now back substitute. */
	data[n - 1] = d.data[n - 1];
	for (int i = n - 2; i >= 0; i--)
		data[i] = d.data[i] - c.data[i] * data[i + 1];
}



template <class T>
void vector1<T>::trisols(const vector1<T> &a,const vector1<T> &b,const vector1<T> c,const vector1<T> &d){
	// a diagonal - 1
	// b diagonal
	// c diagonal + 1
	// d are the solutions
	int n = a.size;
	/* Modify the coefficients. */
	c.data[0] /= b.data[0];	/* Division by zero risk. */
	d.data[0] /= b.data[0];	/* Division by zero would imply a singular matrix. */
	for (int i = 1; i < n; i++){
		long double id = 1.0 / (b.data[i] - c.data[i-1] * a.data[i]);  /* Division by zero risk. */
		c.data[i] *= id;	                         /* Last value calculated is redundant. */
		d.data[i] = (d.data[i] - d.data[i-1] * a.data[i]) * id;
	}

	/* Now back substitute. */
	data[n - 1] = d.data[n - 1];
	for (int i = n - 2; i >= 0; i--)
		data[i] = d.data[i] - c.data[i] * data[i + 1];
}

//cyclic solution (William H Press et al, Numerical Recipes in C++, 3rd Edition, Chapter 2)
// this is a tridiagonal matrix equation solver where the upper right and lower left elements are
// non- zero
template <class T>
void vector1<T>::cyclic(const vector1<T> &a, const vector1<T> b, const vector1<T> &c, const T &alpha, const T &beta, const vector1<T> &r) {
    
    int i,n=a.size;
    vector1<T> u(n),z(n);
    T fact,gamma;
    gamma = -b.data[0];
    b.data[0]=b.data[0]-gamma;
    b.data[n-1]=b.data[n-1]-alpha*beta/gamma;
    
    this->trisols(a,b,c,r);
    u.data[0]=gamma;
    u.data[n-1]=alpha;
    z.trisols(a,b,c,u);
    fact=(data[0]+beta*data[n-1]/gamma)/(1.0+z.data[0]+beta*z.data[n-1]/gamma);
    for(i=0;i<n;i++) data[i] -= fact*z[i];

}

//cyclic solution (William H Press et al, Numerical Recipes in C++, 3rd Edition, Chapter 2)
// this is a tridiagonal matrix equation solver where the upper right and lower left elements are
// non- zero
template <class T>
void cyclic(const vector1<T> &a, const vector1<T> &b, const vector1<T> &c, const T &alpha, const T &beta, const vector1<T> &r,vector1<T> &x) {
    int i,n=a.size;
    vector1<T> u(n),z(n),bb(n);
    T fact,gamma;
    gamma = -b.data[0];
    bb.data[0]=b.data[0]-gamma;
    bb.data[n-1]=b.data[n-1]-alpha*beta/gamma;
    for(int i = 1 ; i < n-1 ; i++) bb.data[i]=b.data[i];
    
    x.trisols(a,bb,c,r);
    u.data[0]=gamma;
    u.data[n-1]=alpha;
    z.trisols(a,bb,c,u);
    fact=(x.data[0]+beta*x.data[n-1]/gamma)/(1.0+z.data[0]+beta*z.data[n-1]/gamma);
    for(i=0;i<n;i++) x.data[i] -= fact*z.data[i];

}

//multiply a vector x by a cylic matrix, where alpha,beta are the elements in the corner
template <class T>
void vector1<T>::cyclicmult(vector1<T> &a, vector1<T> &b, vector1<T> &c, T alpha, T beta, vector1<T> &x) {
    int n=b.size;
    data[0]=b.data[0]*x.data[0] + c.data[0] * x.data[1] + alpha*x.data[n-1];
	for ( int i = 1 ; i < n-1 ;i++ ) {
	data[i] = b.data[i]*x.data[i] + c.data[i] * x.data[i+1] + a.data[i]*x.data[i-1];
	}
    data[n-1]=b.data[n-1]*x.data[n-1]+a.data[n-1]*x.data[n-2]+beta*x.data[0];
}


//output operator, formatted to be understood by mathematica
template <class T>
ostream& operator<<=(ostream &s, const vector1<T> &v) {
 s.precision(10);
int n = v.size;

for (int i = 0 ; i< n ; ++i)
if ( i == n-1) s  <<v.data[i];
else if ( i == 0 ) s  <<v.data[i] << ",";
else s  <<v.data[i] << ",";

return s;   
}

template <class T>
ostream& operator<<(ostream &s, const vector1<T> &v) {
s.precision(10);
int n = v.size;

for (int i = 0 ; i< n ; ++i)
if ( i == n-1) s  <<v.data[i] << "}";
else if ( i == 0 ) s  << "{"   <<v.data[i] << ", ";
else s  <<v.data[i] << ", ";

return s;
}

//format two vectors for mathematica such that the final array looks like:
/*
 { { v[0],u[0]},{v[1],u[1]},...,{v[n-1],u[n-1]}}
 */
// Mainly to print the values of a certain function with it's corresponding position
template <class T>
void print_two_vectors(ofstream &s, const vector1<T>&v, const vector1<T> &u) {
    if ( v.size != u.size) error("vectors diff sizes in joint output");
    int n = v.size;
    for (int i = 0 ; i < n ; ++i)
    if ( i == n-1) s  << "{" << u.data[i] << "," << v.data[i] << "}" << "};\n";
    else if ( i == 0 ) s  << "{"   << "{" << u.data[i] << "," << v.data[i] << "}" << ", ";
    else s  << "{" << u.data[i] << "," << v.data[i] << "}" << ", ";
}


// input elements of a vector
template <class T>
istream& operator>>(istream &s, vector1<T> &v) {
int n = v.size;
cout << "enter " << n << " elements: \n";
for (int i = 0 ; i < n ; ++i)
        { cout << "v[" << i << "] = ";
        s >> v.data[i];
        }
return s;
}

template <class T> //trapezium rule for a discrete data set &v1
long double trap(const vector1<T> &v1, double dt) {
	int n = v1.size;
	double tot = 0;
	tot +=(v1.data[0]+v1.data[n-1])*0.5;
	for ( int i = 1 ; i < n-1 ; i++ ) {
		tot += v1.data[i];
	}
	tot *= dt;
	return tot;
}

template <class T>
T meanish(const vector1<T> &v1) {
	int n = v1.size;
	double tot = 0.0;
	for(int i = 0 ; i < n ; i++) {
	tot += v1.data[i];
	}
	return tot/(double)n;

}

template <class T> // trapezium rule
long double rtrap(const vector1<T> &v1, double dt) { // radial derivative starting at r = 0
	int n = v1.size;
	double tot = 0;
	tot +=(v1.data[n-1]*(n-1)*dt)*0.5;
	for ( int i = 1 ; i < n-1 ; i++ ) {
		tot += v1.data[i]*i*dt; //multiply by r
	}
	tot *= dt;
	return tot;
}

template <class T>
vector1<T> index(int &n) {
    vector1<T> a(n);
    for (int i = 0 ; i < n ; i++ ) {

        a.data[i]=T(i);

    }
    return a;
}

template <class T>
vector1<T> SQRvector(vector1<T> a) { //return a vector with elements squared
    int n = a.size;
    vector1<T> b(n);
    for(int i = 0 ; i < n ; i++) {
        b.data[i]=SQR(a.data[i]);
    }
    return b;
}

//compare two vectors, returhn the number of elements which are different
template <class T>
int compare( vector1<T> &a, vector1<T> &b) {
    double TOLMIN = 1.0E-12; //set tolerance to account for finite precision error
    int comp=0;
    int n=a.size;
    int m=b.size;
    if ( n != m ) error("diff sizes in vector comparions");
    for(int i=0; i<n ; i++) {
        if ( fabs(a.data[i]-b.data[i])<TOLMIN ) comp += 0;
        else comp+=1;
    }
    return comp;
}

// convert stl::vector to a vector1
template <class T>
vector1<T> vectovec(std::vector<T> v) {
    int n = v.capacity();
    vector1<T> u(n);
    for(int i=0;i<n;i++)
        u[i]=v[i];
    return u;

}

template <class T>
int compare(vector1<T> &a, vector1<T> &b, T eps) {
    int comp = 0;
    if( a.size != b.size) error("diff sizes in compare between vectors");
    for(int i = 0 ; i < a.size ; i++)
        if( abs(a[i]-b[i]) > eps ) comp +=1;
    
    return comp;
}

template <class Y, class T>
vector1<T> writefunction(Y &func,  T a, T b, int n) {
    T h = (b-a)/T(n-1);
    vector1<T> res(n+1);
    int ii = 0;
    for(T s = a;  s <= b; s += h) {
        res[ii] = func(s);
        ii++;
    }
    return res;
}

template <class Y>
vector1<Y> complexmultiply(vector1<Y> &a, vector1<Y> &b) {
    // multiply together two complex vectors where the real part is stored in the even indices and the corresponding
    // complex part is stored in the next odd index
    if( a.size != b.size) error("size of vectors does not agree in complex multiply");
    if( a.size % 2 != 0 ) error("vector size odd in complex multiply");
    vector1<Y> res(a.size);
    
    for(int i = 0 ; i < a.size ;) {
        
        res.data[i]= a.data[i]*b.data[i]-a.data[i+1]*b.data[i+1];
        res.data[i+1]=a.data[i]*b.data[i+1]+b.data[i]*a.data[i+1];
        i+=2;
    }
    return res;
    
}

template <class Y>
vector1<Y> smooth(const vector1<Y> &a, int s, bool open) {
    int Nz = a.size;    
    vector1<long double> temp(Nz);
        for(int j = 0 ; j < Nz ; j++) {
            long double hh = 0;
            for(int k = j-s ; k <= j+s ; k++) {
                if( k < 0 ) {
                    if(open) {
                    hh += a.data[0];   // constant outside boundaries.
                    }
                    else {
                    hh += 0;     // 0 outside boundaries
                    }
                }
                else if ( k > Nz-1) {
                    if(open) {
                    hh += a.data[Nz-1];   // constant outside boundaries.
                    }
                    else {
                    hh += 0;     // 0 outside boundaries
                    }                    
                }
                else {
                    hh += a.data[k];
                }
                
            }
            temp[j] =  hh/(2*s);
        }
    return temp;
}

/*vector1<double> randomvec(int n, double r) { // random spherical vector in n dimensions with length r
	MTRand sa;	
	double theta,phi;
	int p;
	if ( n == 1 ) {
			p = sa.randInt(1);
			if ( p == 0 ) p = -1;
			vector1<double> u(1);
			u[0] = p*r;
			return u;
	}
	else if ( n == 2 ) {
			theta = sa.rand(2*pii);
			vector1<double> u(2);
			u[0] = r*cos(theta);
			u[1] = r*sin(theta);
			return u;
	}
	else if ( n == 3 ) {
			vector1<double> u2(3);
			theta = sa.rand(2*pii);
			phi = sa.rand(pii);
			u2[0] = r*sin(theta)*cos(phi);
			u2[1] = r*sin(theta)*sin(phi);
			u2[2] = r*cos(theta);
			return u2;
	}
	else { error("Cannot work in more than 3 dimensions, for now");
	return vector1<double>(1);
	}
	

}*/

#endif





