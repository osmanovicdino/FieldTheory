#include "basic.h"
#include "vector1.h"
#include "vector1.cpp"
#include "matrix2.h"

#ifndef MATRIX_CPP
#define MATRIX_CPP

template <class T>
matrix<T>::matrix() {
ncols = 1;
nrows = 1;
datapoints = ncols*nrows;
mat = new T [1];
mat[0] = T(0);
}

template <class T>
matrix<T>::matrix(int n) {
ncols = n;
nrows = n;
datapoints = ncols*nrows;
mat = new T [datapoints];
for ( int i = 0 ; i < datapoints ; i++ ) mat[i] = T(0);
}

template <class T>
matrix<T>::matrix(int n, int m) {
nrows = n;
ncols = m;
datapoints = n*m;
mat = new T [datapoints];
for ( int i = 0 ; i < datapoints ; i++ ) mat[i] = T(0);
}

template <class T>
matrix<T>::matrix(int n, int m, T a) {
nrows = n;
ncols = m;
datapoints = n*m;
mat = new T [datapoints];
for ( int i = 0 ; i < datapoints ; i++ ) mat[i] = T(a);
}

template <class T>
matrix<T>::matrix(int n, int m, int d, T a[] ) {
    if ( n*m != d ) error("cannot create matrix with these dimensons!");
    nrows = n;
    ncols = m;
    datapoints = n*m;
    mat = new T[datapoints];
    for(int i=0 ; i < datapoints ;i++) mat[i] = a[i];
}

template <class T>
matrix<T>::matrix(int n, int m, vector1<T> &a) {
    nrows = n;
    ncols = m;
    datapoints = n*m;
    if ( a.getsize() != datapoints ) error("incorrect size for vector in matrix constructor");
    mat = new T[datapoints];
    
    for(int i=0;i<datapoints;i++) mat[i]=a[i];
}


template <class T>
matrix<T>::matrix(const matrix<T> &m) {
	nrows = m.nrows;
	ncols = m.ncols;
	datapoints = nrows*ncols;
	mat = new T [datapoints];
	for ( int i = 0 ; i < datapoints ; i++ ) mat[i] = m.mat[i];
}
template <class T>
matrix<T>::matrix(int n, int m , const matrix<T> &a) {
    nrows = n;
    ncols = m;
    if( nrows < a.nrows || ncols < a.ncols ) error("cannot create new matrix from old one when new matrix is smaller");
    datapoints = nrows*ncols;
    mat = new T [datapoints];
    for(int i = 0 ; i < a.nrows ; i++)
        for(int j = 0 ; j < a.ncols ; j++) {
            mat[i*ncols+j]=a.mat[i*a.ncols+j];
        }
    for(int i = a.nrows ; i < nrows ; i++)
        for(int j = 0 ; j < ncols ; j++) {
           // cout << "here we are" << endl;
            mat[i*ncols+j]=0;
        }
    for(int j = a.ncols ; j < ncols ; j++) {
        for(int i = 0 ; i < nrows; i++)
            mat[i*ncols+j]=0;
    }
    
}

template <class T>
matrix<T>::matrix(int n, int m , const matrix<T> &a, T val) {
    nrows = n;
    ncols = m;
    if( nrows < a.nrows || ncols < a.ncols ) error("cannot create new matrix from old one when new matrix is smaller");
    datapoints = nrows*ncols;
    mat = new T [datapoints];
    for(int i = 0 ; i < a.nrows ; i++)
        for(int j = 0 ; j < a.ncols ; j++) {
            mat[i*ncols+j]=a.mat[i*a.ncols+j];
        }
    for(int i = a.nrows ; i < nrows ; i++)
        for(int j = 0 ; j < ncols ; j++) {
            mat[i*ncols+j]=val;
        }
    for(int j = a.ncols ; j < ncols ; j++) {
        for(int i = 0 ; i < nrows; i++)
            mat[i*ncols+j]=val;
    }
    
}




template <class T>
inline matrix<T>::~matrix() {
	delete[] mat;
}

template <class T>
void matrix<T>::resize(int i, int j) {
    delete [] mat;
    nrows = i;
    ncols = j;
    datapoints = nrows*ncols;
    mat = new T [datapoints];
    for( int i = 0 ; i < datapoints ; i++) mat[i]=T(0);

}

template <class T>
inline T& matrix<T>::operator[](int i) {
	if ( i <0  || i > datapoints ) error("index out of range in matrix []");
	return mat[i];
}

template <class T>
inline vector1<T> matrix<T>::getrowvector(const int &i) {
        vector1<T> col(ncols);
        //if (i < 0 || i > nrows - 1 ) { cout << i << endl; cout << nrows << endl; error("index out of range in getrowvector(i)"); }


        for ( int j = 0 ; j < ncols ; j++ ) {
            col.data[j] = mat[i*ncols+j];
        }
        return col;
}

template <class T>
inline vector1<T> matrix<T>::operator()(int i , char h ) {


	if ( h == 'r' ) {
		vector1<T> col(ncols);
		if (i < 0 || i > nrows - 1 ) { cout << i << endl; cout << nrows << endl; error("index out of range in operator()(int,row)"); }


		for ( int j = 0 ; j < ncols ; j++ ) {
			col.data[j] = mat[i*ncols+j];
		}
		return col;
	}

	else if ( h == 'c' ) {
            
		vector1<T> row(nrows);
		if (i < 0 || i > ncols - 1 ) error("index out of range in operator()(int,col)");

		for ( int j = 0 ; j < nrows ; j++ ) {

			row.data[j] = mat[j*ncols+i];
		}

		return row;
	}

	else {
	error("you must select whether you are choosing a row vector or a column vector in operator ()");
	return vector1<T>();
	}
}

template <class T>
inline T& matrix<T>::operator()(int i, int j) {
	/*if ( i > nrows-1 || j > ncols-1 || i < 0 || j < 0  ) { 
            
            cout << "index out of range in matrix operator(int,int)" << endl;
            cout << "(i,j) = (" << i << "," << j <<")" << endl;
            cout << "maxvalues: " << nrows-1 << "," << ncols-1 << endl;
            error(""); }*/
	
	return mat[i*ncols+j];

}

template <class T>
inline T& matrix<T>::operator()(int i) const {
    if ( i > datapoints) error("index of out range in (i) matrix operator");
    return mat[i];
}

template <class T>
inline T& matrix<T>::gpcons(int i, int j) const {
  	// if ( i < 0 || j < 0 || i > nrows-1 || j > ncols-1  ) { 
   //          cout << "index out of range in matrix operator(int,int)" << endl;
   //          cout << "(i,j) = (" << i << "," << j <<")" << endl;
   //          cout << "maxvalues: " << nrows-1 << "," << ncols-1 << endl;

   //          error(""); }
		return mat[i*ncols+j];  
}

template <class T>
inline T& matrix<T>::periodicselect(int i, int j) {
        if( i < 0 ) {
    while ( i < 0 ) {
        i += nrows;
        if(i > 0) break;
    }
    }
    if ( i > nrows-1 ) {
    while ( i > nrows-1) {
        i -= nrows;
        if( i <= nrows -1) break;
    }
    }
    if( j < 0 ) {
    while ( j < 0 ) {
        j += ncols;
        if(j > 0) break;
    }
    }
    if ( j > ncols-1 ) {
    while ( j > ncols-1) {
        j -= ncols;
        if( j <= ncols -1) break;
    }
    }
    return mat[i*ncols+j];
}

template <class T> template <class Y>
matrix<T>::matrix(const matrix<Y> &m) {
	nrows = m.getnrows();
	ncols = m.getncols();
	datapoints = nrows*ncols;
	mat = new T [datapoints];
	for ( int i = 0 ; i < datapoints ; i++ ) mat[i] = T(m(i));
}



template <class T>
matrix<T>& matrix<T>::operator=(const matrix<T> &a) {

    delete mat;
    nrows = a.nrows;
    ncols = a.ncols;
    datapoints = a.datapoints;
    mat = new T [datapoints];
    for(int i = 0  ; i < datapoints ; i++) mat[i] = a.mat[i];
    
    

return *this;
}





template <class T>
matrix<T>& matrix<T>::operator+=(const matrix<T> &a) {
if (nrows != a.nrows || ncols != a.ncols )  error("diff sizes in operator +=");
for (int i = 0 ; i < datapoints ; i++ ) {
        mat[i]+=a.mat[i];
}
return *this;
}

template <class T>
matrix<T>& matrix<T>::operator+=(T &a) {
    for(int i = 0; i < datapoints ; i++) {
        mat[i] += a;
    }
    return *this;
}

template <class T>
matrix<T>& matrix<T>::operator-=(const matrix<T> &a) {
if (nrows != a.nrows || ncols != a.ncols) error("diff sizes in operator -=");
for (int i = 0 ; i < datapoints ; i++ )
        mat[i]-=a.mat[i];

return *this;
}

template <class T>
matrix<T>& matrix<T>::operator/=(T a) {
        for ( int i = 0 ; i < datapoints ; i++ ) {
                mat[i]/=a;
        }
        return *this;
}

template <class T>
matrix<T>& matrix<T>::operator*=(T a) {

        for ( int i = 0; i < datapoints ; i++ ) {

            mat[i] *= a;
        }
        return *this;
}

template <class T>
matrix<T>& matrix<T>::operator*=(const matrix<T> &a) {
        if (nrows != a.nrows || ncols != a.ncols) error("diff size matrics in *= operator");
        const matrix<T> mt = *this;
        T sum;

        for ( int i = 0 ; i < nrows ; i++ ) {
        for ( int k = 0 ; k < nrows ; k++ ) {
                mat[i*ncols+k]=0;
                sum = 0;
        for ( int j = 0 ; j < nrows ; j++ ) {

        sum += a.mat[i*ncols+j]*mt.mat[j*ncols+k];
        }
        mat[i*ncols+k] = sum;
        }
        }
        return *this;
}

template <class T>
void matrix<T>::inverse() {
int i,icol,irow,j,k,l,ll,q,p;
T big,dum,pivinv;
irow=0;
icol=0;
int n = nrows;

vector1<T> indxc = index<T>(n);
vector1<T> indxr = index<T>(n);
vector1<T> ipiv(n);

for ( i = 0 ; i < n ; i++ ) {
        big = 0.0;
        for ( j = 0 ; j < n ; j++ )
                if(ipiv[j] != 1)
                        for ( k = 0 ; k < n ; k++ ) {
                                if (ipiv[k] == 0 ) {
                                        if ( fabs(mat[j*ncols+k]) >= big ) {
                                                big=fabs(mat[j*ncols+k]);
                                                irow = j;
                                                icol = k;
                                        }
                                }
                        }
        ++(ipiv[icol]);
        if ( irow != icol ) {
                for ( l=0;l<n;l++ ) SWAP(mat[irow*ncols+l],mat[icol*ncols+l]);
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if ( mat[icol*ncols+icol] == 0.0 ) error("determinant zero in matrix inverse");
        pivinv=1.0/mat[icol*ncols+icol];
        mat[icol*ncols+icol]=1.0;
        for ( l= 0 ; l<n;l++ ) mat[icol*ncols+l] *= pivinv;
        for ( ll = 0 ; ll < n ; ll++ )
                if ( ll != icol ) {
                        dum=mat[ll*ncols+icol];
                        mat[ll*ncols+icol]=0.0;
                        for ( l = 0 ; l < n ; l++ ) mat[ll*ncols+l] -= mat[icol*ncols+l]*dum;
                }
}

for ( l = n-1 ; l>=0;l--) {
        if ( indxr[l] != indxc[l])
                for ( k = 0; k < n ; k++ ) {
                        q = int(indxr[l]);
                        p = int(indxc[l]);
                        SWAP(mat[k*ncols+q],mat[k*ncols+p]);
                }
}


}

template <class T>
void matrix<T>::transpose() {
        //matrix<T> c = *this;
        SWAP(nrows,ncols);
//        for ( int i = 0 ; i < nrows ; i++ ) {
//                for ( int j = 0 ; j < ncols ; j++ )
//                        mat[i*ncols+j]=c.mat[j*nrows+i];
//                }
}

template <class T>
void matrix<T>::minima(T &a) {
    a = mat[0];
    for ( int i=1 ; i < datapoints; i++) {
        if ( mat[i] < a ) {
                a = mat[i];
                }
    }
}

template <class T>
void matrix<T>::maxima(T &a) {
    a = mat[0];
    for ( int i=1 ; i < datapoints; i++) {
        if ( mat[i] > a ) {
                a = mat[i];
                }
    }
}

template <class T>
void matrix<T>::maxima(int &a, int &b) {
    a = mat[0];
    int ii = 0 ;
    int jj = 0 ;
    for(int i = 0 ; i < nrows ; i++) {
    for ( int j=0 ; j < ncols; j++) {
        if ( mat[i*ncols+j] > a ) {
                a = mat[i*ncols+j];
                ii = i;
                jj = j;
                }
    }
    }
    a = ii;
    b = jj;
}

template <class T>
void matrix<T>::minima(int &a, int &b) {
    a = mat[0];
    int ii = 0 ;
    int jj = 0 ;
    for(int i = 0 ; i < nrows ; i++) {
    for ( int j=0 ; j < ncols; j++) {
        if ( mat[i*ncols+j] < a ) {
                a = mat[i*ncols+j];
                ii = i;
                jj = j;
                }
    }
    }
    a = ii;
    b = jj;
}

template <class T>
void matrix<T>::trisolscol(vector1<T> &a, vector1<T> &b, vector1<T> c, vector1<T> &d, int j) {
    if ( j > ncols || j < 0) error("tridiagonal solve assignment is out of row range");
    if (d.size != nrows ) error("size of vector incompatible with matrix row size");

    int n = a.size;

    c.data[0] /= b.data[0];	/* Division by zero risk. */
    d.data[0] /= b.data[0];	/* Division by zero would imply a singular matrix. */
    for (int i = 1; i < n; i++){
            long double id = 1.0 / (b.data[i] - c.data[i-1] * a.data[i]);  /* Division by zero risk. */
            c.data[i] *= id;	                         /* Last value calculated is redundant. */
            d.data[i] = (d.data[i] - d.data[i-1] * a.data[i]) * id;
    }
    mat[ncols*(n-1)+j] = d.data[n-1];
    for( int i = n-2 ; i >= 0 ; i--)
        mat[ncols*i+j] = d.data[i] - c.data[i]*mat[ncols*(i+1)+j];



}

template <class T>
void matrix<T>::trisolsrow(vector1<T> &a, vector1<T> &b, vector1<T> c, vector1<T> &d, int j) {
    if ( j > nrows || j < 0) error("tridiagonal solve assignment is out of row range");
    if (d.size != ncols ) error("size of vector incompatible with matrix col size");

    int n = a.size;

    c.data[0] /= b.data[0];	/* Division by zero risk. */
    d.data[0] /= b.data[0];	/* Division by zero would imply a singular matrix. */
    for (int i = 1; i < n; i++){
            long double id = 1.0 / (b.data[i] - c.data[i-1] * a.data[i]);  /* Division by zero risk. */
            c.data[i] *= id;	                         /* Last value calculated is redundant. */
            d.data[i] = (d.data[i] - d.data[i-1] * a.data[i]) * id;
    }
    mat[ncols*j+(n-1)] = d.data[n-1];
    for( int i = n-2 ; i >= 0 ; i--)
        mat[ncols*j+i] = d.data[i] - c.data[i]*mat[ncols*j+(i+1)];



}







template <class T>
T matrix<T>::rtrap(T dr, T dz) {

    	int n = nrows;
	int m = ncols;
	long double tot = 0;
	for ( int i = 0 ; i < n ; i++ ) {
		for ( int j = 0 ; j < m ; j++ ) {
		mat[i*m+j] = mat[i*m+j]*(i*dr);
		}
	}



	tot += mat[0]+mat[m-1]+mat[(n-1)*m]+mat[(n-1)*m+(m-1)];

	for ( int i = 1 ; i < n-1 ; i++ ) {
	tot += 2*mat[i*m];
	tot += 2*mat[i*m+(m-1)];
	}
	for ( int j = 1 ; j < n-1 ; j++ ) {
	tot += 2*mat[j];
	tot += 2*mat[(n-1)*m+j];
	}
	for ( int i = 1 ; i < n-1 ; i++ ) {
		for ( int j = 1 ; j < n-1 ; j++ ) {
		tot += 4*mat[i*m+j];
		}
	}
	tot *= 0.25*dr*dz;
	return tot;
}


template <class T>
void matrix<T>::setrow(vector1<T> &v, int i) {
        int n = v.size;
        if ( n != ncols || i > nrows ) error("not same size in setting row vector, or vector range out of bounds");
        for ( int k = 0 ; k < ncols ; k++ )
                mat[ncols*i+k]=v.data[k];

}

template <class T>
void matrix<T>::setrow(T &d, int i) {

        if ( i > nrows ) error("out of bounds in trying to set constant row");
        for ( int k = 0 ; k < ncols ; k++ )
                mat[ncols*i+k]=d;
}

template <class T>
void matrix<T>::setcol(vector1<T> &v, int i) {
        int n = v.size;
        if ( n != nrows || i > ncols ) error("not same size in setting row vector, or vector range out of bounds");
        for ( int k = 0 ; k < nrows ; k++ )
                mat[ncols*k+i]=v.data[k];
}

template <class T>
matrix<T> matrix<T>::flipcols() {
    matrix<T> m(nrows,ncols);
    for(int i=0;i<nrows;i++)
        for(int j = ncols-1; j>=0;j--) {
            m(i,j)=mat[i*ncols+ncols-1-j];
        }
    return m;
}

template <class T>
matrix<T> operator*(const matrix<T> &m1,const matrix<T> &m2) {
if (m1.nrows != m2.nrows || m1.ncols != m2.ncols) error("diff size matrics in * operator");
int numrows = m1.nrows;
int numcols = m2.ncols;
// double *p1 = m1.getdata();
// double *p2 = m2.getdata();
double sum;
matrix<T> mt(numrows,numcols);
for ( int i = 0 ; i < numcols ; i++ ) {
for ( int k = 0 ; k < numcols ; k++ ) {
sum = 0;
for ( int j = 0 ; j < numcols ; ++j ) {
sum += m1.mat[i*numcols+j]*m2.mat[j*numcols+k];
}
mt.mat[i*numcols+k] = sum;
}}
// delete p1;
// delete p2;

return mt;
}

template <class T>
matrix<T> operator>>=( matrix<T> &m1, matrix<T> &m2) {
if (m1.nrows != m2.nrows || m1.ncols != m2.ncols) error("diff size matrics in * operator");
int numrows = m1.nrows;
int numcols = m2.ncols;
T *p1 = m1.getdata();
T *p2 = m2.getdata();
T sum;
matrix<T> mt(numrows,numcols);
for ( int i = 0 ; i < numcols ; i++ ) {
for ( int k = 0 ; k < numcols ; k++ ) {
sum = 0;
for ( int j = 0 ; j < numcols ; ++j ) {
sum += p1[i*numcols+j]*p2[j*numcols+k];
}
mt.mat[i*numcols+k] = sum;
}}
delete p1;
delete p2;

return mt;
}

template <class T>
vector1<T> operator*(matrix<T> &m1, vector1<T> &v) {
int n = v.getsize();
if ( n != m1.getnrows() ) error("diff size in mat * vec operation");
vector1<T> v2(n);
//T sum;
// T *p =m1.getdata();
// T *pp=v.getdat();

#pragma omp parallel for
for ( int i = 0 ; i < n ; i++ ) {
T sum = 0;
for ( int j = 0 ; j < n ; j++ ) {
sum += m1.mat[i*n+j] * v.data[j];
}
v2.data[i] = sum;
}



return v2;
}

template <class T>
matrix<T> operator*(const T a,const matrix<T> &m) {

int numrows = m.getnrows();
int numcols = m.getncols();

matrix<T> p(numrows,numcols);
int n = m.getdatapoints();
for ( int i = 0 ; i < n ; i++ ) {
        p[i] = a*m.mat[i];

}


return p;
}



// template <class T, class Q>
// matrix<T> operator*(Q a, matrix<T> &m) {

//     T ss = T(a);
// int numrows = m.getnrows();
// int numcols = m.getncols();
// T *q = m.getdata();
// matrix<T> p(numrows,numcols);
// int n = m.getdatapoints();
// for ( int i = 0 ; i < n ; i++ ) {
//         p[i] = ss*q[i];
// }
// delete q;
// return p;
// }

// template <class T, class Q>
// matrix<T> operator*(matrix<T> &m, Q a) {
//     return a*m;
// }



template <class T>
ostream& operator<<(ostream &s, const matrix<T> &a) {
int nr = a.nrows;
int nc = a.ncols;
s<< "{";
for (int i = 0; i < nr ; ++i) {
	s << "{";
	for (int j = 0; j < nc ; ++j) {
		j==nc-1 ? s <<a.mat[i*nc+j] : s <<a.mat[i*nc+j] << ",";
	}
	i != nr-1 ? s << "}\n," : s << "}";
}
s << "};" <<endl;

	return s;
}

template <class T>
T scalar(matrix<T> &a1, matrix<T> &a2) {
    T t = 0;
    int n = a1.datapoints;
    for(int i = 0 ; i < n ; i++) {
        t += a1.mat[i]*a2.mat[i];
    }
    return t;
    
    
}

template <class T>
ostream& operator<<=(ostream &s, const matrix<T> &a) {
int nr = a.nrows;
int nc = a.ncols;

for (int i = 0; i < nr ; ++i) {

	for (int j = 0; j < nc ; ++j) {
            j == nc-1 ? s << a.mat[i*nc+j] << endl : s << a.mat[i*nc+j] << ",";
        }

}


	return s;
}

template <class T>
void gnuform(ofstream &s, const matrix<T> &x, const matrix<T> &y, matrix<T> &z) {
    // x,y are the collection of gnuplot points, z is the surface height
    int nc = z.ncols;
    int nr = z.nrows;

    for(int i=0; i < nr ; i++)
        for(int j=0 ; j < nc ; j++)
    s << x(i,j) << "\t" << y(i,j) << "\t" << z(i,j) << endl;

}

template <class T>
void gnuform(ofstream &s, double a, double b, double c, double d, matrix<T> &z) {
    // given a matrix of z values, spreading from a to b in the row direction and c to d in the column direction
    // prepare an ostream object to be plottable by gnuplot surface plot
    int nc = z.ncols;
    int nr = z.nrows;

    double dr = (b-a)/nr;
    double dc = (c-d)/nc;

    for(int i=0; i < nr ; i++) {
        for(int j=0 ; j < nc ; j++) 
    s << a+i*dr << "\t" << c+j*dc << "\t" << z(i,j) << endl;
    s << endl;
        }


}


template <class T>
void gnuform_r2(ofstream &s, double dr, matrix<T> &z, double r =0) {
    int nr = z.nrows;
    int nq = z.ncols;
    double dq = (2.0*pi)/(nq);
    double x,y;

    for(int i=0; i < nr ;i++) {
        for(int j=0 ; j < nq ; j++) {
            x = (r+i*dr)*cos(dq*j);
            y = (r+i*dr)*sin(dq*j);
            s << x << "\t" << y << "\t" << z(i,j) << endl;
        }
    }



}


template <class T>
matrix<T> identitymatrix(int n) {
    matrix<T> a(n,n);
    for(int i = 0 ; i < n ; i++) {
        a(i,i) = 1;        
    }
    return a;
}

template <class T>
T sum_all_elements(const matrix<T> &a) {
    T sum = 0;
    for(int i = 0  ; i < a.getnrows() ; i++) {
        for(int j = 0 ; j < a.getncols() ; j++) {
            sum += a.gpcons(i,j);
        }
    }
    return sum;
    
}

template <class T>
void mathematica_density_r(ofstream &s, double dr, matrix<T> &z, double r=0) {
    int nr = z.nrows;
    int nq = z.ncols;
    double dq = (nq+1)/(2.0*pi);
    double x,y;
    matrix<double> res(z.datapoints,3);
    int q=0;
    for(int i=0; i < nr ; i++)
        for(int j = 0 ; j < nq ;j++) {
            x = (r + i*dr)*cos(2.0*j*dq*pi/(nq+1));
            y = (r + i*dr)*sin(2.0*j*dq*pi/(nq+1));
            
            res(q,0)=x;
            res(q,1)=y;
            res(q,2)=z(i,j);
            q++;
        }
    s << fixed << res;

}

template<class T>
matrix<T> operator-( const matrix<T> &m1) {
    int n = m1.nrows;
    int m = m1.ncols;
    matrix<T> mt(n,m);
    for(int i = 0 ; i < n*m ; i++)
        mt.mat[i]= - m1.mat[i];
    return mt;
        
}

template<class T>
matrix<T> operator-( const matrix<T> &m1,const matrix<T> &m2) {
int i;
int n = m1.nrows;
int m = m1.ncols;
if ( n != m2.nrows || m != m2.ncols ) error("diff size matrices in matrix element subtraction");
matrix<T> mt(n,m);
for (i = 0 ; i < n*m ; i++ ) {
mt.mat[i] =  m1.mat[i]-m2.mat[i];
}


return mt;
}

template<class T>
matrix<T> operator+( const matrix<T> &m1,const matrix<T> &m2) {
int i;
int n = m1.nrows;
int m = m1.ncols;
if ( n != m2.nrows || m != m2.ncols ) error("diff size matrices in matrix element addition");
matrix<T> mt(n,m);
for (i = 0 ; i < n*m ; i++ ) {
mt.mat[i] =  m1.mat[i]+m2.mat[i];
}


return mt;
}

template <class T>
matrix<T> operator&( matrix<T> &m1, matrix<T> &m2) {
int i;
int n = m1.nrows;
int m = m1.ncols;
if ( n != m2.nrows || m != m2.ncols ) {
    cout << n << "-" << m2.nrows << endl;
    cout << m << "-" << m2.ncols << endl;
    error("diff size matrices in matrix element multiplication");
}

matrix<T> mt(n,m);

for (i = 0 ; i < n*m ; i++ ) {
mt.mat[i] =  m1.mat[i]*m2.mat[i];
}

return mt;

}



template <class T, class Q>
matrix<T> operator&( matrix<T> &m1, matrix<Q> &m2) {
int i;
int n = m1.nrows;
int m = m1.ncols;
if ( n != m2.nrows || m != m2.ncols ) error("diff size matrices in matrix element multiplication");

matrix<T> mt(n,m);

for (i = 0 ; i < n*m ; i++ ) {
mt.mat[i] =  m1.mat[i]*m2.mat[i];
}

return mt;

}

template <class T>
matrix<T> operator/( matrix<T> &m1, matrix<T> &m2) {
    int n = m1.nrows;
    int m = m1.ncols;
    if ( n != m2.nrows || m != m2.ncols ) error("diff size matrices in matrix element division");
    matrix<T> mt(n,m);
    for(int i=0;i<n*m;i++) {
        mt.mat[i]=m1.mat[i]/m2.mat[i];
    }
    return mt;
}


template <class T>
void symmeterize(matrix<T> &m) {
    int i,j,n=m.getnrows();
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
        m.mat[j*m.ncols+i]=m.mat[i*m.ncols+j];
        }
    }
}

template <class T>
T rtrap(matrix<T> b, double dr, double dz,double r) {
    	int n = b.nrows;
	int m = b.ncols;
	long double tot = 0;
	for ( int i = 0 ; i < n ; i++ ) {
		for ( int j = 0 ; j < m ; j++ ) {
		b.mat[i*m+j] = b.mat[i*m+j]*(r+i*dr);
		}
	}


	tot += b.mat[0]+b.mat[m-1]+b.mat[(n-1)*m]+b.mat[(n-1)*m+(m-1)];
        //cout << tot << endl;
	for ( int i = 1 ; i < n-1 ; i++ ) {
	tot += 2*b.mat[i*m];
	tot += 2*b.mat[i*m+(m-1)];
	}
	for ( int j = 1 ; j < m-1 ; j++ ) {
	tot += 2*b.mat[j];
	tot += 2*b.mat[(n-1)*m+j];
	}
	for ( int i = 1 ; i < n-1 ; i++ ) {
		for ( int j = 1 ; j < m-1 ; j++ ) {
		tot += 4*b.mat[i*m+j];
		}
	}
	tot *= 0.25*dr*dz;
	return tot;
}


//for polar integrals weight functions are slightly different
template <class T>
T rtrap2(matrix<T> b, double dr, double dz,double r) {
    	int n = b.nrows;
	int m = b.ncols;
	long double tot = 0;
	for ( int i = 0 ; i < n ; i++ ) {
		for ( int j = 0 ; j < m ; j++ ) {
		b.mat[i*m+j] = b.mat[i*m+j]*(r+i*dr);
		}
	}

//	tot += b.mat[0]+b.mat[m-1]+b.mat[(n-1)*m]+b.mat[(n-1)*m+(m-1)];
//
//	for ( int i = 1 ; i < n ; i++ ) {
//	tot += 2*b.mat[i*m];
//	tot += 2*b.mat[i*m+(m-1)];
//	}
	for ( int j = 0 ; j < m ; j++ ) {
	tot += 2*b.mat[j];
	tot += 2*b.mat[(n-1)*m+j];
	}
	for ( int i = 1 ; i < n-1 ; i++ ) {
		for ( int j = 0 ; j < m ; j++ ) {
		tot += 4*b.mat[i*m+j];
		}
	}
	tot *= 0.25*dr*dz;
	return tot;
}

template<class T>
T rtraps(matrix<T> &b, double dr , double dz , double r) {
    //rtrap where we assume there is a zero value just beyond the integration range

    int n =  b.nrows;
    int m =  b.ncols;
    T tot = 0.0;
    
    tot += r*(b.mat[0]+b.mat[m-1]);
    
    for ( int i = 1 ; i < n ; i++ ) {
    tot += 2*(r+i*dr)*b.mat[i*m];
    tot += 2*(r+i*dr)*b.mat[i*m+(m-1)];
    }
    if(r>1E-015) {
    for ( int j = 1 ; j < m-1 ; j++ ) {
    tot += 2*(r)*b.mat[j];
    }
    }
    for ( int i = 1 ; i < n ; i++ ) {
            for ( int j = 1 ; j < m-1 ; j++ ) {
            tot += 4*(r+i*dr)*b.mat[i*m+j];
            }
    }
    tot *= 0.25*dr*dz;
    return tot;

    
}

template <class Y>
double trap(const matrix<Y> &b,double dx,double dy) {
	int n = b.nrows;
	int m = b.ncols;
	long double tot = 0;
	tot += b.mat[0]+b.mat[m-1]+b.mat[(n-1)*m]+b.mat[(n-1)*m+(m-1)];
	for ( int i = 1 ; i < n ; i++ ) {
	tot += 2*b.mat[i*m];
	tot += 2*b.mat[i*m+(m-1)];
	}
	for ( int j = 1 ; j < m ; j++ ) {
	tot += 2*b.mat[j];
	tot += 2*b.mat[(n-1)*m+j];
	}
	for ( int i = 1 ; i < n ; i++ ) {
		for ( int j = 1 ; j < m ; j++ ) {
		tot += 4*b.mat[i*m+j];
		}
	}
	tot *= 0.25*dx*dy;
	return tot;

}



template <class T>
void SQRMatrix(matrix<T> &a) {
    for(int i=0; i < a.datapoints ;i++) {
        a[i]=SQR(a[i]);
    }
}

template <class T>
matrix<T> absmatrix(const matrix<T> &a) {
    matrix<T> res(a.nrows,a.ncols);
    for(int i = 0  ; i < a.datapoints ; i++)
        res.mat[i]=abs(a.mat[i]);
    return res;    
}

template <class T>
matrix<T> posmatrix(const matrix<T> &a) {
    matrix<T> res(a.nrows,a.ncols);
    for(int i = 0  ; i < a.datapoints ; i++) {
        if(a.mat[i] < 0 ) res.mat[i]= 0;
        else res.mat[i] = a.mat[i];
        }
    return res;    
}

template <class T>
matrix<T> logmatrix(const matrix<T> &a) {
    matrix<T> res(a.nrows,a.ncols);
    for(int i = 0  ; i < a.datapoints ; i++)
        res[i]=log(a(i));
    return res;
}

template <class T>
matrix<T> expmatrix(const matrix<T> &a) {
    matrix<T> res(a.nrows,a.ncols);
    for(int i = 0  ; i < a.datapoints ; i++)
        res[i]=exp(a.mat[i]);
    return res;
}

template <class T>
matrix<T> mirror(const matrix<T> &a) {
    int n = a.nrows;
    int m = a.ncols;
    matrix<T> res(n,2*m);
    for(int i = 0  ; i < n ; i++) {
        for(int j = 0 ; j < m ; j++) {
            res(i,j+m)=a.gpcons(i,j);
            res(i,j)=a.gpcons(i,m-j-1);
        }
    }
    
    return res;
}

template <class T>
matrix<T> antimirror(const matrix<T> &a) {
    int n = a.nrows;
    int m = a.ncols;
    matrix<T> res(n,2*m);
    for(int i = 0  ; i < n ; i++) {
        for(int j = 0 ; j < m ; j++) {
            res(i,j+m)=a.gpcons(i,j);
            res(i,j)=-a.gpcons(i,m-j-1);
        }
    }
    
    return res;
}

template <class T>
matrix<T> mirrorc(const matrix<T> &a) {
	int n = a.nrows;
	int m = a.ncols;
    if ( m % 2 == 0 ) error("mirrorc needs an odd matrix");
	matrix<T> res(n,2*m-1);
	for(int i = 0 ; i < n ; i++ ) {
		for(int j = 0 ; j < 2*m-1 ; j++ ) {
		if ( j < m ) res(i,j) = a.gpcons(i,m-1-j);
                else res(i,j) = a.gpcons(i,j-(m-1));
		}
	}
	return res;
}

template <class T>
matrix<T> antimirrorc(const matrix<T> &a) {
	int n = a.nrows;
	int m = a.ncols;
    if ( m % 2 == 0 ) error("mirrorc needs an odd matrix");
	matrix<T> res(n,2*m-1);
	for(int i = 0 ; i < n ; i++ ) {
		for(int j = 0 ; j < 2*m-1 ; j++ ) {
		if ( j < m ) res(i,j) = -a.gpcons(i,m-1-j);
		else res(i,j) = a.gpcons(i,j-(m-1));
		}
	}
	return res;
}

template <class T>
vector1<T> con_vec(matrix<T> &a) {
    vector1<T> u(a.datapoints);
    for(int i=0;i<a.datapoints;i++)
        u.data[i]=a.mat[i];
    return u;
}

template <class T>
T variance(const matrix<T> &m1, const matrix<T> &m2) {
    if(m1.nrows != m2. nrows || m1.ncols != m2.ncols ) error("matrices need to be of same dimension in variance");
    T av = 0;
    for(int i = 0 ; i < m1.datapoints ;i++) {
            av+=SQR(m1.mat[i]-m2.mat[i]);
    }
    
    return sqrt(av);
}

template <class T>
void ludcmp(matrix<T> &a, vector1<int> &indx, T &d) {
const double TINY = 1.0E-20;
int i,imax=0,j,k;
double big,dum,sum,temp;
int n = a.nrows;
int m = a.ncols;
vector1<T> vv(n);
d=1.0;
for(i=0;i<n;i++) {
	big=0.0;
	for(j=0;j<n;j++) {
		temp=fabs(a.mat[i*m+j]);
		if (temp >big ) big = temp;
	}
	if ( big == 0.0 ) error("Singular matrix in routine ludcmp");
	vv.data[i] = 1.0/big;
}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a.mat[i*m+j];
			for (k=0;k<i;k++) sum -= a.mat[i*m+k] + a.mat[k*m+j];
			a.mat[i*m+j]=sum;
		}
		big = 0.0;
		for ( i = j ; i < n ; i++) {
			sum = a.mat[i*m+j];
			for(k=0;k<j;k++) sum -= a.mat[i*m+k] + a.mat[k*m+j];
			a.mat[i*m+j]=sum;
			if ((dum=vv.data[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if ( j != imax ) {
			for(k=0;k<n;k++){
				dum=a.mat[imax*m+k];
				a.mat[imax*m+k]=a.mat[j*m+k];
				a[j*m+k] = dum;
			}
			d=-d;
			vv.data[imax]=vv.data[j];
		}
		indx.data[j]=imax;
		if ( a.mat[j*m+j] == 0.0 ) a.mat[j*m+j] = TINY;
		if ( j != n-1 ) {
			dum = 1.0/(a.mat[j*m+j]);
			for ( i = j+1 ; i < n ; i++ ) a.mat[i*m+j] *= dum;
		}
	}
}

template <class T>
void lubksb(matrix<T> &a, vector1<int> &indx , vector1<T> &b) {
    int i,ii=0,ip,j;
    double sum;
    int n = a.nrows;
    for(i=0;i<n;i++) {
        ip =indx[i];
        sum =b[ip];
        b[ip]=b[i];
        if(ii != 0)
            for (j=ii-1;j<i;j++) sum -=a(i,j)*b[j];
        else if ( sum != 0)
            ii = i+1;
        b[i]=sum;
    }
    for(i=n-1;i>=0;i--) {
        sum = b[i];
        for(j=i+1;j<n;j++) sum -= a(i,j)*b[j];
        b[i]=sum/a(i,i);
    }
}

template <class T>
void qrdcmp(matrix<T> &a, vector1<T> &c, vector1<T> &d, bool &sing) {
int i,j,k;
double scale,sigma,sum,tau;
int n = a.nrows;
int m = a.ncols;
sing=false;
for(k=0;k<n-1;k++){

scale = 0.0;
for(i=k;i<n;i++){ scale = MAX(scale,fabs(a.mat[i*m+k])); }
if (scale == 0.0) {
	sing=true;
	c.data[k]=d.data[k]=0;
}
else {

	for(i=k;i<n;i++)  a.mat[i*m+k] /= scale;
	for( sum = 0.0,i=k;i<n;i++) sum += SQR(a.mat[i*m+k]);
	sigma=SIGN(sqrt(sum),a.mat[k*m+k]);
	a.mat[k*m+k]+=sigma;
	c.data[k]=sigma*a.mat[k*m+k];
	d.data[k]=-scale*sigma;
	for(j=k+1;j<n;j++) {
		for(sum=0.0,i=k;i<n;i++) sum += a.mat[i*m+k]*a.mat[i*m+j];
		tau=sum/c.data[k];
		for(i=k;i<n;i++) a(i,j) -= tau*a.mat[i*m+k];
	}
}
}
d.data[n-1]=a.mat[(n-1)*m+(n-1)];
if(d.data[n-1] == 0.0 ) sing =true;
}

template <class T>
void rsolv(matrix<T> &a,vector1<T> &d, vector1<T> &b)
{
	int i,j;
	double sum;
	int n = a.nrows;
	int m = a.ncols;
	b.data[n-1] /= d.data[n-1];
	for ( i=n-2;i>=0;i--) {
		for ( sum = 0.0, j=i+1;j<n;j++) sum += a.mat[i*m+j]*b.data[j];
		b[i] = (b[i]-sum)/d[i];
	}
}

template <class T>
void rotate(matrix<T> &r, matrix<T> &qt, const int i, const T a, const T b) {
	int j;
	double c,fact,s,w,y;
	int n = r.nrows;
	int m = r.ncols;
	c=0.0;
	if ( a==0.0) {
		c=0.0;
		s=(b>=0.0 ? 1.0 : -1.0 );
	}
	else if ( fabs(a) > fabs(b)) {
		fact = b/a;
		c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
		s=fact*c;
	}
	else {
		fact = a/b;
		s= SIGN(1.0/sqrt(1.0+(fact*fact)),b);
		c=fact*s;
	}
	for ( j = i; j < n ; j++ ) {
		y=r.mat[i*m+j];
		w=r.mat[(i+1)*m+j];
		r.mat[i*m+j]=c*y-s*w;
		r.mat[(i+1)*m+j]=s*y+c*w;
	}
	for ( j = 0 ; j < n ; j++ ) {
		y = qt.mat[i*m+j];
		w = qt.mat[(i+1)*m+j];
		qt.mat[i*m+j]=c*y-s*w;
		qt.mat[(i+1)*m+j] = s*y+c*w;
	}
}

template <class T>
void qrupdt(matrix<T> &r, matrix<T> &qt, vector1<T> &u, vector1<T> &v) {
	int i,k;
	int n = u.size;
	int m = r.ncols;
	for ( k =n-1;k>=0;k--)
		if ( u[k] != 0.0) break;
	if ( k < 0 ) k =0;
	for ( i = k-1 ; i >= 0.0;i--) {
		rotate(r,qt,i,u.data[i],-u.data[i+1]);
		if (u.data[i] == 0.0)
			u.data[i]=fabs(u.data[i+1]);
		else if ( fabs(u.data[i]) > fabs(u[i+1]))
			u.data[i] = fabs(u.data[i])*sqrt(1.0+SQR(u.data[i+1]/u.data[i]));
		else u.data[i] = fabs(u.data[i+1])*sqrt(1.0+SQR(u.data[i]/u.data[i+1]));
	}
	for ( i = 0 ; i < n ; i++ ) r.mat[i] += u.data[0]*v.data[i];
	for ( i = 0 ; i < k ; i++ )
		rotate(r,qt,i,r.mat[i*m+i],-r.mat[(i+1)*m+i]);
}

template <class T>
T maximum( const matrix<T> &b) {
    T a = b.mat[0];
    for ( int i=1 ; i < b.datapoints; i++) {
        if ( b.mat[i] > a ) {
                a = b.mat[i];
                }
        
        }
    return a;
    
}

template <class T>
T minimum( const matrix<T> &b) {
    T a = b.mat[0];
    for ( int i=1 ; i < b.datapoints; i++) {
        if ( b.mat[i] < a ) {
                a = b.mat[i];
                }
        
        }
    return a;
    
}

template <class T>
void outfunc(const matrix<T> &a, string s) {
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if(!myfileg.is_open()) error("failed to open file");
    myfileg <<= a;
    myfileg.close();
}


template <class T>
void outfunc2(const matrix<T> &a, string s) {
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str(),ios::app);
    if(!myfileg.is_open()) error("failed to open file");
    myfileg <<= a;
    myfileg.close();
}

template <class T>
void outfunc(const vector1<T> &a, string s) {
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if(!myfileg.is_open()) error("failed to open file");
    myfileg <<= a;
    myfileg.close();
    
}

template <class T>
matrix<T> importcsv(string filename, T temp, bool &error) {
    ifstream data(filename.c_str());
    if(!data.is_open()) {cout << "file not found during import" << endl; error= true; return matrix<long double>(100,100);  }
    string line;
    vector<vector<T> > a;
    while(getline(data,line))
    {
        
        stringstream  lineStream(line);
        string        cell;
        vector<T> b;
        while(getline(lineStream,cell,',')) { 
            b.push_back(atof(cell.c_str()));
            
        }
        a.push_back(b);
    }
    
    int nr = a.size();
    int nz = a[0].size();
    
    matrix<T> qq(nr,nz);
    for(int i = 0 ; i < nr ; i++) {
        for(int j = 0 ; j < nz ; j++) {
            qq(i,j)=a[i][j];
        }
    }
    error = false;
    return qq;
    
}

matrix<double> periodic(matrix<double> &a, matrix<double> &b, double l) {
if(a.getnrows() != b.getnrows() ) error("size of matrices must be the same");

matrix<double> res(a.getnrows(),3);
for(int i = 0 ; i < a.getnrows() ; i++ ) {
double dx  = a(i,0)-b(i,0);
if(abs(dx)>l/2) dx = dx-SIGN(l,dx);
double dy  = a(i,1) - b(i,1);
if(abs(dy)>l/2) dy = dy-SIGN(l,dy);
double dz  = a(i,2) - b(i,2);
if(abs(dz)>l/2) dz = dz-SIGN(l,dz);
res(i,0) = dx;
res(i,1) = dy;
res(i,2) = dz;

}
return res;

}

vector1<double> distancedistribution(matrix<double> &a, matrix<double> &b, double l) {
matrix<double> res = periodic(a,b,l);
vector1<double> disSQR(res.getnrows());
for(int i = 0 ; i < res.getnrows() ; i++ ) {
disSQR[i] = SQR(res(i,0))+SQR(res(i,1))+SQR(res(i,2));
}
return disSQR;
}

#endif

