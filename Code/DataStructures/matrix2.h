/* 
 * File:   matrix.h
 * Author: Dino
 *
 * Created on 28 November 2010, 10:11
 */

/* A matrix class, despite the fact that there are many matrix classes already avaliable, I thought it would
 be easier to create the class myself, this also makes checking for errors easier. Inspiration was taken from
 BH Flowers book "Numerical Computing in C++ and also from Numerical Recipes in C++ (Press et al) */

#ifndef MATRIX_H
#define	MATRIX_H

#include "vector1.h"
#include "vector1.cpp"

#include "Complex.h"

//we'd like these classes below to have access to the matrix class
struct cube;
class MD;


template <class T>
class matrix  {
friend struct cube;
friend class MD;

private:
	int ncols; // no. of cols
	int nrows; // no . of rows
	T *mat; //where the data is stored
	int datapoints; //the total number of data nrows*ncols
public:
	matrix(); // null constructor
        matrix(int);
	matrix(int, int); // constructor with int rows and int columns
	matrix(int, int, T); // int rows and columns, with elements initialized to T

	matrix(const matrix&); //copy constructor
        matrix(int,int,const matrix&); //imbed small matrix in large one, with rest of values 0 
        matrix(int,int,const matrix&,T); //imbed small matrix in large one, with rest of values T 
        template <class Y>
        matrix(const matrix<Y>&);

        matrix(int,int,int, T[]);
        matrix(int, int, vector1<T>&);

	~matrix(); // destructor

	inline T& operator[](int); //return the data value loacted at the data point
        inline T& operator()(int) const;
	inline vector1<T> operator()(int,char); //return the int row (char = 'r' ) or column (char = 'c')
	inline T& operator()(int,int); // return the int,int element
    inline T& gpcons(int,int) const;
    inline vector1<T> getrowvector(const int&);
    inline T& periodicselect(int,int); //periodic select along one axis
        
        matrix& operator+() { return *this; } //unary plus

	matrix& operator=(const matrix&); //assignment operator
        
        matrix& operator=(const matrix&) const;
	matrix& operator+=(const matrix&); // add a matrix to another
        matrix& operator+=(T&);
	matrix& operator-=(const matrix&); //minus a matrix from another
	matrix& operator/=(T); //divide a matrix by T
	matrix& operator*=(T); //multiply a matrix by T
	matrix& operator*=(const matrix&); // multiply the matrix by another

    matrix* clone() const { return new matrix(*this); } 

	int getdatapoints() { return datapoints; }
    int getdatapoints() const { return datapoints; }
    int getnrows() { return nrows; }
    int getnrows() const { return nrows; }
    int getNsafe() {return nrows; }
	int getncols() { return ncols; }
    int getncols() const { return ncols; }
    void resize(int,int); //resize matrix to int,int

   // void removeoverlaps(int attempts) //remove overlaps greater than size d in the data
    
    T* getdata() { 
        T *y = new T [datapoints];
        for(int i= 0 ; i < datapoints ; i++ )
            y[i]=mat[i];
        return y;
    }
    T* getdata() const {
       // cout << "here" << endl;  
        T *y = new T [datapoints];
        for(int i= 0 ; i < datapoints ; i++ )
            y[i]=mat[i];
       // cout << "ok" << endl;
        return y;          
    }

	void inverse(); //inverse matrix
	void transpose(); //transpose matrix
        T rtrap(T dr, T dz);
        void minima(T&);
        void maxima(T&);
	void minima(int&,int&);        
	void maxima(int&,int&);
        void trisolscol(vector1<T>&,vector1<T>&,vector1<T>,vector1<T>&,int);
        void trisolsrow(vector1<T>&,vector1<T>&,vector1<T>,vector1<T>&,int);
        void cylsolscol(const vector1<T>&, const vector1<T>, const vector1<T>&, const T&, const T&, const vector1<T>&,int);
        void cylsolsrow(const vector1<T>&, const vector1<T>, const vector1<T>&, const T&, const T&, const vector1<T>&,int);        
        void inversetri(vector1<T>&,vector1<T>&,vector1<T>);
	void setrow(vector1<T>&,int); // set the int row to be  vector
        void setrow(T&,int); // set the int row to a be a constant T
	void setcol(vector1<T>&,int); // set the int column to be vector
        matrix<T> flipcols();
        
        template <class Y>
        friend Y scalar(matrix<Y>&, matrix<Y>&);

        template <class Y>
	friend void ludcmp(matrix<Y>&,vector1<int>&,Y&); //LU Decompse matrix

        // LU Solve a matrix equation where the first matrix is the LU Decompisition, and vector is the vector passed
        // to the ludcmp
        template <class Y>
        friend void lubksb(matrix<Y>&, vector1<int>&,vector1<Y>&);

        // QR decompose a matrix
        template <class Y>
	friend void qrdcmp(matrix<Y>&,vector1<Y>&,vector1<Y>&,bool&);

        // Solve for a QR decomposed matrix
        template <class Y>
	friend void rsolv(matrix<Y>&,vector1<Y>&,vector1<Y>&);

        template <class Y>
	friend void rotate(matrix<Y>&,matrix<Y>&,const int, const Y, const Y);
	
        template <class Y>
        friend void qrupdt(matrix<Y>&,matrix<Y>&,vector1<Y>&,vector1<Y>&);

        // matrix multiplication
        template <class Y>
	friend matrix<Y> operator*(const matrix<Y> &m1,const matrix<Y> &m2);
        
        template <class Y>
        friend matrix<Y> operator>>=(matrix<Y> &m1,matrix<Y> &m2);

        //matrix vector multiplication
        // template <class Y>
        // friend vector1<Y> operator*(matrix<Y> &m1, vector1<Y> &v1);

        template <class Q, class Y>
        friend vector1<Q> operator*(matrix<Y> &m1, vector1<Q> &v1);

        template <class Y>
        friend matrix<Y> operator-(const matrix<Y>&);
        
        //matrix subtraction
        template <class Y>
        friend matrix<Y> operator-(const matrix<Y> &m1, const matrix<Y> &m2);
        
        template<class Y>
        friend matrix<Y> operator+(const matrix<Y> &m1,const matrix<Y> &m2);
        //multiplication by a scalar of a different type
    
        // template <class Y, class Q>
        // friend matrix<Y> operator*(Q, matrix<Y> &m1 );

        // template <class Y, class Q>
        // friend matrix<Y> operator*(matrix<Y> &m1, Q);

        // multiply matrix by a constant
        template <class Y>
	friend matrix<Y> operator*(const Y,const matrix<Y> &m1);

        template <class Q, class Y>
        friend matrix<Y> operator*(const Q, const matrix<Y> &m1);

        // return a matrix whose values are equal to the multiplication of the correspong multiplication
        // of the same (row index, column index) of the two matrices m1 and m2
        template <class Y>
	friend matrix<Y> operator&(matrix<Y> &m1, matrix<Y> &m2);

        template <class Y, class Q>
        friend matrix<Y> operator&(matrix<Y> &m1, matrix<Q> &m2);

        //matrix division of one matrix by another (i.e the elements of matrix m1 divided by corresponding element
        // of m2.
        template <class Y>
        friend matrix<Y> operator/(matrix<Y> &m1, matrix<Y> &m2);


        // output stream for matrix (mathematica formatted)
        template <class Y>
	friend ostream& operator<<(ostream&, const matrix<Y>&);

        template <class Y>
	friend ostream& operator<<=(ostream&, const matrix<Y>&);


        // output in gnu form
        template <class Y>
        friend void gnuform(ofstream&,const matrix<Y>&,const matrix<Y>&, matrix<Y>&);

        // other output in gnuform
        template <class Y>
        friend void gnuform(ofstream&, double a, double b, double c, double d, matrix<Y>&);

        //gnuform for a radial plot
        template <class Y>
        friend void gnuform_r(ofstream&, double dr, matrix<Y>&,double r);

        template <class Y>
        friend void gnuform_r2(ofstream&, double dr, matrix<Y>&,double r);

        template <class Y>
        friend void mathematica_density_r(ofstream&, double dr, matrix<Y>&, double r);

        //identity matrix
        template <class Y>
        friend matrix<Y> identitymatrix(int);
        
        template <class Y>
        friend Y sum_all_elements(const matrix<Y>&);

        // 2D trapezium rule for a set of data points matrix, and seperation double, double
        template <class Y>
	friend double trap(const matrix<Y>&, double, double);

        // 2D radial trapezium rule for a set of data points, and seperation double,double
        template <class Y>
	friend Y rtrap(matrix<Y>, double, double,double r);

        template <class Y>
        friend Y rtrap2(matrix<Y>, double, double, double r);
        
        template <class Y>
        friend Y rtraps(matrix<Y>&,double,double,double r);

        template<class Y, class H>
        friend H qromb2(Y &func, H a, H b, int decdigs);        

        template<class Y, class H>
        friend H qromo2(Y &func, H a, H b, int decdigs);
        //convert a matrix to a vector
        template <class Y>
        friend vector1<Y> con_vec(matrix<Y>&);

        template <class Y>
        friend void symmeterize(matrix<Y>&); // make a matrix symmetric

        template <class Y>
        friend void SQRMatrix(matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> absmatrix(const matrix<Y>&);

        template <class Y>
        friend matrix<Y> posmatrix(const matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> logmatrix(const matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> expmatrix(const matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> mirror(const matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> antimirror(const matrix<Y>&);    

        template <class Y>
        friend matrix<Y> mirrorc(const matrix<Y>&);

        template <class Y>
        friend matrix<Y> antimirrorc(const matrix<Y>&);
		
        template <class Y>
        friend Y variance(const matrix<Y>&, const matrix<Y>&);
        
        template <class Y>
        friend Y maximum(const matrix<Y>&);
        template <class Y>
        friend Y minimum(const matrix<Y>&);
        
        template <class Y>
        friend matrix<Y> convertmatrix(matrix<Y>&a,int,int);

        template <class Y>
        friend matrix<Y> poly_convertmatrix(matrix<Y>&a,int,int);

        template <class Y>
        friend void ft2D(matrix<Y> &m, matrix<Y> &reout, matrix<Y> &imout);

        friend int FFT2D(matrix<COMPLEX> &c, int nx, int ny, int dir);

        friend matrix<double> realmatrix(matrix<complex<double> > &a);
        friend matrix<double> imaginarymatrix(matrix<complex<double> > &a);
        friend matrix<complex<double> > combinematrix(matrix<double> &re, matrix<double> &im);
};



#endif	/* MATRIX_H */




