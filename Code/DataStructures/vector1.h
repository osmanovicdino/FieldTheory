/* 
 * File:   vector1.h
 * Author: Dino
 *
 * Created on 28 November 2010, 09:16
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>
using namespace std;


#ifndef VECTOR1_H
#define	VECTOR1_H


template <class T>
class matrix;

template <class T>
class s_matrix;

template <class T>
class d_matrix;

typedef double (*enfunc)(double);
typedef double (*func)(double, double);
typedef double (*func2)(double,double,double);




//vector class, supports basic geometric functions, as well as being useful for storing data
//ideas and functions in class developed by author as well as taken from numerous sources
template <class T>
class vector1 {
friend class s_matrix<T>;
friend class d_matrix<T>;
friend class matrix<T>;
private:
	int size; //size of the vector
	T * data; //pointer to data members


public:
	vector1(); //null constructor
	vector1(int); //create a vector of int elements
	vector1(int,const T*); //create a vector of int elements, with element values in the pointer
	vector1(int, T[]); //create vector with elements in array
	vector1(int,T); //create a vector initialized to the value T
	vector1(const vector1&);
	vector1(const vector1&,int); //vector1& with n extra elements
	~vector1(); //destructor

	void delete_element(int); //delete an element;
	inline T& operator[](int i) {return data[i];}
	inline T& operator()(int i) {return data[i];} // equivalent to above
    inline T& gpcons(int) const;
    inline T& periodicselect(int); //move boundry conditions periodically.
    vector1& operator+(); // unary plus
	vector1& operator=(const vector1&); // assignment operator
	vector1& operator+=(const vector1&); // add another vector
	vector1& operator-=(const vector1&); // minus another vector
	vector1& operator+=(T); // add a constant Te to all elements
	vector1& operator-=(T); // same as above but minus
	vector1& operator*=(T); // multiply by a factor T
        
        template<class Y>
        vector1& operator*=(Y x)
        {
                for (int i = 0; i < size; ++i)
                {
                        data[i] *= x;
                }
                return *this;
        }
    //vector1& operator*=(matrix<T>&);
	vector1& operator/=(T); // divide by a double
        bool operator==(const vector1&);


        template <class Y>
	   friend vector1<Y> operator-(const vector1<Y>&); // unary minus a = -a

        template <class Y>
        friend vector1<Y> operator+(const vector1<Y>& , const vector1<Y>&); // add 2 vectors

        template <class Y>
        friend vector1<Y> operator-(const vector1<Y>&, const vector1<Y>&); // minus two vectors

        template <class Y>
        friend vector1<Y> operator*(const vector1<Y>&, Y); // multiply a vector by a constant

        template <class Y, class Q>
        friend vector1<Y> operator*(const vector1<Y>&, Q); //multiply a vector by a value Q (not necessarily the same type)

        template <class Y>
        friend vector1<Y> operator*(Y, const vector1<Y>&); // reverse
	        
        template <class Y>
        friend vector1<Y> operator*(matrix<Y> &m1, vector1<Y> &v1);

        template <class Y>
        friend vector1<Y> operator/(const vector1<Y> &a, const vector1<Y> &b);
	
        template <class Y>         
        friend vector1<Y> operator&(const vector1<Y>&, const vector1<Y>&); // multiply each element in a vector by the corresponding element in another vector

        template <class Y>
        friend bool operator==(const vector1<Y>&, const vector1<Y>&); //tests


	int getsize() const { return size; }
	T* getdat() {return data; }
	void setval(T a); //set all the values in the vector to a
	void swap(int,int); // swap
        void resize(int); // resize a vector to size int
        void resize_keep(int); //keep the first int elements of the vector

	void trisols(vector1&,vector1&,vector1,vector1&); // solution to a tridiagonal matrix equation with 3 vectors being the diagonals and one the solution
        void trisols(const vector1&,const vector1&, const vector1, const vector1&); // overloaded function
        void trimult(vector1&, vector1&, vector1&, vector1&); // A * vector where A is tridiagonal matrix
        void cyclic(const vector1&, const vector1, const vector1&, const T&, const T&, const vector1&); // solves tridiagonal where there are elements in the corners
        void cyclicmult(vector1&, vector1&, vector1&,T alpha, T beta, vector1&); // multiply a vector by a cyclic matrix
          
        template <class Y>
        friend void cyclic(const vector1<Y>&, const vector1<Y>&, const vector1<Y>&, const Y&, const Y&, const vector1<Y>&, vector1<Y>&); // solves tridiagonal where there are elements in the corners
      
        
        template <class Y>
        friend double scalar(const vector1<Y>&, const vector1<Y>&); //scalar product

        template <class Y>
        friend vector1<Y> unitvector(const vector1<Y>&, const vector1<Y>&);

        template <class Y>
        friend vector1<Y> unitvector(const vector1<Y>&, const vector1<Y>&,double&);

        template <class Y>
	    friend double norm(const vector1<Y>&, const vector1<Y>&); // works out distances
	
        template <class Y>
        friend double norm2(const vector1<Y>&); // magnitude

        template <class Y>
        friend ostream& operator<<=(ostream&,const vector1<Y>&);
        
        template <class Y>
	   friend ostream& operator<<(ostream&, const vector1<Y>&); //output (formatted for mathematica)
	
        template <class Y>
        friend istream& operator>>(istream&, vector1<Y>&); //input (rarely used)

        template< class Y>
        friend void print_two_vectors(ofstream&, const vector1<Y>&, const vector1<Y>&);

        // assign a int row for a matrix, with elements of vector
        template <class Y>
	   friend void setrow(vector1<Y>&,int);

        template <class Y>
	friend void setcol(vector1<Y>&,int); // assignments of rows, columns

        template <class Y>
	friend Y minval(const vector1<Y>&); // find minimum value of vector

        template <class Y>
	friend Y maxval(const vector1<Y>&); // find maximum value of vector

	template <class Y>
	friend Y minvalabs(const vector1<Y>&); // find minimum value of vector

        template <class Y>
	friend Y maxvalabs(const vector1<Y>&); // find maximum value of vector

	template <class Y>
	friend Y elementmultiply(const vector1<Y>&);

        template <class Y>
	friend int minindex(const vector1<Y>&); // find index of minimum

        template <class Y>
	friend int maxindex(const vector1<Y>&); // index of maximum

        // from a vector, return the maximum value, as well as the index where it occurs.
        template <class Y>
	friend void maxima(vector1<Y>&,T&,int&);


        template <class Y>
	friend vector1 symmetrize(vector1<Y>&,vector1<Y>&,vector1<Y>&); // symmetrize a unsymmetrix tridiagonal matrix


	template <class Y>
	friend Y meanish(const vector1<Y>&);
        // trapezium rule, perform integration of discrete set of data
        template <class Y>
	friend long double trap(const vector1<Y>&,double);

        // perform a radial integral, starting from 0, with spacing double;
        template <class Y>
	friend long double rtrap(const vector1<Y>&,double);

        
        // template<class Y, class H>
        // friend H qromb2(Y &func, H a, H b, int decdigs = 10);
        
        // template <class Y, class H>
        // friend H qromb(Y &func, H a, H b, const H eps =  1.0e-10);

        // return the vector with each element logarithmed
        template <class Y>
	friend vector1<Y> logvector(vector1<Y>&);

        // return the vector with each element exponentiated
        template <class Y>
	friend vector1<Y> expvector(vector1<Y>&);

        // return a vector of length n with elements (1,2,3,...,n)
        template <class Y>
        friend vector1<Y> index(int &n);

        // return the vector with each element squared
        template <class Y>
        friend vector1<Y> SQRvector(vector1<Y>);

        // compare two vectors and return the number of elements which are not the same;
        template <class Y>
        friend int compare(vector1<Y>&,vector1<Y>&);

        template <class Y>
        friend void trisolscol(vector1<T>&,vector1<T>&,vector1<T>,vector1<T>&,int);

        template <class Y>
        friend void trisolsrow(vector1<T>&,vector1<T>&,vector1<T>,vector1<T>&,int);

        template <class Y>
        friend vector1<Y> vectovec(std::vector<Y>);


        template <class Y>
        friend int compare(vector1<Y>&,vector1<Y>&,Y);
        
        template <class Y, class P>
        friend vector1<P> writefunction(Y &func,P,P,int);
        
        template <class Y>
        friend vector1<Y> complexmultiply(vector1<Y>&,vector1<Y>&);
        
        template <class Y>
        friend void four1(vector1<Y>&,const int);
        template <class Y>
        friend void realft(vector1<Y>&,const int);
        
        template <class Y>
        friend vector1<Y> smooth(const vector1<Y>&, int, bool);
};



#endif	/* VECTOR1_H */






