#ifndef VEC_VEC_H
#define	VEC_VEC_H
#include "vector1.h"
#include "vector1.cpp"

#include "matrix.h"
#include "matrix.cpp"

template <class T>
class vec_vec {
	private:
		int n; // the size of the array of vectors;
		vector<vector1<T> > *dat;
	public:
		vec_vec();
		vec_vec(int);
		vec_vec(int,int);
		vec_vec(const vec_vec&);
		vec_vec(const matrix<T>&);


		~vec_vec();

		int getn() { return n; }

		inline vector1<T>& operator[](int);
		inline T& operator()(int,int);
		
		vec_vec& operator=(const vec_vec&);

		void remove(const int&);
		void setvalue(const int&, const vector1<T>&);
		void append(const vector1<T>&);

		T max(); //return the max value in vec_vec


		template <class Y>
		friend ostream& operator<<(ostream&, const vec_vec<Y>&);
		template <class Y>
		friend ostream& operator<<=(ostream&, const vec_vec<Y>&);


};
#include "vec_vec.cpp"

#endif