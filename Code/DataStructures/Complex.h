#ifndef COMPLEX_H
#define COMPLEX_H

struct COMPLEX
{
    double real;
    double imag;
    COMPLEX()
    {
        real = 0.0;
        imag = 0.0;
    }

    template <class T>
    COMPLEX(T a)
    {
        real = (double)a;
        imag = 0.0;
    }

    friend ostream &operator<<(ostream &s, const COMPLEX &a)
    {
        s << "(" << a.real << "," << a.imag << "i)";
        return s;
    }
    double getreal() {
        return this->real;
    }
    COMPLEX& operator+=(const COMPLEX &m1) {
        real += m1.real;
        imag += m1.imag;
        return *this;
    }

    friend COMPLEX operator+(const COMPLEX &m1, const COMPLEX &m2)
    {
        COMPLEX res;
        res.real = m1.real + m2.real;
        res.imag = m1.imag + m2.imag;
        return res;
    }

    friend COMPLEX operator*(const COMPLEX &m1, const COMPLEX &m2) {
       COMPLEX res;
       res.real = m1.real*m2.real - m1.imag*m2.imag;
       res.imag = m1.imag*m2.real + m1.real*m2.imag;
       return res; 
    }
    friend COMPLEX operator*(double m2, const COMPLEX &m1)
    {
        COMPLEX res;
        res.real = m1.real * m2;
        res.imag = m1.imag * m2;
        return res;
    }

    COMPLEX& operator=(const COMPLEX &v)
    { //assigment operator

        real = v.real;
        imag = v.imag;
        return *this;
    }
};

#endif /* COMPLEX_H */
