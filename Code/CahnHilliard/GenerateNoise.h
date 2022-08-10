#ifndef GENERATENOISE_H
#define GENERATENOISE_H
//Routines to Generate Noise in K-Space
#include "chbuilder.h"

struct modelAnoise {
    double prefactor;
    modelAnoise(double prefactorr) : prefactor(prefactorr) {}
    double operator()(double k1, double k2)
    {
        return prefactor;
    }
};

struct modelBnoise
{
    double prefactor;
    modelBnoise(double prefactorr) : prefactor(prefactorr) {}
    double operator()(double k1, double k2)
    {
        return prefactor*(k1 * k1 + k2 * k2);
    }
};

struct modelABnoise
{
    double prefactor1;
    double prefactor2;
    modelABnoise(double prefactorr1, double prefactorr2) : prefactor1(prefactorr1), prefactor2(prefactorr2) {}
    double operator()(double k1, double k2)
    {
        return prefactor1 * (k1 * k1 + k2 * k2) + prefactor2;
    }
};



template<class T>
struct GenNoise {
    CH_builder params;
    T **random_field;

    GenNoise();
    GenNoise(const CH_builder &p);
    GenNoise(const GenNoise<T> &a);
    GenNoise &operator=(const GenNoise<T> &); // assignment operator

    ~GenNoise();

    template<class Q>
    void GenFields(Q &func, vector1<double> &strs, double samprat); //create the random fields in k space, with psd of func, and size params;
};

#include "GenerateNoise.cpp"

#endif /* GENERATENOISE_H */
