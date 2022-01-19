
#ifndef CAHNHILLIARD_H
#define CAHNHILLIARD_H

#include "chbuilder.h"
#include "utilityOutputs.h"
#include "fftw3.h"
#include "weights.h"
#include "chemistry.h"
#include "GenerateNoise.h"
#include "updateRules.h"
#include "field_wrapper.h"






//In pseudocode

// class Reactions;
// class Weights;
// class UpdateRules;
/* 


class Weight_Wrapper {
private:
    int n;
    int N0;
    int N1;
    double **calculated_weights;
public:
    void Calculate_Weights();
};

class Update_Rules {
private:
    int n;
    int N0;
    int N1;
    double **updated;
public:
    void Calculate_Updates();
};


//template <typename Q, typename Q2>
struct CH {

Reaction_Wrapper a;
Weight_Wrapper b;
Update_Rules c;

// we should view these as a container, which take as the input the fields (with some parameters
// and produce the real space output of the things to be transformed

void CalculateWeights();
void CalculateChemistries();
void DoTheFouriers();
void DoTheUpdateRules();
void DoTheInverseFouriers();

//then we go back to the beginning


};
*/

template <class T>
struct CH { //standard CH, everything real, cosine boundary conditions

CH_builder myp;


T **fields; //tehse are my fields

Field_Wrapper<T, T> chems;
Field_Wrapper<T, T> weigs;
Field_Wrapper<T, T> transformed1;
Field_Wrapper<T, T> transformed2;
Field_Wrapper<T, T> transformed3;

Rule_Wrapper<T,T,T,T> rules;
Field_Wrapper<T,T> reverse_transform;

CH(const CH_builder &p);

~CH();

void set_field(const matrix<T>&,int);

void set_field(T**);

void set_chems(const Field_Wrapper<T, T> &a) { chems = a; }

void set_weights(const Field_Wrapper<T, T> &a) { weigs = a; }

void set_rules(const Rule_Wrapper<T,T,T,T> &a) { rules = a; }

void print_all_results(string s1);

virtual void Update();

};

struct CHN : public CH<complex<double> > { //define new cahn hilliard method for those things which have order parameter dependent mobility 
    //double **fields;

    CHN(const CH_builder &p) : CH(p) {}
    Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > newrules;

    void set_new_rules(const Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > &a) { newrules = a; }

    void Update();
};

template <class Q>
struct CHWithNoise : public CH<complex<double>>
{ // define new cahn hilliard method for those things which have order parameter dependent mobility
    // double **fields;
    Q &func;
    GenNoise<complex<double> > mynoise;
    double *str;

    CHWithNoise(const CH_builder &p, Q &funcc2, double *strs) : CH(p), mynoise(GenNoise<complex<double> >(p)), func(funcc2) {

        str = new double [p.number_of_fields];
        for(int i = 0  ; i < p.number_of_fields ; i++) {
            str[i] = strs[i];
        }

    }

    //void set_new_rules(const Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>> &a) { newrules = a; }

    void Update();
};

#include "cahnhilliard.cpp"

#endif