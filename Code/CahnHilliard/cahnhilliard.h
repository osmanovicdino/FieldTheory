
#ifndef CAHNHILLIARD_H
#define CAHNHILLIARD_H

#include "chbuilder.h"
#include "utilityOutputs.h"
#include "fftw3.h"
#include "weights.h"
#include "chemistry.h"
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

struct CH {

CH_builder myp;


double **fields; //tehse are my fields

Field_Wrapper chems;
Field_Wrapper weigs;
Field_Wrapper transformed1;
Field_Wrapper transformed2;
Field_Wrapper transformed3;


Rule_Wrapper rules;
Field_Wrapper reverse_transform;

CH(const CH_builder &p);

~CH();

void set_field(const matrix<double>&,int);

void set_field(double**);

void set_chems(const Field_Wrapper &a) { chems = a;}

void set_weights(const Field_Wrapper &a) {weigs = a;}

void set_rules(const Rule_Wrapper &a) {rules = a;}

void print_all_results(string s1);

void Update();






};

#include "cahnhilliard.cpp"

#endif