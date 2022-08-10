
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

void check_field(bool &a) {
    if(fields[0][0]!=fields[0][0]) a = false;
}

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

struct CHSubFrac : public CH<complex<double>>
{
    Field_Wrapper<complex<double>, complex<double>> initial_cond;
    // Field_Wrapper< complex<double>, complex<double> > fractional_weight;
    Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>> fractional_weight;

    void set_init_cond(const Field_Wrapper<complex<double>, complex<double>> &a) { initial_cond = a; }
    void set_frac_weight(const Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>> &a) { fractional_weight = a; }

    CHSubFrac(const CH_builder &p) : CH(p)
    {
    }

    void Update();
};

struct CHFrac : public CH<complex<double> > 
{
    Field_Wrapper < complex<double>, complex<double> > initial_cond;
    Field_Wrapper < complex<double>, complex<double> > old_fields; //need to set old fields as well

    // Field_Wrapper< complex<double>, complex<double> > fractional_weight;
    Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > fractional_weight;

    void set_init_cond(const Field_Wrapper<complex<double>, complex<double> > &a) {initial_cond = a; }
    void set_old_fields(const Field_Wrapper<complex<double>, complex<double> > &a) { old_fields = a; }
    void set_frac_weight(const Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > &a) { fractional_weight = a; }

    CHFrac(const CH_builder &p) : CH(p) 
    {

    }

    void Update();
};




struct CHFracDt : public CH<complex <double> > {

    Field_Wrapper< complex<double>, complex<double> > old_fieldsft;
    Field_Wrapper<complex<double>, complex<double> > old_weightsft;

    Rule_Wrapper< complex<double>,complex<double>,complex<double>,complex<double> > ri;
    Rule_Wrapper< complex<double>,complex<double>,complex<double>,complex<double> > rA;
    Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>> rB;
    Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>> rC;

    CHFracDt(const CH_builder &p) : CH(p),
                                                         old_fieldsft(Field_Wrapper<complex<double>, complex<double>>(p)),
                                                         old_weightsft(Field_Wrapper<complex<double>, complex<double>>(p)),
                                                         ri(Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>>(p)),
                                                         rA(Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>>(p)),
                                                         rB(Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>>(p)),
                                                         rC(Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double>>(p))
    {
        for (int i = 0; i < myp.number_of_fields; i++)
        {
            IdentityWeight<complex<double> > fw;
            old_fieldsft.add_method(fw, i);
            old_weightsft.add_method(fw, i);
            
        }
    }

    void AddBundleMethod( FracRuleBundle &CDD, int i) {
        ri.add_method(*(CDD.I), i);
        rA.add_method(*(CDD.A), i);
        rB.add_method(*(CDD.B), i);
        rC.add_method(*(CDD.C), i);
    }

    void setupInitial() {

        
        
        transformed1.Calculate_Results(fields);
        

        //do the cut off
        int cut_off_k = 100;

        for (int i = 0; i < myp.number_of_fields; i++)
        {
            for (int i1 = 0; i1 < myp.N1; i1++)
            {
                for (int j = 0; j < myp.N2; j++)
                {
                    double k1, k2;
                    if (i1 <= myp.N1 / 2)
                    {
                        k1 = i1;
                    }
                    else
                    {
                        k1 = (i1 - myp.N1);
                    }
                    if (j <= myp.N2 / 2)
                    {
                        k2 = j;
                    }
                    else
                    {
                        k2 = (j - myp.N2);
                    }

                    //double tempor = SQR(k1) + SQR(k2);
                    if(SQR(k1)+SQR(k2)>SQR(cut_off_k ) ) {
                        transformed1.calculated_reactions[i][i1 * myp.N2 + j] = 0.0;//cut off
                    }
                    // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
                }
            }
        }
        //cout << "cut off done" << endl;

        reverse_transform.Calculate_Results(transformed1.calculated_reactions);



        set_field(reverse_transform.calculated_reactions);
        weigs.Calculate_Results(fields);

        transformed2.Calculate_Results(weigs.calculated_reactions);

        old_fieldsft.Calculate_Results(transformed1.calculated_reactions);

        old_weightsft.Calculate_Results(transformed2.calculated_reactions);

        ri.Calculate_Results(transformed1.calculated_reactions,transformed2.calculated_reactions,fields);


        
        for (int i = 0; i < myp.number_of_fields; i++)
        {
            for(int j = 0  ; j < myp.get_total() ; j++) {
                ri.calculated_reactions[i][j] = transformed1.calculated_reactions[i][j] - ri.calculated_reactions[i][j];
            }
        }



        

        reverse_transform.Calculate_Results(old_fieldsft.calculated_reactions);


        transformed1.Calculate_Results(reverse_transform.calculated_reactions);
        old_fieldsft.Calculate_Results(transformed1.calculated_reactions);

        // reverse_transform.GetMaximas();
        // reverse_transform.GetMaximasIndex();
        // reverse_transform.GetMinimas();
        // reverse_transform.GetMinimasIndex();
        // cout << endl;
        // cout << "old field added" << endl;
        // pausel();

        

        weigs.Calculate_Results(reverse_transform.calculated_reactions);

        transformed2.Calculate_Results(weigs.calculated_reactions);

        old_weightsft.Calculate_Results(transformed2.calculated_reactions);


    }

    void Update() {
        weigs.Calculate_Results(fields);
        chems.Calculate_Results(fields);
        transformed1.Calculate_Results(fields);
        transformed2.Calculate_Results(weigs.calculated_reactions);
        transformed3.Calculate_Results(chems.calculated_reactions);

        rA.Calculate_Results(transformed1.calculated_reactions,transformed2.calculated_reactions,transformed1.calculated_reactions);


        rB.Calculate_Results(old_fieldsft.calculated_reactions, transformed1.calculated_reactions, old_weightsft.calculated_reactions);
        // string filename7 = "test5";
        // outfunc(rA.calculated_reactions[0], filename7, myp);
        // string filename8 = "test6";
        // outfunc(rB.calculated_reactions[0], filename8, myp);
        // string filename3 = "oldf";
        // outfunc(transformed1.calculated_reactions[0], filename3, myp);
        // string filename4 = "test1";
        // outfunc(rA.calculated_reactions[0], filename4, myp);
        // string filename5 = "test2";
        // outfunc(rB.calculated_reactions[0], filename5, myp);
        // pausel();
        rB += rA;
        rC.Calculate_Results(ri.calculated_reactions, rB.calculated_reactions,transformed3.calculated_reactions);

        // string filename3 = "test1";
        // outfunc(old_fieldsft.calculated_reactions[0], filename3, myp);
        // string filename4 = "test2";
        // outfunc(old_weightsft.calculated_reactions[0], filename4, myp);
        // string filename5 = "test3";
        // outfunc(transformed1.calculated_reactions[0], filename5, myp);
        // string filename6 = "test4";
        // outfunc(transformed2.calculated_reactions[0], filename6, myp);


        // string filename9 = "test7";
        // outfunc(rC.calculated_reactions[0], filename9, myp);
        // pausel();

        old_fieldsft.Calculate_Results(transformed1.calculated_reactions);

        old_weightsft.Calculate_Results(transformed2.calculated_reactions);

        reverse_transform.Calculate_Results(rC.calculated_reactions);
        reverse_transform.GetMaximas();
        reverse_transform.GetMaximasIndex();
        reverse_transform.GetMinimas();
        reverse_transform.GetMinimasIndex();
        cout << endl;
        set_field(reverse_transform.calculated_reactions);

    }


};

struct CHC : public CH<complex<double>>
{
    matrix<double> epsilon_couplings;

    matrix<double> epsilon_couplingsSQR;
    double diffusion;
    double epsilon;
    double c0;
    double c1;
    double cons1,cons2,cons3,cons4;
    double cons3s;
    double temp1;

    double alpha;
    double dt;

    Field_Wrapper<complex<double>, complex<double>> oldfieldFT;
    Field_Wrapper<complex<double>, complex<double>> oldfieldNLW;
    //Field_Wrapper<complex<double>, complex<double>> oldweightFT;

    Field_Wrapper<complex<double>, complex<double>> InitWeight;
    

    vector<matrix<double> > inverses; //inverse of each update matrix for each value of k1,k2
    vector<matrix<double> > baremat; //each matrix for each value of k1,k2

    CHC(const CH_builder &p);

    void set_interaction(double val, int i, int j);
    void set_diffusion(double diff) {diffusion = diff; cout << "diffusion set to: " << diff << endl;}
    void set_epsilon(double epss) {epsilon = epss;
        cout << "epsilon set to: " << epss << endl;
    }
    void set_c0_c1(double c00, double c11) {c0 =  c00; c1 = c11;
        double nu = 1.0;
        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
        cons3s = cons3-1;
        cout << "c0 set to: " << c00 << endl;
        cout << "c1 set to: " << c11 << endl;
    }
    void set_temp1(double temp11) {temp1 = temp11;
        cout << "temp set to: " << temp1 << endl;
    }

    void set_alpha(double alphaa) {alpha = alphaa;
        cout << "alpha set to: " << alpha << endl;
    }
    void set_dt(double dtt) { dt= dtt;
        cout << "dt set to: " << dt << endl;
    }

    void setup_matrices();

    matrix<double> create_D_mat_split(double k1, double k2) {
        int field_no = myp.number_of_fields;
        matrix<double> dmat(field_no,field_no);

        for(int i = 0 ; i < field_no ; i++) {
            dmat(i, i) += - diffusion *temp1 *(SQR(k1) + SQR(k2));
        }

        for(int i = 0 ; i < field_no ; i++) {
            for(int j  = 0 ; j < field_no ; j++) {
                dmat(i, j) += -diffusion * temp1 * (SQR(k1) + SQR(k2)) *epsilon_couplings(i, j);
            }
        }
        dmat(0, 0) += -diffusion * cons3s * temp1 * (SQR(k1) + SQR(k2)) - diffusion * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2));
        return dmat;
    }

    void calculate_non_linear_weight(complex<double> **);

    void calculate_non_linear_weightSQR(complex<double> **);

    void calculate_initial_weight();
    void calculate_initial_weightSQR();

    void Update();

    template <class Q>
    void UpdateNoise(Q &func,GenNoise<complex<double> > &,vector1< double>&);

    void UpdateSQR();
};

struct CHD : public CH<complex<double>>
{

    double chi_12;

    double chi_13;
    double chi_23;

    double diffusion1;
    double diffusion2;

    double x0; // minima of the density for comp1
    double y0; //minima of the density for comp2

    double ml1;
    double ml2;
    double ml3;
    double ml4;

    double sl1;
    double sl2;
    double sl3;
    double sl4;

    double dt;
    double temp1;
    double epsilon;

    void set_temp1(double temp11)
    {
        temp1 = temp11;
        cout << "temp set to: " << temp1 << endl;
    }

    void set_epsilon(double epsilonn)
    {
        epsilon = epsilonn;
        cout << "eps set to: " << epsilon << endl;
    }
    void set_dt(double dtt)
    {
        dt = dtt;
        cout << "dt set to: " << dt << endl;
    }

    vector<matrix<double>> inverses; // inverse of each update matrix for each value of k1,k2
    vector<matrix<double>> baremat;

    CHD(const CH_builder &p);
    void setup_matrices();

    void set_interaction_and_diffusion(double x12, double x13, double x23, double D1, double D2);


    matrix<double> create_D_mat_split(double k1, double k2) {
        int field_no = myp.number_of_fields;
        matrix<double> dmat(field_no, field_no);

        dmat(0, 0) = 1 + dt * temp1 * (SQR(k1) + SQR(k2)) * ml1 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl1;
        dmat(0, 1) = dt * temp1 * (SQR(k1) + SQR(k2)) * ml2 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl2;
        dmat(1, 0) = dt * temp1 * (SQR(k1) + SQR(k2)) * ml3 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl3;
        dmat(1, 1) = 1 + dt * temp1 * (SQR(k1) + SQR(k2)) * ml4 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl4;

        return dmat;
    }

    void setupInitial(vector1<int> cut_off_k) {
    transformed1.Calculate_Results(fields);

    for (int i = 0; i < myp.number_of_fields; i++)
    {
        for (int i1 = 0; i1 < myp.N1; i1++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                double k1, k2;
                if (i1 <= myp.N1 / 2)
                {
                    k1 = i1;
                }
                else
                {
                    k1 = (i1 - myp.N1);
                }
                if (j <= myp.N2 / 2)
                {
                    k2 = j;
                }
                else
                {
                    k2 = (j - myp.N2);
                }

                // double tempor = SQR(k1) + SQR(k2);
                if (SQR(k1) + SQR(k2) > SQR(cut_off_k[i]))
                {
                    transformed1.calculated_reactions[i][i1 * myp.N2 + j] = 0.0; // cut off
                }
                // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
            }
        }
    }
    // cout << "cut off done" << endl;

    reverse_transform.Calculate_Results(transformed1.calculated_reactions);

    set_field(reverse_transform.calculated_reactions);
    }

    void calculate_non_linear_weight(complex<double> **);

    void Update();


};

#include "cahnhilliard.cpp"
#include "cahnhilliardcombo.cpp"
#include "cahnhilliarddouble.cpp"

#endif