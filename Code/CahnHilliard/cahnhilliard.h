
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

void set_field(T*,int);

void set_chems(const Field_Wrapper<T, T> &a) { chems = a; }

void set_weights(const Field_Wrapper<T, T> &a) { weigs = a; }

void set_rules(const Rule_Wrapper<T,T,T,T> &a) { rules = a; }

void print_all_results(string s1);
void print_some_results(string s1,vector1<bool>&ps);

void check_field(bool &a) {
    if(fields[0][0]!=fields[0][0]) a = false;
}

void high_k_correct(int cut_offf);

virtual void Update();

};

template <>
void CH<complex<double>>::set_field(complex<double> **orig)
{
if (myp.dimension == 2)
{
    for (int k = 0; k < myp.number_of_fields; k++)
    {
        for (int i = 0; i < myp.N1; i++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                fields[k][i * myp.N2 + j] = {orig[k][i * myp.N2 + j].real(), 0.};
                // fields[k][i * myp.N2 + j] = {orig[k][i * myp.N2 + j],0.};
            }
        }
    }
}
else
{
    for (int k = 0; k < myp.number_of_fields; k++)
    {
        for (int i = 0; i < myp.N1; i++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                for (int lk = 0; lk < myp.N3; lk++)
                {
                    fields[k][i * myp.N2 * myp.N3 + j * myp.N3 + lk] = {orig[k][i * myp.N2 * myp.N3 + j * myp.N3 + lk].real(), 0.};
                    // fields[k][i * myp.N2 + j] = orig[k][i * myp.N2 + j].real();
                }
            }
        }
    }
}
}

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
    vector1<double> str;

    CHWithNoise(const CH_builder &p, Q &funcc2, vector1<double> &strs) : CH(p), str(strs), mynoise(GenNoise<complex<double> >(p)), func(funcc2) {



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



typedef complex<double> cd;
template<typename T>
struct CHC : public CH<T>
{
    matrix<double> epsilon_couplings;

    vector1<bool> phase_separators;

    matrix<double> epsilon_couplingsSQR;
    vector1<double> diffusion;
    vector1<double> epsilon;
    vector1<double> c0;
    vector1<double> c1;
    vector1<double> cons1,cons2,cons3,cons4;
    vector1<double> cons3s;
    double temp1;

    double alpha;
    double dt;

    Field_Wrapper<T, T> oldfieldFT;
    Field_Wrapper<T, T> oldfieldNLW;
    //Field_Wrapper<complex<double>, complex<double>> oldweightFT;

    Field_Wrapper<T, T> InitWeight;
    

    vector<matrix<double> > inverses; //inverse of each update matrix for each value of k1,k2
    vector<matrix<double> > baremat; //each matrix for each value of k1,k2

    

    CHC(const CH_builder &p);

    void set_phase_separators(vector1<bool> &ps) {phase_separators = ps;}
    void set_interaction(double val, int i, int j);
    void set_diffusion(double diff,int i) {diffusion[i] = diff; cout << "diffusion set to: " << diff << endl;}
    void set_epsilon(double epss, int i) {epsilon[i] = epss;
        cout << "epsilon set to: " << epss << endl;
    }
    void set_c0_c1(double c00, double c11, int i, double nu1 = 1.0) {
        if(!phase_separators[i]) error("trying to set phase separation parameters for non-phase separating model");
        c0[i] =  c00; c1[i] = c11;
        double nu = nu1;
        cons1[i] = 4 * nu;
        cons2[i] = (-6 * c0[i] * nu - 6 * c1[i] * nu);
        cons3[i] = (2 * c0[i] * c0[i] * nu + 8 * c0[i] * c1[i] * nu + 2 * c1[i] * c1[i] * nu);
        cons4[i] = -2 * c0[i] * c0[i] * c1[i] * nu - 2 * c0[i] * c1[i] * c1[i] * nu;
        cons3s[i] = cons3[i] - 1;
        cout << "c0 set to: " << c00 << endl;
        cout << "c1 set to: " << c11 << endl;
    }
    void set_temp1(double temp11)
    {
        temp1 = temp11;
        cout << "temp set to: " << temp1 << endl;
    }

    void set_alpha(double alphaa) {alpha = alphaa;
        cout << "alpha set to: " << alpha << endl;
    }
    void set_dt(double dtt) { dt= dtt;
        cout << "dt set to: " << dt << endl;
    }

    void setup_matrices();

    double getk(int k1, int k2, int k3, int N1, int N2, int N3, int &rel);
    double getk(int k1, int k2, int N1, int N2, int &rel);
    double getk(int k1, int N1, int &rel);

    matrix<double> create_D_mat_split(double k1)
    {
        int field_no = (this->myp).number_of_fields;
        matrix<double> dmat(field_no, field_no);

        for (int i = 0; i < field_no; i++)
        {
            dmat(i, i) += -diffusion[i] * temp1 * (SQR(k1) );
        }

        for (int i = 0; i < field_no; i++)
        {
            for (int j = 0; j < field_no; j++)
            {
                dmat(i, j) += -diffusion[i] * temp1 * (SQR(k1) ) * epsilon_couplings(i, j);
            }
        }
        for (int i = 0; i < field_no; i++)
        {
            if (phase_separators[i])
                dmat(i, i) += -diffusion[i] * cons3s[i] * temp1 * (SQR(k1)) - diffusion[i] * SQR(epsilon[i]) * SQR(temp1) * SQR(SQR(k1));
        }
        return dmat;
    }

    matrix<double> create_D_mat_split(double k1, double k2) {
        int field_no = (this->myp).number_of_fields;
        matrix<double> dmat(field_no,field_no);

        for(int i = 0 ; i < field_no ; i++) {
            dmat(i, i) += - diffusion[i] *temp1 *(SQR(k1) + SQR(k2));
        }

        for(int i = 0 ; i < field_no ; i++) {
            for(int j  = 0 ; j < field_no ; j++) {
                dmat(i, j) += -diffusion[i] * temp1 * (SQR(k1) + SQR(k2)) *epsilon_couplings(i, j);
            }
        }
        for(int i = 0  ; i < field_no ; i++) {
        if(phase_separators[i])
        dmat(i, i) += -diffusion[i] * cons3s[i] * temp1 * (SQR(k1) + SQR(k2)) - diffusion[i] * SQR(epsilon[i]) * SQR(temp1) * SQR(SQR(k1) + SQR(k2));
        }
        return dmat;
    }

    matrix<double> create_D_mat_split(double k1, double k2, double k3)
    {
        int field_no = (this->myp).number_of_fields;
        matrix<double> dmat(field_no, field_no);

        for (int i = 0; i < field_no; i++)
        {
            dmat(i, i) += -diffusion[i] * temp1 * (SQR(k1) + SQR(k2)+SQR (k3));
        }

        for (int i = 0; i < field_no; i++)
        {
            for (int j = 0; j < field_no; j++)
            {
                dmat(i, j) += -diffusion[i] * temp1 * (SQR(k1) + SQR(k2) + SQR(k3)) * epsilon_couplings(i, j);
            }
        }
        for (int i = 0; i < field_no; i++)
        {
            if (phase_separators[i])
                dmat(i, i) += -diffusion[i] * cons3s[i] * temp1 * (SQR(k1) + SQR(k2)+SQR(k3)) - diffusion[i] * SQR(epsilon[i]) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)+SQR(k3));
        }
        return dmat;
    }

    void calculate_non_linear_weight(T **);

    void calculate_non_linear_weightSQR(T **);

    void calculate_initial_weight(int);
    void calculate_initial_weightSQR();

    void Update();

    void UpdateTD(matrix<double> &fieldtoset, int i);

    template <class Q>
    void UpdateNoise(Q &func,GenNoise<T> &,vector1< double>&);

    void UpdateSQR();


};



struct CHD : public CH<double>
{
    matrix<double> epsilon_couplings;

    matrix<double> epsilon_couplingsSQR;

    double chi_12;

    double chi_13;
    double chi_23;

    // double diffusion1;
    // double diffusion2;

    matrix<double> diffusion_matrix;

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

    void set_interaction(double val, int i, int j) {
       if((i==0 && j==1)||(i==1&&j==0) ) {
        error("trying to set interaction in the double phase separating system between components 0 and 1, however, these are the components that are defined to be phase separating, and thus already have an interaction");
       }
       else{
        epsilon_couplings(i, j) = val;
        epsilon_couplings(j, i) = val;
       }
    }

    vector<matrix<double>> inverses; // inverse of each update matrix for each value of k1,k2
    vector<matrix<double>> baremat;

    CHD(const CH_builder &p);
    void setup_matrices();

    void set_interaction_and_diffusion(double x12, double x13, double x23, matrix<double> D1);

    matrix<double> create_D_mat_split(double k1, double k2) {
        int field_no = myp.number_of_fields;
        matrix<double> dmat(field_no, field_no);

         dmat(0, 0) = 1 + dt * temp1 * (SQR(k1) + SQR(k2)) * ml1 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl1;
         dmat(0, 1) = dt * temp1 * (SQR(k1) + SQR(k2)) * ml2 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl2;
         dmat(1, 0) = dt * temp1 * (SQR(k1) + SQR(k2)) * ml3 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl3;
         dmat(1, 1) = 1 + dt * temp1 * (SQR(k1) + SQR(k2)) * ml4 - dt * SQR(epsilon) * SQR(temp1) * SQR(SQR(k1) + SQR(k2)) * sl4;


        //  for (int i = 0; i < field_no; i++)
        //  {
        //      for (int j = 0; j < field_no; j++)
        //      {
        //          dmat(i, j) += dt * diffusion_matrix(i, i) * temp1 * (SQR(k1) + SQR(k2)) * epsilon_couplings(i, j);
        //      }
        // }



        // for(int i = 2 ; i < field_no ; i++) {
        //     dmat(i, i) += 1. + dt* diffusion_matrix(i, i) * temp1 * (SQR(k1) + SQR(k2));
        // }




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
                // if (i1 <= myp.N1 / 2)
                // {
                    k1 = i1;
                // }
                // else
                // {
                //     k1 = (i1 - myp.N1);
                // }
                // if (j <= myp.N2 / 2)
                // {
                    k2 = j;
                // }
                // else
                // {
                //     k2 = (j - myp.N2);
                // }

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

    void calculate_non_linear_weight(double **);

    void Update();

    void Update_With_Chem(); //turn on the chemical reactions
};

#include "cahnhilliard.cpp"
#include "cahnhilliardcombo.cpp"
#include "cahnhilliarddouble.cpp"

#endif