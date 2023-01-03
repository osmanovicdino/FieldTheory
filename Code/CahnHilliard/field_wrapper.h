#ifndef FIELD_WRAPPER_H
#define FIELD_WRAPPER_H

#include "chbuilder.h"
#include "weights.h"
#include "chemistry.h"
#include "updateRules.h"

template <class T, class Q> //maps from data type Q to T
struct Field_Wrapper
{

    CH_builder params;
    T **calculated_reactions;
    Weight<T,Q> **chems;

    Field_Wrapper();
    Field_Wrapper(const CH_builder &p);
    Field_Wrapper(const Field_Wrapper<T,Q> &a);
    Field_Wrapper &operator=(const Field_Wrapper<T,Q> &); //assignment operator

    ~Field_Wrapper();

    void add_method(Weight<T,Q>&, int); //adds a method of manipulating field i

    void Calculate_Results(Q **fields);

    void Add_Noise(GenNoise<T> &a);

    void set_field(const matrix<T> &mymat, int k)
    {
        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                calculated_reactions[k][i * params.N2 + j] = mymat.gpcons(i, j);
            }
        }
    }
    void set_field(T **orig)
    {
        for (int k = 0; k < params.number_of_fields; k++)
        {
            for (int i = 0; i < params.N1; i++)
            {
                for (int j = 0; j < params.N2; j++)
                {
                    calculated_reactions[k][i * params.N2 + j] = orig[k][i * params.N2 + j];
                }
            }
        }
    }
    void means(vector1<T> &v);
    void means();
    void Check_fields();
    void rescale();
    void GetMaximas();
    void GetMaximas(vector1<T> &vv);
    void GetMinimas();
    void GetMinimas(vector1<T> &vv);
    void GetMaximasIndex();
    void GetMinimasIndex();
    
};

template <class T, class T1, class T2, class T3> //output class T, input classes T1,T2,T3
struct Rule_Wrapper
{

    CH_builder params;
    T **calculated_reactions;
    updateRules<T,T1,T2,T3> **chems;

    Rule_Wrapper();
    Rule_Wrapper(const CH_builder &p);
    Rule_Wrapper(const Rule_Wrapper &a);
    Rule_Wrapper &operator=(const Rule_Wrapper &); //assignment operator
    Rule_Wrapper &operator+=(const Rule_Wrapper&);

    ~Rule_Wrapper();

    void add_method(updateRules<T, T1, T2, T3> &, int); //adds a method of manipulating field i

    void Calculate_Results(T1 **fields, T2 **weis, T3 **reacs);

    void Check_fields();

    void GetMaximas();
    void GetMinimas();
    void GetMaximasIndex();
    void GetMinimasIndex();
    void GetTotal();
    
    void Check_negatives(bool &he);
};

#include "updateRules.cpp"
#include "updateRulesFrac.cpp"
#include "field_wrapper.cpp"
#include "rule_wrapper.cpp"

#endif /* FIELD_WRAPPER_H */
