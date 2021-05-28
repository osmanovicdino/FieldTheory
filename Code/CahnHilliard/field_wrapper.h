#ifndef FIELD_WRAPPER_H
#define FIELD_WRAPPER_H

#include "chbuilder.h"
#include "weights.h"
#include "chemistry.h"
#include "updateRules.h"

struct Field_Wrapper
{

    CH_builder params;
    double **calculated_reactions;
    Weight **chems;

    Field_Wrapper();
    Field_Wrapper(const CH_builder &p);
    Field_Wrapper(const Field_Wrapper &a);
    Field_Wrapper &operator=(const Field_Wrapper &); //assignment operator

    ~Field_Wrapper();

    void add_method(Weight&, int); //adds a method of manipulating field i

    void Calculate_Results(double **fields);
};


struct Rule_Wrapper {

    CH_builder params;
    double **calculated_reactions;
    updateRules **chems;

    Rule_Wrapper();
    Rule_Wrapper(const CH_builder &p);
    Rule_Wrapper(const Rule_Wrapper &a);
    Rule_Wrapper &operator=(const Rule_Wrapper &); //assignment operator

    ~Rule_Wrapper();

    void add_method(updateRules &, int); //adds a method of manipulating field i

    void Calculate_Results(double **fields, double **weis, double **reacs);
};

#include "field_wrapper.cpp"
#include "rule_wrapper.cpp"

#endif /* FIELD_WRAPPER_H */
