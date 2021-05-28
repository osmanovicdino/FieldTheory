#ifndef CHBUILDER_H
#define CHBUILDER_H

struct CH_builder
{

    int number_of_fields;
    int N1;
    int N2;

    int get_total() const { return N1*N2; }
};


#endif /* CHBUILDER_H */
