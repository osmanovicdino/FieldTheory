#ifndef CHEMISTRY_H
#define CHEMISTRY_H


// class Chemistry
// {
// public:
//     virtual void operator()(double **a, double **fields, int j, int n) = 0;
//     virtual Chemistry *clone() const = 0;
//     virtual void print() {cout << "the base class\n"; }
// };

class NoWeight : public Weight {
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
                      a[j][i] = 0.0;
    }
    NoWeight *clone() const
    {
        return new NoWeight(*this);
    }
    void print() { cout << "no chemistry\n"; }
};

class InvasionAntiInvasionLinear : public Weight {
    // A simple invader minus catalytic releaser chemistry 
    //WARNING: WE DO NO  CHECK FOR OUT OF BOUND ERRORS DUE TO PERFORMANCE DECREASES
    //IT IS YOUR RESPONSIBILITY TO MAKE SURE THAT EVERYTHING IS HALAL
private:
    double rate1;
    double rate2;
public:
    InvasionAntiInvasionLinear(double r1, double r2) : rate1(r1), rate2(r2) {}

    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful 
        int end = p.get_total();

        for(int i = 0 ; i < end ; i++)    {
            a[j][i] = -rate1*fields[2][i]*fields[0][i]+rate2*fields[1][i]*fields[3][i];
        }
    }

    void print() { cout << "Invasion Anti Invasion Linear\n rate1: " << rate1 << " rate2: " << rate2 << endl; }

    InvasionAntiInvasionLinear *clone() const
    {
        return new InvasionAntiInvasionLinear(*this);
    }
};

class InvasionSelfRelease : public Weight
{
    // A simple invader minus catalytic releaser chemistry
    //WARNING: WE DO NO  CHECK FOR OUT OF BOUND ERRORS DUE TO PERFORMANCE DECREASES
    //IT IS YOUR RESPONSIBILITY TO MAKE SURE THAT EVERYTHING IS HALAL
private:
    double rate1;
    double rate2;
    int val1;
    int val2;
    int val3;

public:
    InvasionSelfRelease(double r1, double r2, int val11, int val22, int val33) : rate1(r1), rate2(r2), val1(val11), val2(val22), val3(val33) {}

    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();

        for (int i = 0; i < end; i++)
        {
            a[j][i] = -rate1 * fields[val1][i] * fields[val3][i] + rate2 * fields[val1][i] * fields[val2][i];
        }
    }

    void print() { cout << "Invasion Anti Invasion Linear\n rate1: " << rate1 << " rate2: " << rate2 << endl; }

    InvasionSelfRelease *clone() const
    {
        return new InvasionSelfRelease(*this);
    }
};

//is there anything unique to Chemistry in reaction wrapper?




#include "chemistry.cpp"

#endif /* CHEMISTRY_H */


