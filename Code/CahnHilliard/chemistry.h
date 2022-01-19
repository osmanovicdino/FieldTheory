#ifndef CHEMISTRY_H
#define CHEMISTRY_H


// class Chemistry
// {
// public:
//     virtual void operator()(double **a, double **fields, int j, int n) = 0;
//     virtual Chemistry *clone() const = 0;
//     virtual void print() {cout << "the base class\n"; }
// };

template<class T, class Q>
class NoWeight : public Weight<T, Q> {
public:
    void operator()(T **a, Q **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
                      a[j][i] = T(0.0); //return nothing
    }
    NoWeight<T,Q> *clone() const
    {
        return new NoWeight<T,Q>(*this);
    }
    void print() { cout << "no chemistry\n"; }
};

template <class T>
class GenericChemistry: public Weight<T,T>
{
private:
double rate1;
vector1<int> pows;
int n;

public: 
GenericChemistry(double rate11, vector1<int> powss) : pows(powss) { n = pows.getsize(); rate1 = rate11; }
void operator()(T **a, T **fields, int j, const CH_builder &p)
{
    // as there is crosstalk here it is required for you to be careful
    int end = p.get_total();

    for (int i = 0; i < end; i++)
    {
        double tot = 1.0;
        for(int k = 0  ; k < n ; k++) {
            tot *= Power(fields[k][i].real(),pows[k]);
        }
        a[j][i] = rate1 * tot;
    }
}
GenericChemistry *clone() const
{
    return new GenericChemistry(*this);
}
};

template <class T> //allow for complex data
class InvasionDecay : public Weight<T, T>
{
private:
    double rate1;
    int which1;

public:
    InvasionDecay(double r1, int which11) : rate1(r1), which1(which11) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();

        for (int i = 0; i < end; i++)
        {
            a[j][i] = -rate1 * fields[which1][i];
        }
    }

    void print() { cout << "Invasion Linear\n rate1: " << rate1 << endl; }

    InvasionDecay *clone() const
    {
        return new InvasionDecay(*this);
    }
};

template <class T> //allow for complex data
class InvasionLinear : public Weight<T,T>
{
private:
    double rate1;
    int which1;
    int which2;

public:
    InvasionLinear(double r1, int which11, int which22) : rate1(r1), which1(which11), which2(which22) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();

        for (int i = 0; i < end; i++)
        {
            a[j][i] = -rate1 * fields[which1][i] * fields[which2][i];
        }
    }

    void print() { cout << "Invasion Linear\n rate1: " << rate1 << endl; }

    InvasionLinear *clone() const
    {
        return new InvasionLinear(*this);
    }
};

template <class T>
class InvasionSquared : public Weight<T, T>
{
private:
    double rate1;
    int which1;
    int which2;

public:
    InvasionSquared(double r1, int which11, int which22) : rate1(r1), which1(which11), which2(which22) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();

        for (int i = 0; i < end; i++)
        {
            a[j][i] = -rate1 * SQR(fields[which1][i]) * fields[which2][i];
        }
    }

    void print() { cout << "Invasion Squared\n rate1: " << rate1 << endl; }

    InvasionSquared *clone() const
    {
        return new InvasionSquared(*this);
    }
};

template <class T>
class InvasionLinearReversibleA : public Weight<T, T>
{
private: // this is just a reversible reaction up to a larger one, containing the loss and the gain due to the formation and loss
    double rate1;
    double rate2;
    int which1;
    int which2;
    int which3;

public:
    InvasionLinearReversibleA(double r1, double r2, int which11, int which22, int which33) : rate1(r1), rate2(r2), which1(which11), which2(which22), which3(which33) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();

        // cout << j << endl;
        // cout << which1 << " " << which2 << " " << which3 << endl;

        // cout << end << endl;

        for (int i = 0; i < end; i++)
        {
           
            a[j][i] = -rate1 * fields[which1][i] * fields[which2][i] + rate2 * fields[which3][i];
        }
        
    }

    void print() { cout << "InvasionLinearReversibleA\n rate1: " << rate1 << endl; }

    InvasionLinearReversibleA *clone() const
    {
        return new InvasionLinearReversibleA(*this);
    }
};

template <class T>
class InvasionLinearReversibleB : public Weight<T,T> //above and below
{
private:
    double rate1;
    double rate2;
    double rate3;
    double rate4;
    int which_am_i;
    int which_below;
    int which_above;
    int which_invader;

public:
    InvasionLinearReversibleB(double r1, double r2, double r3, double r4, int whicha, int which11, int which22, int which33) : rate1(r1), rate2(r2), rate3(r3), which_am_i(whicha) , which_below(which11), which_above(which22), which_invader(which33) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        // as there is crosstalk here it is required for you to be careful
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            
            a[j][i] = -rate1 * fields[which_am_i][i] * fields[which_invader][i] + rate2 * fields[which_below][i] * fields[which_invader][i] - rate3 * fields[which_am_i][i] + rate4 * fields[which_above][i];
        }
    }

    void print() { cout << "InvasioLinearReversibleB\n rate1: " << rate1 << endl; }

    InvasionLinearReversibleB *clone() const
    {
        return new InvasionLinearReversibleB(*this);
    }
};


template <class T>
class MultipleReactions : public Weight<T,T>
{

private:
Weight<T,T>** multichem;
int no_chem;

public:
    MultipleReactions(int noo) : no_chem(noo) {
        multichem = new Weight<T,T> *[no_chem];

        for(int i = 0  ; i < no_chem ; i++) {
        NoWeight<T,T> *chem1 = new NoWeight<T,T>;
        multichem[0] = chem1;
        }
    }

    MultipleReactions(const MultipleReactions &a) {
       no_chem = a.no_chem;
       multichem = new Weight<T,T> *[no_chem];
       for (int i = 0; i < no_chem; i++)
       {
           Weight<T,T> *chem1 = (a.multichem[i])->clone();
           multichem[i] = chem1;
       }
    }

    MultipleReactions& operator=(const MultipleReactions &a) {
        for (int i = 0; i < no_chem; i++)
        {
            delete multichem[i];
            }
        delete multichem;

        no_chem = a.no_chem;
        multichem = new Weight<T,T> *[no_chem];
        for (int i = 0; i < no_chem; i++)
        {
            Weight<T,T> *chem1 = (a.multichem[i])->clone();
            multichem[i] = chem1;
        }
        return *this;
    }

    ~MultipleReactions() {
        for (int i = 0; i < no_chem; i++)
        {
            delete multichem[i];
        }
        delete multichem;
    }


    void operator()(T **a, T **fields, int j, const CH_builder &p) {
        T **b;
        int Ng = p.get_total();
        b = new T *[p.number_of_fields];

 

        for(int i = 0 ; i < no_chem ; i++)
        b[i] = (T *)fftw_malloc(Ng * sizeof(T));

        for(int i = 0 ; i < no_chem ; i++) {


            (multichem[i])->operator()(b,fields,i,p);
            //cout << "chem done" << endl;
            int end = p.get_total();

            for(int k = 0 ; k < end; k++)
            a[j][k] += b[i][k];
        }

        for (int i = 0; i < no_chem; i++)
        {
            fftw_free(b[i]);
        }
        delete b;
    }

    void add_chemical_reaction(Weight<T,T> &a, int j)
    {
        multichem[j] = a.clone();
    }

    MultipleReactions *clone() const {
        return new MultipleReactions(*this);
    }
};

template <class T>
class InvasionAntiInvasionLinear : public Weight<T,T>
{
    // A simple invader minus catalytic releaser chemistry 
    //WARNING: WE DO NO  CHECK FOR OUT OF BOUND ERRORS DUE TO PERFORMANCE DECREASES
    //IT IS YOUR RESPONSIBILITY TO MAKE SURE THAT EVERYTHING IS HALAL
private:
    double rate1;
    double rate2;
public:
    InvasionAntiInvasionLinear(double r1, double r2) : rate1(r1), rate2(r2) {}

    void operator()(T **a, T **fields, int j, const CH_builder &p)
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

template  <class T>
class InvasionSelfRelease : public Weight<T, T>
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

    void operator()(T **a, T **fields, int j, const CH_builder &p)
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


