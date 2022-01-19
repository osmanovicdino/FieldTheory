#ifndef UPDATERULES_H
#define UPDATERULES_H

template<class result, class input1, class input2, class input3>
struct updateRules {
CH_builder par;
    //void what_am_I() { cout << "I am the base class for Update Rules" << endl; }
// updateRules() {par.N1=1; par.N2=1; par.number_of_fields=1;}
updateRules(const CH_builder &para) : par(para) {}
void setpar(const CH_builder &para);
virtual void operator()(result **, input1**,input2**,input3**,int,const CH_builder &) = 0;
//virtual void trylower(result **, input1 **, input2 **, input3 **, int, const CH_builder &, double newdt) = 0;
virtual updateRules *clone() const = 0;
virtual void print() = 0;

};

template <class result, class input1, class input2, class input3>
void updateRules<result,input1,input2,input3>::setpar(const CH_builder &para)
{
    par = para;
}

#endif /* UPDATERULES_H */
