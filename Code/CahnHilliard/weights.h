#ifndef WEIGHTS_H
#define WEIGHTS_H


template <class T, class Q> //goes to T, from Q
class Weight {
    public:
        virtual void operator()(T **a, Q **fields, int j, const CH_builder &p) = 0;
        virtual Weight *clone() const = 0;
        virtual void print() { cout << "the base Weight class\n"; }
};

template <class T>
class IdentityWeight : public Weight<T, T>
{ // real to real
private:
    double eps;

public:
    IdentityWeight()  {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = fields[j][i];
        }
    }
    IdentityWeight *clone() const
    {
        return new IdentityWeight(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffusiveWeight : public Weight<T, T> { //real to real
private:
    double eps;
public:
    DiffusiveWeight(double epss) : eps(epss) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p) {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = eps*fields[j][i];
        }
    }
    DiffusiveWeight *clone() const
    {
        return new DiffusiveWeight(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffDiffusiveWeight : public Weight<T,T> {
private:
    int which;
    double chi;
public:
    DiffDiffusiveWeight(double chii, int whichh) : which(whichh), chi(chii) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = chi* fields[which][i];
        }
    }
    DiffDiffusiveWeight *clone() const
    {
        return new DiffDiffusiveWeight(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffDiffusiveWeightGenericN : public Weight<T, T>
{
private:
    vector1<int> whichs;
    vector1<double> chis;

public:
    DiffDiffusiveWeightGenericN(vector1<double> chii, vector1<int> whichh) : whichs(whichh), chis(chii) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            T as = 0.0;
            for (int k = 0; k < whichs.getsize(); k++)
                as += chis[k] * (fields[whichs[k]][i]);
            a[j][i] = as;
        }
    }
    DiffDiffusiveWeightGenericN *clone() const
    {
        return new DiffDiffusiveWeightGenericN(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffDiffusiveWeightGenericNSQR : public Weight<T, T>
{
private:
    vector1<int> whichs;
    vector1<double> chis;

public:
    DiffDiffusiveWeightGenericNSQR(vector1<double> chii, vector1<int> whichh) : whichs(whichh), chis(chii) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            T as = 0.0;
            for (int k = 0; k < whichs.getsize(); k++)
                as += chis[k] * SQR(fields[whichs[k]][i]);
            a[j][i] = as*fields[j][i];
        }
    }
    DiffDiffusiveWeightGenericNSQR *clone() const
    {
        return new DiffDiffusiveWeightGenericNSQR(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffDiffusiveWeightDouble : public Weight<T, T>
{
private:
    int which1;
    int which2;
    double chi1;
    double chi2;

public:
    DiffDiffusiveWeightDouble(double chii1, double chii2, int whichh1, int whichh2) : which1(whichh1), chi1(chii1), which2(whichh2), chi2(chii2) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = chi1 * fields[which1][i] + chi2 * fields[which2][i];
        }
    }
    DiffDiffusiveWeightDouble *clone() const
    {
        return new DiffDiffusiveWeightDouble(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template<class T>
class DiffDiffusiveWeightSQR : public Weight<T,T>
{
private:
    int which;
    double chi;
public:
    DiffDiffusiveWeightSQR(double chii, int whichh) : which(whichh), chi(chii) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = chi * SQR(fields[which][i]) * fields[j][i];
        }
    }
    DiffDiffusiveWeightSQR *clone() const
    {
        return new DiffDiffusiveWeightSQR(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class DiffDiffusiveWeightSQRandLinear : public Weight<double, double>
{
private:
    int which;
    double chi;
    double chi2;

public:
    DiffDiffusiveWeightSQRandLinear(double chii, double chii2, int whichh) : which(whichh), chi(chii),chi2(chii2) {}
    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = chi * SQR(fields[which][i]) * fields[j][i] + chi2 * fields[which][i];
            
        }
    }
    DiffDiffusiveWeightSQRandLinear *clone() const
    {
        return new DiffDiffusiveWeightSQRandLinear(*this);
    }
    void print() { cout << "the Diffusive weight class\n"; }
};

template <class T>
class CahnHilliardWeight : public Weight<T,T>
{
private:
    double c0;
    double c1;
    double nu;
    double cons1;
    double cons2;
    double cons3;
    double cons4;

public:
    CahnHilliardWeight(double c00, double c11, double nuu) {
        c0 = c00;
        c1 = c11;
        nu = nuu;

        cons1 = 4*nu;
        cons2 = (-6*c0*nu-6*c1*nu);
        cons3 = (2*c0*c0*nu+8*c0*c1*nu+2*c1*c1*nu);
        cons4 = -2*c0*c0*c1*nu-2*c0*c1*c1*nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p) {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
            {
                a[j][i] = cons1 * CUB(fields[j][i]) + cons2 * SQR(fields[j][i]) + cons3*(fields[j][i]) + cons4;
            }
    }
    CahnHilliardWeight *clone() const
    {
        return new CahnHilliardWeight(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

template<class T>
class CahnHilliardWithCouplingWeight : public Weight<T,T>
{

private:
    double c0;
    double c1;
    double nu;
    double ei1;
    double ei2;
    double ei3;
    double cons1;
    double cons2;
    double cons3;
    double cons4;
    int which1;
    int which2;
    int which3;

public:
    CahnHilliardWithCouplingWeight(double c00, double c11, double nuu, double ei11, double ei22, double ei33, int which11, int which22, int which33)
    {
        c0 = c00;
        c1 = c11;
        nu = nuu;
        ei1 = ei11;
        ei2 = ei22;
        ei3 = ei33;
        which1 = which11;
        which2 = which22;
        which3 = which33;

        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = ei1 * (fields[which1][i]) + ei2 * (fields[which2][i]) + ei3 * (fields[which3][i]) + cons1 * CUB(fields[0][i]) + cons2 * SQR(fields[0][i]) + cons3 * (fields[0][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeight *clone() const
    {
        return new CahnHilliardWithCouplingWeight(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

template <class T>
class CahnHilliardWithCouplingWeightGenericN : public Weight<T, T>
{

private:
    double c0;
    double c1;
    double nu;
    double cons1;
    double cons2;
    double cons3;
    double cons4;
    vector1<double> coups;
    vector1<int> whichs;


public:
    CahnHilliardWithCouplingWeightGenericN(double c00, double c11, double nuu, vector1<double> epss, vector1<int> whichss) : coups(epss), whichs(whichss)
    {
        c0 = c00;
        c1 = c11;
        nu = nuu;

        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            T as = 0.0;
            for(int k = 0  ; k < whichs.getsize() ; k++)
                as += coups[k] * (fields[whichs[k]][i]);
            // if(i == 0 || i == 1000) {
            //     cout << as << endl;
            //     cout << coups << endl;
            //     cout << whichs << endl;
            //     for(int k = 0  ; k < whichs.getsize() ; k++) 
            //     cout << (fields[whichs[k]][i]) << ",";
            //     cout << endl;
            //     pausel();
            // }
            
            
            a[j][i] = as + cons1 * CUB(fields[j][i]) + cons2 * SQR(fields[j][i]) + cons3 * (fields[j][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeightGenericN *clone() const
    {
        return new CahnHilliardWithCouplingWeightGenericN(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

template <class T>
class CahnHilliardWithCouplingWeightGenericNSQR : public Weight<T, T>
{

private:
    double c0;
    double c1;
    double nu;
    double cons1;
    double cons2;
    double cons3;
    double cons4;
    vector1<double> coups;
    vector1<int> whichs;

public:
    CahnHilliardWithCouplingWeightGenericNSQR(double c00, double c11, double nuu, vector1<double> epss, vector1<int> whichss) : coups(epss), whichs(whichss)
    {
        c0 = c00;
        c1 = c11;
        nu = nuu;

        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            T as = 0.0;
            for (int k = 0; k < whichs.getsize(); k++)
                as += coups[k] * SQR(fields[whichs[k]][i]);
            // if(i == 0 || i == 1000) {
            //     cout << as << endl;
            //     cout << coups << endl;
            //     cout << whichs << endl;
            //     for(int k = 0  ; k < whichs.getsize() ; k++)
            //     cout << (fields[whichs[k]][i]) << ",";
            //     cout << endl;
            //     pausel();
            // }

            a[j][i] = as * fields[j][i] + cons1 * CUB(fields[0][i]) + cons2 * SQR(fields[0][i]) + cons3 * (fields[0][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeightGenericNSQR *clone() const
    {
        return new CahnHilliardWithCouplingWeightGenericNSQR(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

template<class T>
class CahnHilliardWithCouplingWeightSQR : public Weight<T,T>
{

private:
    double c0;
    double c1;
    double nu;
    double ei1;
    double ei2;
    double cons1;
    double cons2;
    double cons3;
    double cons4;
    int which1;
    int which2;

public:
    CahnHilliardWithCouplingWeightSQR(double c00, double c11, double nuu, double ei11, double ei22, int which11, int which22)
    {
        c0 = c00;
        c1 = c11;
        nu = nuu;
        ei1 = ei11;
        ei2 = ei22;
        which1 = which11;
        which2 = which22;

        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = ei1 * SQR(fields[which1][i]) * (fields[j][i]) + ei2 * SQR(fields[which2][i]) * (fields[j][i]) + cons1 * CUB(fields[j][i]) + cons2 * SQR(fields[j][i]) + cons3 * (fields[j][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeightSQR *clone() const
    {
        return new CahnHilliardWithCouplingWeightSQR(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

template <class T>
class CahnHilliardWithCouplingWeightSQRandLinear : public Weight<T,T>
{

private:
    double c0;
    double c1;
    double nu;
    double ei1;
    double ei2;
    double ei2l;
    double cons1;
    double cons2;
    double cons3;
    double cons4;
    int which1;
    int which2;

public:
    CahnHilliardWithCouplingWeightSQRandLinear(double c00, double c11, double nuu, double ei11, double ei22, double ei2ll, int which11, int which22)
    {
        c0 = c00;
        c1 = c11;
        nu = nuu;
        ei1 = ei11;
        ei2 = ei22;
        ei2l =  ei2ll;
        which1 = which11;
        which2 = which22;

        cons1 = 4 * nu;
        cons2 = (-6 * c0 * nu - 6 * c1 * nu);
        cons3 = (2 * c0 * c0 * nu + 8 * c0 * c1 * nu + 2 * c1 * c1 * nu);
        cons4 = -2 * c0 * c0 * c1 * nu - 2 * c0 * c1 * c1 * nu;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = ei1 * SQR(fields[which1][i]) * (fields[j][i]) + ei2 * SQR(fields[which2][i]) * (fields[j][i]) + ei2l * fields[which2][i] + cons1 * CUB(fields[0][i]) + cons2 * SQR(fields[0][i]) + cons3 * (fields[0][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeightSQRandLinear *clone() const
    {
        return new CahnHilliardWithCouplingWeightSQRandLinear(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};



template<class T> //we want a template here in case we start with complex data
class CoupledPhaseSeparatingSystem : public Weight<T, T>
{
private:
    int no_phases;
    matrix<double> chi_between_phases; //from this, get the chis
    vector1<double> spec_volumes;


    public:
    CoupledPhaseSeparatingSystem(int no_phasess,matrix<double> &myint, vector1<double> &spec_volumess) : 
    chi_between_phases(matrix<double>(no_phasess,no_phasess)) ,
    spec_volumes(spec_volumess),
    no_phases(no_phasess)
    {
        

        if(myint.getnrows() != no_phases  || myint.getncols() != no_phases  )
        {
            cout << no_phases << endl;
            cout << myint.getnrows() << endl;
            cout << myint.getncols() << endl;
            error("error in specifying the interaction matrix in PhaseSeparatingSystem");
        }

        if(spec_volumes.getsize() != no_phases) error("spec volumes incorrect");
        for(int i = 0  ; i < no_phases  ; i++) {
            for(int j = 0  ; j < no_phases  ; j++) {

                chi_between_phases(i,j) = 0.5*(2*myint(i,j) - myint(i,i) - myint(j,j));
            }
        }
    }
    int get_number_phase() {
        return no_phases;
    }

    matrix<double> getchi() {
        return chi_between_phases;
    }

    void operator()(T **a, T **fields, int j, const CH_builder &p) {
        int end = p.get_total();
        for(int i = 0  ; i < end ; i++) {
            T totf = 0.0;
            // for(int k = 0 ; k < no_phases ; k++) {
            //     totf+=fields[k][i];
            // }
            T orang = (1. / spec_volumes[j]) + (1. / spec_volumes[j]) * log(fields[j][i]);
            T val2 = 0.0;
            for(int k  = 0 ; k < no_phases ; k++) {
               val2 += fields[k][i] * (chi_between_phases(j,k) );
            }
            T utan = val2;
            // if(j > 0) {
            //     cout << chi_between_phases << endl;
            //     cout << i << endl;
            //     cout << orang << endl;
            //     cout << (1. / spec_volumes[j]) - (1. / spec_volumes[no_phases]) << endl;
            //     cout << (1. / spec_volumes[j]) * log(fields[j][i])  << endl;
            //     cout << (1. / spec_volumes[no_phases]) * log(1-totf) << endl;
            //     cout << totf << endl;
            //     cout << (1./spec_volumes[no_phases]) << endl;
            //     cout << utan << endl;
            //     pausel();
            // }

            a[j][i] = orang + utan;
        }
    }

    void print() {
        cout << "coupled phase separating system" << endl;
    }
    CoupledPhaseSeparatingSystem *clone() const
    {
        return new CoupledPhaseSeparatingSystem(*this);
    }
};

class FourierWeightForward1D : public Weight<complex<double>, complex<double>>
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {

        fftw_plan p2(fftw_plan_dft_1d(p.N1, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_FORWARD, FFTW_ESTIMATE));

        fftw_execute(p2);
        double corr = 1. / (p.N1);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
        // fftw_cleanup();
    }
    FourierWeightForward1D *clone() const
    {
        return new FourierWeightForward1D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightBackward1D : public Weight<complex<double>, complex<double>>
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_dft_1d(p.N1, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p2);

        fftw_destroy_plan(p2);
    }
    FourierWeightBackward1D *clone() const
    {
        return new FourierWeightBackward1D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightForward2D : public Weight<complex<double>, complex<double> >
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {
        
        fftw_plan p2(fftw_plan_dft_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_FORWARD, FFTW_ESTIMATE));

        fftw_execute(p2);
        double corr = 1. / (p.N1);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
        // fftw_cleanup();
    }
    FourierWeightForward2D *clone() const
    {
        return new FourierWeightForward2D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightBackward2D : public Weight<complex<double>, complex<double>>
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_dft_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p2);
        double corr = 1. / (p.N1);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightBackward2D *clone() const
    {
        return new FourierWeightBackward2D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightForward3D : public Weight<complex<double>, complex<double>>
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {
        if(p.dimension!= 3) error("change dimension of CH builder");

        fftw_plan p2(fftw_plan_dft_3d(p.N1, p.N2, p.N3, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_FORWARD, FFTW_ESTIMATE));

        fftw_execute(p2);
        double corr = 1. / (p.N1);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
        // fftw_cleanup();
    }
    FourierWeightForward3D *clone() const
    {
        return new FourierWeightForward3D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightBackward3D : public Weight<complex<double>, complex<double>>
{
private:
public:
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    {
        if (p.dimension != 3)
            error("change dimension of CH builder");
        fftw_plan p2;

        p2 = fftw_plan_dft_3d(p.N1, p.N2, p.N3, reinterpret_cast<fftw_complex *>(fields[j]), reinterpret_cast<fftw_complex *>(a[j]), FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p2);
        double corr = 1. / (p.N2 * p.N3);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightBackward3D *clone() const
    {
        return new FourierWeightBackward3D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class FourierWeightForward2D_r2c : public Weight<complex<double>, double>
{
private:
public:
    void operator()(complex<double> **a, double **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_dft_r2c_2d(p.N1, p.N2, fields[j], reinterpret_cast<fftw_complex *> (a[j]), 1);

        fftw_execute(p2);
        double corr = 1./(p.N1*4.);
        int end = p.get_total();
        for(int i = 0  ; i < end ; i++) {
            a[j][i]*=corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightForward2D_r2c *clone() const
    {
        return new FourierWeightForward2D_r2c(*this);
    }
    void print() {cout << "FourierWeightForwards" << endl;}
};

class FourierWeightBackward2D_c2r : public Weight<double, complex<double>>
{
private:
public:
    void operator()(double **a, complex<double> **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_dft_c2r_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *> (fields[j]), a[j], 1);

        fftw_execute(p2);
        int end = p.get_total();
        double corr = 1. / (p.N1);
        for (int i = 0; i < end ; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightBackward2D_c2r *clone() const
    {
        return new FourierWeightBackward2D_c2r(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};

class CosineWeightForward1D : public Weight<double, double>
{
private:
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_r2r_1d(p.N1, fields[j], a[j], FFTW_REDFT10, 1);

        fftw_execute(p2);
        double corr = 1. / (p.N1 * 2.);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    CosineWeightForward1D *clone() const
    {
        return new CosineWeightForward1D(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class CosineWeightBackward1D : public Weight<double, double>
{
private:
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_r2r_1d(p.N1, fields[j], a[j], FFTW_REDFT01, 1);

        fftw_execute(p2);


        fftw_destroy_plan(p2);
    }
    CosineWeightBackward1D *clone() const
    {
        return new CosineWeightBackward1D(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};

class CosineWeightForward : public Weight<double, double>
{
private:
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_r2r_2d(p.N1, p.N2, fields[j], a[j], FFTW_REDFT10, FFTW_REDFT10, 1);

        fftw_execute(p2);
        double corr = 1. / (p.N1 * 4.);
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    CosineWeightForward *clone() const
    {
        return new CosineWeightForward(*this);
    }
    void print() { cout << "FourierWeightForwards" << endl; }
};

class CosineWeightBackward : public Weight<double, double>
{
private:
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        fftw_plan p2;

        p2 = fftw_plan_r2r_2d(p.N1, p.N2, fields[j], a[j], FFTW_REDFT01, FFTW_REDFT01, 1);

        fftw_execute(p2);
        int end = p.get_total();
        double corr = 1. / (p.N1);
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    CosineWeightBackward *clone() const
    {
        return new CosineWeightBackward(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};

class CosineWeightBackwardGradient : public Weight<double, double>
{

    // Given a value of the cosine transform, compute the real space gradient, this means applying an
    // inverse sin transform to (-k f(k)) where f(k) is our original function 
private:
int dir; //a value 0 or 1 depending on x or y;
public:
    void set_direction(int dirr) {
        if(dirr != 0 && dirr != 1) {
            error("fourier backwards gradient needs to be an integer of 0 or 1");
        }
        dir = dirr;
    }

    void multiplythrough(double *fields, const CH_builder &p) {
        //negative factors here arise from cosine/sin transform asymmetry
        if(dir == 0 ) {
            for(int i  = 0 ; i < p.N1 ; i++) {
                for (int j = 0; j < p.N2 ; j++) {
                    fields[i * p.N2 + j] = -i * fields[i * p.N2 + j];
                }
            }
        }
        else if( dir = 1) {

            for (int i = 0; i < p.N1; i++)
            {
                for (int j = 0; j < p.N2; j++)
                {
                    fields[i * p.N2 + j] = -j * fields[i * p.N2 + j];
                }
            }
        }
        else {
            error("error in Gradient specification in Fourier Gradient");
            /*error in through multiply*/
        }
    }
    void operator()(double **a, double **fields, int j, const CH_builder &p)
    { //warning, this operator destroys fields
        
        double *newfield;
        newfield = (double *)fftw_malloc(p.N1 * p.N2 * sizeof(double));
        int t = p.get_total();
        for(int i = 0  ; i < t ; i++ ) {
            newfield[i] = fields[j][i];
        }

        multiplythrough(newfield,p);
        
        fftw_plan p2;


        //inverse sin or cosine transform depending on which direction we are considering
        if(dir == 0 ) {
        p2 = fftw_plan_r2r_2d(p.N1, p.N2, newfield, a[j], FFTW_RODFT01, FFTW_REDFT01, 1);
        }
        else{
        p2 = fftw_plan_r2r_2d(p.N1, p.N2, newfield, a[j], FFTW_REDFT01, FFTW_RODFT01, 1);
        }
        fftw_execute(p2);
        int end = p.get_total();
        double corr = 1. / (p.N1);
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
        fftw_free(newfield);
    }
    CosineWeightBackwardGradient *clone() const
    {
        return new CosineWeightBackwardGradient(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};


class FourierWeightBackwardGradient2D : public Weight<complex<double>, complex<double> >
{

    // Given a value of the cosine transform, compute the real space gradient, this means applying an
    // inverse sin transform to (-k f(k)) where f(k) is our original function
private:
    int dir; //a value 0 or 1 depending on x or y;
public:
    void set_direction(int dirr)
    {
        if (dirr != 0 && dirr != 1)
        {
            error("fourier backwards gradient needs to be an integer of 0 or 1");
        }
        dir = dirr;
    }

        void multiplythrough(complex<double> *fields, const CH_builder &p)
    {
        //negative factors here arise from cosine/sin transform asymmetry
        if (dir == 0)
        {
            for (int i = 0; i < p.N1; i++)
            {
                for (int j = 0; j < p.N2; j++)
                {
                    // fftw_complex mod;
                    // mod[0] = i;
                    // mod[1] = j;
                    double k1;
                    if (i <= p.N1/2)
                    {
                        k1 = (2 * pii) * i;
                    }
                    else
                    {
                        k1 = (2 * pii) * (i -  1.0*p.N1);
                    }
                    // if (j <= p.N2/2)
                    // {
                    //     k2 = (2 * pii) * j;
                    // }
                    // else
                    // {
                    //     k2 = (2 * pii) * (j - p.N2);
                    // }

                    complex<double> mod= {0, k1 }; //gradient, so imaginary
                    fields[i*p.N2 +j] *= mod;
                    
                    //fields[i * p.N2 + j] =  fields[i * p.N2 + j];
                }
            }
        }
        else if (dir == 1)
        {

            for (int i = 0; i < p.N1; i++)
            {
                for (int j = 0; j < p.N2; j++)
                {
                   
                    double k2;
                    if (j <= p.N2/2)
                    {
                        k2 = (2 * pii) * j;
                    }
                    else
                    {
                        k2 = (2 * pii) * (j - 1.0*p.N2);
                    }
                    complex<double> mod = {0, k2}; //gradient, so imaginary



                    fields[i * p.N2 + j] *= mod;
                    
    
                    
                    //fields[i * p.N2 + j] =  fields[i * p.N2 + j];
                }
            }
        }
        else
        {
            error("error in Gradient specification in Fourier Gradient");
            /*error in through multiply*/
        }
    }

    //a is the output, fields is the input
    void operator()(complex<double> **a, complex<double> **fields, int j, const CH_builder &p)
    { //warning, this operator destroys fields

        complex<double> *newfield;
        newfield = (complex<double> *)fftw_malloc(p.N1 * p.N2 * sizeof(complex<double>));
        int t = p.get_total();
        for (int i = 0; i < t; i++)
        {
            newfield[i] = fields[j][i];
            
        }

        multiplythrough(newfield, p);

        fftw_plan p2;

        //inverse sin or cosine transform depending on which direction we are considering

        p2 = fftw_plan_dft_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *>(newfield), reinterpret_cast<fftw_complex *> (a[j]), FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p2);
        int end = p.get_total();
        double corr = 1. / (p.N1);
        for (int i = 0; i < end; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
        fftw_free(newfield);
    }
    FourierWeightBackwardGradient2D *clone() const
    {
        return new FourierWeightBackwardGradient2D(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};

#endif /* WEIGHTS_H */
