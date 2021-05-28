#ifndef WEIGHTS_H
#define WEIGHTS_H


class Weight {
    public:
        virtual void operator()(double **a, double **fields, int j, const CH_builder &p) = 0;
        virtual Weight *clone() const = 0;
        virtual void print() { cout << "the base Weight class\n"; }
};

class DiffusiveWeight : public Weight {
private:
    double eps;
public:
    DiffusiveWeight(double epss) : eps(epss) {}
    void operator()(double **a, double **fields, int j, const CH_builder &p) {
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

class DiffDiffusiveWeight : public Weight {
private:
    int which;
    int chi;
public:
    DiffDiffusiveWeight(double chii, int whichh) : which(whichh), chi(chii) {}
    void operator()(double **a, double **fields, int j, const CH_builder &p)
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

class DiffDiffusiveWeightSQR : public Weight {
private:
    int which;
    int chi;
public:
    DiffDiffusiveWeightSQR(double chii, int whichh) : which(whichh), chi(chii) {}
    void operator()(double **a, double **fields, int j, const CH_builder &p)
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

class CahnHilliardWeight : public Weight {
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

    void operator()(double **a, double **fields, int j, const CH_builder &p) {
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

class CahnHilliardWithCouplingWeight : public Weight {

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
    CahnHilliardWithCouplingWeight(double c00, double c11, double nuu, double ei11, double ei22, int which11, int which22)
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

    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = ei1 * (fields[which1][i]) + ei2 * (fields[which2][i])  + cons1 * CUB(fields[0][i]) + cons2 * SQR(fields[0][i]) + cons3 * (fields[0][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeight *clone() const
    {
        return new CahnHilliardWithCouplingWeight(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};

class CahnHilliardWithCouplingWeightSQR : public Weight
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

    void operator()(double **a, double **fields, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            a[j][i] = ei1 * SQR(fields[which1][i]) * (fields[j][i]) + ei2 * SQR(fields[which2][i]) * (fields[j][i]) + cons1 * CUB(fields[0][i]) + cons2 * SQR(fields[0][i]) + cons3 * (fields[0][i]) + cons4;
        }
    }
    CahnHilliardWithCouplingWeightSQR *clone() const
    {
        return new CahnHilliardWithCouplingWeightSQR(*this);
    }
    void print() { cout << "the Cahn Hilliard weight class\n"; }
};


class FourierWeightForward : public Weight {
private:
public:
    void operator()(double **a, double **fields, int j, const CH_builder &p) {
        fftw_plan p2;

        p2 = fftw_plan_r2r_2d(p.N1, p.N2, fields[j], a[j], FFTW_REDFT10, FFTW_REDFT10, 1);

        fftw_execute(p2);
        double corr = 1./(p.N1*4.);
        int end = p.get_total();
        for(int i = 0  ; i < end ; i++) {
            a[j][i]*=corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightForward *clone() const
    {
        return new FourierWeightForward(*this);
    }
    void print() {cout << "FourierWeightForwards" << endl;}
};

class FourierWeightBackward : public Weight
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
        for (int i = 0; i < end ; i++)
        {
            a[j][i] *= corr;
        }

        fftw_destroy_plan(p2);
    }
    FourierWeightBackward *clone() const
    {
        return new FourierWeightBackward(*this);
    }
    void print() { cout << "FourierWeightBackwards" << endl; }
};

#endif /* WEIGHTS_H */
