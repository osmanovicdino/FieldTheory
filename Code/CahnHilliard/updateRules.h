#ifndef UPDATERULES_H
#define UPDATERULES_H

struct updateRules {
CH_builder par;
    //void what_am_I() { cout << "I am the base class for Update Rules" << endl; }
// updateRules() {par.N1=1; par.N2=1; par.number_of_fields=1;}
updateRules(const CH_builder &para) : par(para) {}
virtual void operator()(double **, double**,double**,double**,int,const CH_builder &) = 0;
virtual updateRules *clone() const = 0;

};

struct NormalDiffusion : updateRules {
//double **calculated;


double dt;
double Di;
double *upd1;
double temp1;

    NormalDiffusion(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules(params) {
        Di = Dii;
        dt  = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for(int i = 0 ; i < params.N1 ; i++) {
            for(int j = 0  ; j < params.N2 ; j++) {
                upd1[i*params.N2 + j] = 1./(1. + dt*Di*temp1*(i*i+j*j));
            }
        }


    }

    NormalDiffusion(const NormalDiffusion &a) : updateRules(a.par) {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
            }
        }
    }

    NormalDiffusion &operator=(const NormalDiffusion &a)
    {
        fftw_free(upd1);
        
        par = a.par;
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
            }
        }
        return *this;
    }

    ~NormalDiffusion() {
            fftw_free(upd1);
    }

    void operator()(double **res, double **f, double **w, double **r, int j , const CH_builder &p) {
        int end = p.get_total();
        for(int i  = 0 ; i < end ; i++) {
            res[j][i] = upd1[i]*f[j][i]+dt*upd1[i]*r[j][i];
        }   
    }

    NormalDiffusion *clone() const {
        return new NormalDiffusion(*this);
    }
};

struct DiffusionWithInteraction : updateRules
{
    //double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;

    DiffusionWithInteraction(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                upd1[i * params.N2 + j] = dt * Di * temp1 * (i * i + j * j);
                upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * (i * i + j * j));
            }
        }
    }

    DiffusionWithInteraction(const DiffusionWithInteraction &a) : updateRules(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }
    }

    DiffusionWithInteraction &operator=(const DiffusionWithInteraction &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);

        par = a.par;
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }
        return *this;
    }

    ~DiffusionWithInteraction()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(double **res, double **f, double **w, double **r, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + dt * upd2[i] * r[j][i];
        }
    }

    DiffusionWithInteraction *clone() const
    {
        return new DiffusionWithInteraction(*this);
    }
};

struct DiffusionWithSurfaceTension : updateRules {
    double **calculated;

    // double **matup1;
    // double **matup2;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double eps;

    DiffusionWithSurfaceTension(const CH_builder &params, double dtt, double Dii, double temp11, double epss) : updateRules(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        eps = epss;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                upd1[i * params.N2 + j] =  dt * Di * temp1 * (i * i + j * j);
            }
        }
        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                double tempor = i * i + j * j;
                upd2[i * params.N2 + j] = 1./(1.+dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));

                // if (upd2[i * params.N2 + j] > 1950.) {
                //     cout << 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(i * i + j * j)) << endl;
                //     cout << Di << endl;
                //     cout << dt << endl;
                //     cout << eps << endl;
                //     cout << temp1 << endl;
                //     cout << i << endl;
                //     cout << j << endl;
                //     cout << SQR(i*i+j*j) << endl;
                //     pausel();
                // }
            }
        }
    }

    DiffusionWithSurfaceTension(const DiffusionWithSurfaceTension &a) : updateRules(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
            }
        }
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }
    }

    DiffusionWithSurfaceTension& operator=(const DiffusionWithSurfaceTension &a) {
        //cout << "called = " << endl;
        fftw_free(upd1);
        fftw_free(upd2);
        par = a.par;
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = dt * Di * temp1 * (i * i + j * j);
            }
        }
        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                double tempor = i * i + j * j;
                upd2[i * par.N2 + j] = 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));
                // upd2[i * par.N2 + j] = 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(i * i + j * j));

            }
        }
        return *this;

    }


    ~DiffusionWithSurfaceTension() {
        fftw_free(upd1);
        fftw_free(upd2);
        }

    void operator()(double **res, double **f, double **w, double **r, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] -upd2[i]*upd1[i]*w[j][i] + dt * upd2[i] * r[j][i];
            // if(res[j][i]>3000.) {
            //     cout << i << endl;
            //     cout << upd2[i] << endl;
            //     cout << upd1[i] << endl;
            //     cout << f[j][i] << endl;

            //     cout << w[j][i] << endl;
            //     cout << r[j][i] << endl;
            //     pausel();
            // }
        }
    }

    DiffusionWithSurfaceTension *clone() const
    {
        return new DiffusionWithSurfaceTension(*this);
    }
};





#endif /* UPDATERULES_H */
