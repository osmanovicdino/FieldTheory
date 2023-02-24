#ifndef UPDATERULES_CPP
#define UPDATERULES_CPP


template <class T, class T1, class T2, class T3>
struct NoRule : updateRules<T,T1,T2,T3> { // don't perform any update

    NoRule(const CH_builder &params) : updateRules<T,T1,T2,T3>(params)
    {
    }

    NoRule(const NoRule &a) : updateRules<T, T1, T2, T3>(a.par)
    {
    }

    NoRule& operator=(const NoRule &a) {
        this->setpar(a.par);
        return *this;
    }


    void operator()(T **res, T1 **f, T2 **w, T3 **r, int j, const CH_builder &p)
    {
        // do nothing
    }

    NoRule *clone() const
    {
        return new NoRule(*this);
    }

    void print()
    {
        cout << "No rule defined! Fix this!" << endl;
    }


};

template <class T>
struct NormalDiffusion : updateRules<T,T,T,T> {
//double **calculated;


double dt;
double Di;
double *upd1;
double temp1;

    NormalDiffusion(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules<T,T,T,T>(params) {
        Di = Dii;
        dt  = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for(int i = 0 ; i < params.N1 ; i++) {
            for(int j = 0  ; j < params.N2 ; j++) {

                double k1, k2;
                if (i <= params.N1 / 2)
                {
                    k1 = i;
                }
                else
                {
                    k1 = (i - params.N1);
                }
                if (j <= params.N2 / 2)
                {
                    k2 = j;
                }
                else
                {
                    k2 = (j - params.N2);
                }

                double tempor = SQR(k1) + SQR(k2);

                upd1[i*params.N2 + j] = 1./(1. + dt*Di*temp1*tempor);
            }
        }


    }

    NormalDiffusion(const NormalDiffusion &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
            }
        }
    }

    NormalDiffusion &operator=(const NormalDiffusion &a)
    {
        fftw_free(upd1);
        
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
            }
        }
        return *this;
    }

    ~NormalDiffusion() {
            fftw_free(upd1);
    }

    void operator()(T **res, T **f, T **w, T **r, int j , const CH_builder &p) {
        int end = p.get_total();
        for(int i  = 0 ; i < end ; i++) {
            res[j][i] = upd1[i]*f[j][i]+dt*upd1[i]*r[j][i];
        }   
    }

    // void trylower(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    // {
    //     int end = p.get_total();
    //     for (int i = 0; i < end; i++)
    //     {
    //         res[j][i] = upd1[i] * f[j][i] + dt * upd1[i] * r[j][i];
    //     }
    // }

    NormalDiffusion *clone() const {
        return new NormalDiffusion(*this);
    }

    void print() {
        cout << "Normal Diffusion" << endl;
    }


};

template <class T>
struct NormalDiffusion3D : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double temp1;

    NormalDiffusion3D(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * params.N3 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                for (int k = 0; k < params.N3; k++)
                {
                    double k1, k2, k3;
                    if (i <= params.N1 / 2)
                    {
                        k1 = i;
                    }
                    else
                    {
                        k1 = (i - params.N1);
                    }
                    if (j <= params.N2 / 2)
                    {
                        k2 = j;
                    }
                    else
                    {
                        k2 = (j - params.N2);
                    }
                    if (k <= params.N3 / 2)
                    {
                        k3 = k;
                    }
                    else
                    {
                        k3 = (k - params.N3);
                    }

                    double tempor = SQR(k1) + SQR(k2) + SQR(k3);

                    upd1[i * params.N2 * params.N3 + j*params.N3 + k] = 1. / (1. + dt * Di * temp1 * tempor);
                }
            }
        }
    }

    NormalDiffusion3D(const NormalDiffusion3D &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                    upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
    }

    NormalDiffusion3D &operator=(const NormalDiffusion3D &a)
    {
        fftw_free(upd1);

        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                    upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
        return *this;
    }

    ~NormalDiffusion3D()
    {
        fftw_free(upd1);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd1[i] * f[j][i] + dt * upd1[i] * r[j][i];
        }
    }

    // void trylower(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    // {
    //     int end = p.get_total();
    //     for (int i = 0; i < end; i++)
    //     {
    //         res[j][i] = upd1[i] * f[j][i] + dt * upd1[i] * r[j][i];
    //     }
    // }

    NormalDiffusion3D *clone() const
    {
        return new NormalDiffusion3D(*this);
    }

    void print()
    {
        cout << "Normal Diffusion" << endl;
    }
};

template<class T>
struct DiffusionWithInteraction : updateRules<T,T,T,T>
{
    //double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;

    DiffusionWithInteraction(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules<T,T,T,T>(params)
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
                double k1, k2;
                if (i <= params.N1 / 2)
                {
                    k1 = i;
                }
                else
                {
                    k1 = (i - params.N1);
                }
                if (j <= params.N2 / 2)
                {
                    k2 = j;
                }
                else
                {
                    k2 = (j - params.N2);
                }

                double tempor = SQR(k1) + SQR(k2);

                upd1[i * params.N2 + j] = dt * Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
            }
        }
    }

    DiffusionWithInteraction(const DiffusionWithInteraction &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }
    }

    DiffusionWithInteraction &operator=(const DiffusionWithInteraction &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }

        return *this;
    }

    ~DiffusionWithInteraction()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
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
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct DiffusionWithInteraction3D : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;

    DiffusionWithInteraction3D(const CH_builder &params, double dtt, double Dii, double temp11) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * params.N3 * sizeof(double));
        upd2 = (double *)fftw_malloc(params.N1 * params.N2 * params.N3 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                for(int k = 0 ;  k < params.N3 ; k++) {
                    double k1, k2, k3;
                    if (i <= params.N1 / 2)
                    {
                        k1 = i;
                    }
                    else
                    {
                        k1 = (i - params.N1);
                    }
                    if (j <= params.N2 / 2)
                    {
                        k2 = j;
                    }
                    else
                    {
                        k2 = (j - params.N2);
                    }
                    if (k <= params.N3 / 2)
                    {
                        k3 = k;
                    }
                    else
                    {
                        k3 = (k - params.N3);
                    }

                    double tempor = SQR(k1) + SQR(k2) + SQR(k3);

                    upd1[i * params.N2 * params.N3 + j * params.N3 + k] = dt * Di * temp1 * tempor;
                    upd2[i * params.N2 * params.N3 + j * params.N3 + k] = 1. / (1. + dt * Di * temp1 * tempor);
                }
            }
        }
    }

    DiffusionWithInteraction3D(const DiffusionWithInteraction3D &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                    upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                    upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
    }

    DiffusionWithInteraction3D &operator=(const DiffusionWithInteraction3D &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for(int k = 0  ; k < a.par.N3 ; k++) {
                    upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                    upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }

        return *this;
    }

    ~DiffusionWithInteraction3D()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + dt * upd2[i] * r[j][i];
        }
    }

    DiffusionWithInteraction3D *clone() const
    {
        return new DiffusionWithInteraction3D(*this);
    }
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct DiffusionWithSurfaceTension : updateRules<T,T,T,T>
{
    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double eps;

    DiffusionWithSurfaceTension(const CH_builder &params, double dtt, double Dii, double temp11, double epss) : updateRules<T,T,T,T>(params)
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

                double k1, k2;
                if (i <= params.N1 / 2)
                {
                    k1 = i;
                }
                else
                {
                    k1 =  (i - params.N1);
                }
                if (j <= params.N2 / 2)
                {
                    k2 = j;
                }
                else
                {
                    k2 = (j - params.N2);
                }

                double tempor = SQR(k1)+SQR(k2);

                upd1[i * params.N2 + j] =  dt * Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1./(1.+dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));

            }
        }
        // for (int i = 0; i < params.N1; i++)
        // {
        //     for (int j = 0; j < params.N2; j++)
        //     {
        //         double tempor = i * i + j * j;
        //         upd2[i * params.N2 + j] = 1./(1.+dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));

        //         // if (upd2[i * params.N2 + j] > 1950.) {
        //         //     cout << 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(i * i + j * j)) << endl;
        //         //     cout << Di << endl;
        //         //     cout << dt << endl;
        //         //     cout << eps << endl;
        //         //     cout << temp1 << endl;
        //         //     cout << i << endl;
        //         //     cout << j << endl;
        //         //     cout << SQR(i*i+j*j) << endl;
        //         //     pausel();
        //         // }
        //     }
        // }
    }

    DiffusionWithSurfaceTension(const DiffusionWithSurfaceTension &a) : updateRules<T,T,T,T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
            }
        }
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }
    }

    DiffusionWithSurfaceTension& operator=(const DiffusionWithSurfaceTension &a) {
        //cout << "called = " << endl;
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
            }
        }
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }

        return *this;

    }


    ~DiffusionWithSurfaceTension() {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
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
    void print()
    {
        cout << "Diffusion With Surface Tension" << endl;
    }
};

template <class T>
struct DiffusionWithSurfaceTension3D : updateRules<T, T, T, T>
{
    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double eps;

    DiffusionWithSurfaceTension3D(const CH_builder &params, double dtt, double Dii, double temp11, double epss) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        eps = epss;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * params.N3 * sizeof(double));
        upd2 = (double *)fftw_malloc(params.N1 * params.N2 * params.N3 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                for (int j = 0; j < params.N2; j++)
                {
                    for(int k  = 0 ; k < params.N3 ; k++) 
                    {

                        double k1, k2, k3;
                        if (i <= params.N1 / 2)
                        {
                            k1 = i;
                        }
                        else
                        {
                            k1 = (i - params.N1);
                        }
                        if (j <= params.N2 / 2)
                        {
                            k2 = j;
                        }
                        else
                        {
                            k2 = (j - params.N2);
                        }
                        if (k <= params.N3 / 2)
                        {
                            k3 = k;
                        }
                        else
                        {
                            k3 = (k - params.N3);
                        }

                        double tempor = SQR(k1) + SQR(k2) + SQR(k3);

                        upd1[i * params.N2 * params.N3 + j * params.N3 + k] = dt * Di * temp1 * tempor;
                        upd2[i * params.N2 * params.N3 + j * params.N3 + k] = 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));
                    }
                }
            }
        }
        // for (int i = 0; i < params.N1; i++)
        // {
        //     for (int j = 0; j < params.N2; j++)
        //     {
        //         double tempor = i * i + j * j;
        //         upd2[i * params.N2 + j] = 1./(1.+dt * Di * SQR(eps) * SQR(temp1) * SQR(tempor));

        //         // if (upd2[i * params.N2 + j] > 1950.) {
        //         //     cout << 1. / (1. + dt * Di * SQR(eps) * SQR(temp1) * SQR(i * i + j * j)) << endl;
        //         //     cout << Di << endl;
        //         //     cout << dt << endl;
        //         //     cout << eps << endl;
        //         //     cout << temp1 << endl;
        //         //     cout << i << endl;
        //         //     cout << j << endl;
        //         //     cout << SQR(i*i+j*j) << endl;
        //         //     pausel();
        //         // }
        //     }
        // }
    }

    DiffusionWithSurfaceTension3D(const DiffusionWithSurfaceTension3D &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * a.par.N3 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for(int k = 0 ; k < a.par.N3 ; k++) {
                upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
    }

    DiffusionWithSurfaceTension3D &operator=(const DiffusionWithSurfaceTension3D &a)
    {
        // cout << "called = " << endl;
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd1[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                for (int k = 0; k < a.par.N3; k++)
                {
                upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k] = a.upd2[i * a.par.N2 * a.par.N3 + j * a.par.N3 + k];
                }
            }
        }

        return *this;
    }

    ~DiffusionWithSurfaceTension3D()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + dt * upd2[i] * r[j][i];
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

    DiffusionWithSurfaceTension3D *clone() const
    {
        return new DiffusionWithSurfaceTension3D(*this);
    }
    void print()
    {
        cout << "Diffusion With Surface Tension" << endl;
    }
};

// struct CompleteRuleFractional
// {
//     Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > fractional_update;
//     Rule_Wrapper<complex<double>, complex<double>, complex<double>, complex<double> > fractional_weight;


// };

// struct MultiPhase : updateRules {

//     int no_pha;
//     double dt;
//     matrix<double> Chi;
//     vector1<double> Di; //vector of diffusion constants
//     double *upd1;
//     double *upd2;
//     double **upd3;
//     double temp1;
//     double surface_width;

//     MultiPhase(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, double surface_widthh, matrix<double> Chii) : updateRules(params), Di(Dii), Chi(Chii)
//     {
//         no_pha = no_phaa;
//         dt = dtt;
//         temp1 = temp11;
//         surface_width = surface_widthh;

//         upd3 = new double * [no_pha];

//         upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
//         upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
//         for(int j = 0  ; j < no_pha ; j++)
//         upd3[j] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

//         for (int i = 0; i < params.N1; i++)
//         {
//             for (int j = 0; j < params.N2; j++)
//             {
//                 upd1[i * params.N2 + j] = dt * Di[0] * temp1 * (i * i + j * j);
//             }
//         }
//         for (int i = 0; i < params.N1; i++)
//         {
//             for (int j = 0; j < params.N2; j++)
//             {
//                 double tempor = i * i + j * j;
//                 upd2[i * params.N2 + j] = dt * Di[0] * SQR(surface_width) * SQR(temp1) * SQR(tempor);

//             }
//         }

//         double vol = 0.0;
//         Chi.maxima(vol);
//         if(vol < 1E-12)
//         vol = 1.0;
//         vol *= 0.5;

//         for(int k = 0 ; k < no_pha+1 ; k++) {
//             Chi(k, k) = vol;
//         }

//         for(int k = 0  ; k < no_pha ; k++) {
//         for (int i = 0; i < params.N1; i++)
//         {
//             for (int j = 0; j < params.N2; j++)
//             {
//                 double tempor = i * i + j * j;
//                 upd3[k][i * params.N2 + j] = 1. / (1. + dt * (vol + Chi(k,no_pha)) * Di[k] * SQR(surface_width) * SQR(temp1) * SQR(tempor));
//             }
//         }
//         }

//     }

//     MultiPhase(const MultiPhase &a) : updateRules(a.par), Di(a.Di), Chi(a.Chi)
//     {
//         no_pha = a.no_pha;
//         surface_width = a.surface_width;
//         dt = a.dt;
//         temp1 = a.temp1;

//         upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double) );
//         upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

//         upd3 = new double *[no_pha];
//         for (int j = 0; j < no_pha; j++)
//             upd3[j] = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

//         for (int i = 0; i < par.N1; i++)
//         {
//             for (int j = 0; j < par.N2; j++)
//             {
//                 upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
//                 upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
//                 for (int k = 0; k < no_pha; k++)
//                 {
//                     upd3[k][i * par.N2 + j] = a.upd3[k][i * par.N2 + j];
//                 }
//             }
//         }

//     }
//     MultiPhase &operator=(const MultiPhase &a)
//     {
//         fftw_free(upd1);
//         fftw_free(upd2);
//         for (int i = 0; i < no_pha; i++)
//         {
//             fftw_free(upd3[i]);
//         }
//         delete upd3;
//         par = a.par;
//         Di = a.Di;
//         Chi = a.Chi;
//         no_pha = a.no_pha;
//         surface_width = a.surface_width;
//         dt = a.dt;
//         temp1 = a.temp1;

//         upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
//         upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
//         upd3 = new double *[no_pha];
//         for (int j = 0; j < no_pha; j++)
//             upd3[j] = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

//         for (int i = 0; i < par.N1; i++)
//         {
//             for (int j = 0; j < par.N2; j++)
//             {
//                 upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
//                 upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
//                 for (int k = 0; k < no_pha; k++)
//                 {
//                     upd3[k][i * par.N2 + j] = a.upd3[k][i * par.N2 + j];
//                 }
//             }
//         }

//         return *this;
//     }

//     ~MultiPhase()
//     {
//         fftw_free(upd1);
//         fftw_free(upd2);
//         for (int i = 0; i < no_pha; i++)
//         {
//             fftw_free(upd3[i]);
//         }
//         delete upd3;
//     }

//     void operator()(double **res, double **f, double **w, double **r, int j, const CH_builder &p)
//     {
//         //cout << "done rule " << j << endl;
//         int end = p.get_total();
//         for (int i = 0; i < end; i++)
//         {
//             double val  = 0.0;

//             for(int k = 0  ; k < no_pha ; k++ ) {
//                 // if(k == no_pha) {
//                 //     double valt = 0.0;
//                 //     for(int k2 = 0 ; k2 < no_pha; k2++) {
//                 //         valt += f[k2][i];
//                 //     }
//                 //     val += Chi(j,k) * (-valt);
//                 // }
//                 double exval = 0.0;
//                 if(k != j) exval += Chi(j,no_pha);
//                 val += (Chi(j,k) - exval)*f[k][i];
//             }

//             res[j][i] = upd3[j][i] * f[j][i] - upd3[j][i] * upd1[i] * (Di[j] / Di[0]) * w[j][i] + upd3[j][i] * upd2[i] * (Di[j] / Di[0]) * val + dt * upd3[j][i] * r[j][i];

// /*             if(   abs((res[j][i]-f[j][i])/f[j][i])>100.0 && j > 0 )
//              {
//                 if ( abs((res[j][i]-f[j][i])/f[j][i]) >0.1 ) {
//                     cout << Chi << endl;
//                     cout << res[j][i] << endl;
//                     cout << f[j][i] << endl;
//                     cout << j << " " << i << endl;
//                     cout << upd1[i] << endl;
//                     cout << upd2[i] << endl;
//                     cout << upd3[j][i] << endl;
//                     cout << w[j][i] << endl;
//                     cout << val << endl;
//                     for(int k = 0 ; k < no_pha ; k++) {
//                         cout << f[k][i] << " ";
//                     }
//                     cout << endl;
//                     cout << r[j][i] << endl;
//                     cout << endl;
//                     cout << upd3[j][i] * upd1[i] * (Di[j] / Di[0]) * w[j][i] << endl;
//                     cout << upd3[j][i] * upd2[i] * (Di[j] / Di[0]) * val << endl;
//                     pausel();
//                 }
//             } */
//             // if(res[j][i]>3000.) {
//             //     cout << i << endl;
//             //     cout << upd2[i] << endl;
//             //     cout << upd1[i] << endl;
//             //     cout << f[j][i] << endl;

//             //     cout << w[j][i] << endl;
//             //     cout << r[j][i] << endl;
//             //     pausel();
//             // }
//         }

//     }

//     MultiPhase *clone() const
//     {
//         return new MultiPhase(*this);
//     }
// };

struct MultiPhase : updateRules<double, double, double, double>
{

    int no_pha;
    double dt;
    matrix<double> Chi;
    vector1<double> Di; //vector of diffusion constants
    double *upd1;
    //double *upd2;
    matrix<double> *upd3;
    //vector1<int> *upd2;
    double temp1;
    matrix<double> surface_width;

    int tot_matrices;

    MultiPhase(const CH_builder &params) : updateRules(params), surface_width(matrix<double>(1,1)), Di(vector1<double>(1)), Chi(matrix<double>(1, 1))
    {
        no_pha = 1;
        dt = 1.;
        temp1 = 1.;
        //surface_width = 1.;
        tot_matrices = 1;

        //upd3 = new double *[no_pha];
        upd3 = new matrix<double>[tot_matrices];
        //upd2 = new vector1<int>[params.N1 * params.N2];

        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        matrix<double> temp(3,3);
        upd3[0] = temp; // = dt * Di[0] * SQR(surface_width) * SQR(temp1) * SQR(tempor);

    }

    

    MultiPhase(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, matrix<double> surface_widthh, matrix<double> Chii) : updateRules(params), Di(Dii), Chi(Chii), surface_width(surface_widthh)
    {
        no_pha = no_phaa;
        dt = dtt;
        temp1 = temp11;
        //surface_width = surface_widthh;
        tot_matrices = params.N1 * params.N2;

        //upd3 = new double *[no_pha];
        upd3 = new matrix<double>[tot_matrices];
        //upd2 = new vector1<int>[params.N1 * params.N2];

        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        
        

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                upd1[i * params.N2 + j] = dt * Di[0] * temp1 * (i * i + j * j);
            }
        }
        int nr = Chi.getnrows() -1;
        int nc = Chi.getncols() -1;

        double min;
        Chi.minima(min);

        matrix<double> test(nr, nc);

        for (int i1 = 0; i1 < nr; i1++)
        {
            for (int j1 = 0; j1 < nc; j1++)
            {

                test(i1, j1) = SQR(surface_width(i1, j1)) * (Chi(i1, j1) - Chi(i1, no_pha));
            }
        }

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                double tempor = i * i + j * j;

                matrix<double> comb(nr, nc);

                for(int i1 = 0 ; i1 < nr ; i1++ ) {
                    for(int j1 = 0 ; j1 < nc ; j1++) {
                        
                        comb(i1,j1) = KroneckerDelta(i1,j1) - Di[i1]*SQR(surface_width(i1,j1))*SQR(temp1)*SQR(tempor)*dt*(Chi(i1,j1)-Chi(i1,no_pha) );
                    }
                }
                comb.inverse();
                //pausel();
                
                //upd2[i * params.N2 + j] = indx;
                upd3[i * params.N2 + j] = comb;// = dt * Di[0] * SQR(surface_width) * SQR(temp1) * SQR(tempor);

            }
        }



        // cout << Di[1] * SQR(surface_width) * SQR(temp1) * dt << endl;;

        // double vol = 0.0;
        // Chi.maxima(vol);
        // if (vol < 1E-12)
        //     vol = 1.0;
        // vol *= 0.5;

        // for (int k = 0; k < no_pha + 1; k++)
        // {
        //     Chi(k, k) = vol;
        // }

        // for (int k = 0; k < no_pha; k++)
        // {
        //     for (int i = 0; i < params.N1; i++)
        //     {
        //         for (int j = 0; j < params.N2; j++)
        //         {
        //             double tempor = i * i + j * j;
        //             upd3[k][i * params.N2 + j] = 1. / (1. + dt * (vol + Chi(k, no_pha)) * Di[k] * SQR(surface_width) * SQR(temp1) * SQR(tempor));
        //         }
        //     }
        // }
    }

    MultiPhase(const MultiPhase &a) : updateRules(a.par), Di(a.Di), Chi(a.Chi), surface_width(a.surface_width)
    {
        no_pha = a.no_pha;
        dt = a.dt;
        temp1 = a.temp1;

        tot_matrices = a.tot_matrices;

        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));
        //upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

        //upd2 = new vector1<int>[par.N1 * par.N2];

        upd3 = new matrix<double>[tot_matrices];


        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                //upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }
        for(int k = 0  ; k < tot_matrices ; k++) {
            upd3[k] = a.upd3[k];

        }
    }
    MultiPhase &operator=(const MultiPhase &a)
    {
        fftw_free(upd1);
        delete upd3;
        //delete upd2;
        par = a.par;
        Di = a.Di;
        Chi = a.Chi;
        no_pha = a.no_pha;
        surface_width = a.surface_width;
        dt = a.dt;
        temp1 = a.temp1;
        tot_matrices = a.tot_matrices;

        upd1 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

       // upd2 = new vector1<int>[par.N1 * par.N2];
        upd3 = new matrix<double>[tot_matrices];


        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                //upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }
        for (int k = 0; k < tot_matrices; k++)
        {
            upd3[k] = a.upd3[k];

        }
        return *this;
    }

    ~MultiPhase()
    {
        fftw_free(upd1);
        //delete upd2;
        delete upd3;
    }

    void operator()(double **res, double **f, double **w, double **r, int j, const CH_builder &p)
    {
        //cout << "done rule " << j << endl;
        if(j > 0) {
            //we skip all the ones where j is greater than zero as this entire problem is to be solved collectively

        }
        else{
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {

            vector1<double> b(no_pha);
            // vector1<double> c(no_pha);

            // for(int k  = 0 ; k < no_pha ; k++)
            //     c[k] = f[k][i];

            // c = upd4[i] * c;

            for(int k = 0 ; k < no_pha ; k++) {
                b[k] = f[k][i]  - upd1[i] *  (Di[k]) * w[k][i] + dt * r[k][i];
            }



            //vector1<double> c = b;

            b = upd3[i] * b;

            //lubksb(upd3[i],upd2[i],b);
            for (int k = 0; k < no_pha; k++)
            {
                res[k][i] = b[k];
            }
/* 
            int i1 = floor(i / par.N1);
            int j1 = i % (par.N1);
            if(i1 == 0 && j1 == 10) {
                cout << upd3[i] << endl;
                cout << upd1[i] << endl;

                for (int k = 0; k < no_pha; k++)
                    cout << f[k][i] << " ";
                cout << endl;
                for (int k = 0; k < no_pha; k++)
                    cout << w[k][i] << " ";
                cout << endl;

                for(int k = 0 ; k < no_pha ; k++)
                    cout << res[k][i] << " ";
                cout << endl;


                pausel();
            } */
            /*             int i1 = floor(i / par.N1);
            int j1 = i % (par.N1);

             if(i1 == 0 && j1 == 32) {
                cout << upd3[i] << endl;
                cout << upd1[i] << endl;
                cout << i << endl;
                cout << b << endl;
                int i1 = floor(i / par.N1);
                int j1 = i % (par.N1);
                int nr = Chi.getnrows() - 1;
                int nc = nr;
                matrix<double> comb(nr, nc);
                double tempor = i1 * i1 + j1 * j1;
                for (int i2 = 0; i2 < nr; i2++)
                {
                    for (int j2 = 0; j2 < nc; j2++)
                    {
                        
                        comb(i2, j2) = KroneckerDelta(i2, j2) - Di[i2] * SQR(surface_width(i2,j2)) * SQR(temp1) * SQR(tempor) * dt * (Chi(i2, j2) - Chi(i2, no_pha));
                    }
                }
                cout << endl;
                cout << comb << endl;

                cout << SQR(temp1) * SQR(tempor) * dt << endl;
                for (int k = 0; k < no_pha; k++)
                {
                    cout << f[k][i] << ",";
                }
                cout << endl;
                for (int k = 0; k < no_pha; k++)
                {
                    cout << upd1[i]*w[k][i] << ",";
                }
                cout << endl;
                pausel();
            } */

            /*             for(int k = 0 ; k < no_pha ; k++) {
                if(abs((b[k]-c[k])/c[k])>100.0) {
                    cout << upd3[i] << endl;
                    cout << i << endl;
                    cout << b << endl;
                    cout << c << endl;
                    int i1 = floor(i / par.N1);
                    int j1 = i % (par.N1);
                    int nr = Chi.getnrows()-1;
                    int nc = nr;
                    matrix<double> comb(nr, nc);

                    for (int i1 = 0; i1 < nr; i1++)
                    {
                        for (int j1 = 0; j1 < nc; j1++)
                        {
                            double tempor = i * i + j * j;
                            comb(i1, j1) = KroneckerDelta(i1, j1) - Di[i1] * SQR(surface_width) * SQR(temp1) * SQR(tempor) * dt * (Chi(i1, j1) - Chi(i1, no_pha));
                        }
                    }
                    cout << endl;
                    cout << comb << endl;
                    for (int k = 0; k < no_pha; k++)
                    {
                        cout << f[k][i] << ",";
                    }
                    cout << endl;
                    for (int k = 0; k < no_pha; k++)
                    {
                        cout << w[k][i] << ",";
                    }
                    cout << endl;
                    pausel();
                }
            } */
            
            //res[j][i] = upd3[j][i] * f[j][i] - upd3[j][i] * upd1[i] * (Di[j] / Di[0]) * w[j][i] + upd3[j][i] * upd2[i] * (Di[j] / Di[0]) * val + dt * upd3[j][i] * r[j][i];
            }
        }
    }

    MultiPhase *clone() const
    {
        return new MultiPhase(*this);
    }

    void print()
    {
        cout << "MultiPhase" << endl;
    }
};

struct MultiPhaseKR : updateRules<complex<double>, complex<double> , complex<double>, complex<double> >
{

    int no_pha;
    double dt;
    matrix<double> Chi;
    vector1<double> Di; //vector of diffusion constants
    complex<double> *upd1;
    Field_Wrapper<complex<double>, complex<double> > revk1; //this will be a complex to real transform
    Field_Wrapper<complex<double>, complex<double> > revk2;

    Field_Wrapper<complex<double>, complex<double> > realspace1;
    Field_Wrapper<complex<double>, complex<double> > realspace2; //do not calculate with this one

    Field_Wrapper<complex<double>, complex<double> > revf1;
    Field_Wrapper<complex<double>, complex<double> > revf2;
    //double *upd2

    //vector1<int> *upd2;
    double temp1;
    matrix<double> surface_width;

    MultiPhaseKR(const CH_builder &params) : updateRules(params), revk1(params), revk2(params), realspace1(params), realspace2(params), revf1(params), revf2(params), surface_width(matrix<double>(1, 1)), Di(vector1<double>(1)), Chi(matrix<double>(1, 1))
    {
        no_pha = 1;
        dt = 1.;
        temp1 = 1.;
        //surface_width = 1.;

        //upd3 = new double *[no_pha];
        //upd2 = new vector1<int>[params.N1 * params.N2];

        upd1 = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
    }

    MultiPhaseKR(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, matrix<double> surface_widthh, matrix<double> Chii) : updateRules(params), revk1(params), revk2(params), revf1(params), revf2(params), realspace1(params), realspace2(params), Di(Dii), Chi(Chii), surface_width(surface_widthh)
    {
        if (Di.getsize() != no_phaa)
            error("invalid diffusion constant in MultiPhaseKR");
        if (Chii.getnrows() != no_phaa || Chii.getncols() != no_phaa)
            error("invalid Chi matrix in MultiPhaseKR");
        no_pha = no_phaa;
        dt = dtt;
        temp1 = temp11;
        //surface_width = surface_widthh;

        //upd3 = new double *[no_pha];
        //upd2 = new vector1<int>[params.N1 * params.N2];

        upd1 = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for (int i1 = 0; i1 < params.N1; i1++)
        {
            for (int j1 = 0; j1 < params.N2; j1++)
            {
                double k1, k2;
                if (i1 <= params.N1/2)
                {
                    k1 = (2 * pii) * i1;
                }
                else
                {
                    k1 = (2 * pii) * (i1 - params.N1);
                }
                if (j1 <= params.N2/2)
                {
                    k2 = (2 * pii) * j1;
                }
                else
                {
                    k2 = (2 * pii) * (j1 - params.N2);
                }

                upd1[i1 * params.N2 + j1]=  { -SQR(temp1) * SQR(2*pii)* (SQR(k1)+SQR(k2) ), 0.0 }; //fourier transform of k^2 operator
                
            }
        }
        int nr = Chi.getnrows() - 1;
        int nc = Chi.getncols() - 1;

        double min;
        Chi.minima(min);

        FourierWeightBackwardGradient2D a1;
        a1.set_direction(0);
        FourierWeightBackwardGradient2D a2;
        a2.set_direction(1);

        for(int k = 0 ; k < no_pha ; k++) {
            revk1.add_method(a1, k);
            revk2.add_method(a2, k);
        }

        FourierWeightForward2D a3;
        for (int k = 0; k < no_pha; k++)
        {
            revf1.add_method(a3, k);
            revf2.add_method(a3, k);
        }
    }

    MultiPhaseKR(const MultiPhaseKR &a) : updateRules(a.par), Di(a.Di), Chi(a.Chi), revk1(a.revk1), revk2(a.revk2), revf1(a.revf1), revf2(a.revf2), realspace1(a.realspace1), realspace2(a.realspace2), surface_width(a.surface_width)
    {
        no_pha = a.no_pha;
        dt = a.dt;
        temp1 = a.temp1;

        upd1 = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));
        //upd2 = (double *)fftw_malloc(par.N1 * par.N2 * sizeof(double));

        //upd2 = new vector1<int>[par.N1 * par.N2];


        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];

                //upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }

    }
    MultiPhaseKR &operator=(const MultiPhaseKR &a)
    {
        fftw_free(upd1);
        //delete upd2;
        par = a.par;
        Di = a.Di;
        Chi = a.Chi;
        no_pha = a.no_pha;
        surface_width = a.surface_width;
        dt = a.dt;
        temp1 = a.temp1;

        upd1 = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));

        revk1 = a.revk1;
        revk2 = a.revk2;

        revf1 = a.revf1;
        revf2 = a.revf2;

        realspace1 = a.realspace1;
        realspace2 = a.realspace2;

        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j]= a.upd1[i * par.N2 + j];
                //upd2[i * par.N2 + j] = a.upd2[i * par.N2 + j];
            }
        }

        return *this;
    }

    ~MultiPhaseKR()
    {
        fftw_free(upd1);
        //delete upd2;
    }

    void operator()(complex<double> **res, complex<double> **fr, complex<double> **w, complex<double> **fk, int j, const CH_builder &p)
    {
        // w is k space
        //f is real space in this framing

        //cout << "done rule " << j << endl;
        if(j>0) {
            //everything is done collectively in this updateRule
        }
        else{

            //the weights are calculated in k space,
            //all of our revs are calculated

            //this should now be real space gradient
            //cout << "older" << fk[0][300000] << endl;
            //we add the effect of the surface tension manually:

            int end = p.get_total();


            for(int k = 0  ; k < no_pha ; k++) {
                
                    for(int i  = 0 ; i < end ; i++) {
                        double val = 0.0;
                        for (int l = 0; l < no_pha; l++)
                        {
                            //the -1. arises from the cosine transform of the laplacian
                            val += SQR(surface_width(k, l)) * Chi(k, l) * (upd1[i].real() * fk[l][i].real());
                        }
                        w[k][i] += val;
                }
            } //we add the fourier weight here by hand
            



            // cout << "done 1" << endl;
            // string filename1 = "xrealspacex";


            revk1.Calculate_Results(w);
            revk2.Calculate_Results(w); //from the real space results, compute 

            // string filename1 = "xrealspacex";
            // string filename2 = "xrealspacey";
            // cout << "done 2" << endl;

            for(int k = 0 ; k < no_pha ; k++) {
                for (int i = 0; i < end; i++) {
                    double val1 =0.0;
                    double val2 =0.0;
                    for (int l = 0; l < no_pha; l++)
                    {
                        val1 += Di[k] * fr[k][i].real() * (KroneckerDelta(k, l) - fr[l][i].real()) * revk1.calculated_reactions[l][i].real();
                        val2 += Di[k] * fr[k][i].real() * (KroneckerDelta(k, l) - fr[l][i].real()) * revk2.calculated_reactions[l][i].real();
                    }
                    
                    realspace1.calculated_reactions[k][i] += val1;
                    realspace2.calculated_reactions[k][i] += val2;
                }
            }

            // cout << "done 3" << endl;





            revf1.Calculate_Results(realspace1.calculated_reactions);
            revf2.Calculate_Results(realspace2.calculated_reactions);

            //without a weight, these should be the 2*pi*k*orig;




            // cout << "done 4" << endl;

            for(int k = 0  ; k < no_pha ; k++)
            {
                for(int i1 = 0 ; i1 < p.N1 ; i1++)
                {
                    for (int j1 = 0; j1 < p.N2; j1++)
                    {

                        double k1, k2;
                        if (i1 <= p.N1/2)
                        {
                            k1 = (2 * pii) * i1;
                        }
                        else
                        {
                            k1 = (2 * pii) * (i1 - p.N1 );
                        }
                        if (j1 <= p.N2/2)
                        {
                            k2 = (2 * pii) * j1;
                        }
                        else
                        {
                            k2 = (2 * pii) * (j1 - p.N2 );
                        }
                        complex<double> myfac1 = {0.0,  k1};
                        complex<double> myfac2 = {0.0,  k2};

                        //we need a negative sign here for the two differentials we are doing

                        res[k][i1 * p.N2 + j1] = SQR(temp1) * SQR(2 * pii)  *(myfac1 * revf1.calculated_reactions[k][i1 * p.N2 + j1] + myfac2 * revf2.calculated_reactions[k][i1 * p.N2 + j1]);
                    }
                }
            }


            // int ss2 = 0;
            // int ds2 = 0;
            // for (int k = 0; k < no_pha; k++)
            // {
            //     for (int i = 0; i < end; i++)
            //     {
            //         if (k == 0)
            //         {
            //             if (sign(res[k][i].real()) == sign(fk[k][i].real()))
            //             {
            //                 ss2 += 1;
            //             }
            //             else
            //                 ds2 += 1;
            //         }
            //     }
            // }

            // cout << ss2 << endl;
            // cout << ds2 << endl;
            // pausel();

            // cout << "done 5" << endl;
            // cout << endl;

        }
    }

    MultiPhaseKR *clone() const
    {
        return new MultiPhaseKR(*this);
    }

    void print()
    {
        cout << "MultiPhaseKR" << endl;
    }
};

struct MultiPhaseKR_Update : updateRules<complex<double>, complex<double>, complex<double>, complex<double> > //use results of above class;
{

    int no_pha;
    double dt;
    matrix<double> Chi;
    vector1<double> Di; //vector of diffusion constants
    complex<double> *upd1;
    complex<double> **upd2;
    //double *upd2
    double A;
    //vector1<int> *upd2;
    double temp1;
    double surface_width;

    MultiPhaseKR_Update(const CH_builder &params) : updateRules(params), Di(vector1<double>(1)), Chi(matrix<double>(1, 1))
    {
        no_pha = 1;
        dt = 1.;
        temp1 = 1.;
        surface_width = 1.;
        double val;
        Chi.maxima(val);
        if(val < 1E-10)
            A= 0.5;
        else 
            A = 0.5*val;
        //upd3 = new double *[no_pha];
        //upd2 = new vector1<int>[params.N1 * params.N2];

        upd1 = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
        upd2 = new complex<double> *[no_pha];
        for (int j = 0; j < no_pha; j++)
            upd2[j] = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
    }

    MultiPhaseKR_Update(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, double surface_widthh, matrix<double> Chii) : updateRules(params),  Di(Dii), Chi(Chii)
    {
        if(Di.getsize() != no_phaa)
            error("invalid diffusion constant in MultiPhaseKR_Update");
        if (Chii.getnrows() != no_phaa || Chii.getncols() != no_phaa)
            error("invalid Chi matrix in MultiPhaseKR_Update");
        no_pha = no_phaa;
        dt = dtt;
        temp1 = temp11;
        surface_width = surface_widthh;
        double val;
        Chi.maxima(val);
        if (val < 1E-10)
            A = 0.5;
        else
            A = 0.5 * val;

        //         upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        //         upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        //         
        //         upd3[j] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        upd1 = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
        upd2 = new complex<double> *[no_pha];
        for (int j = 0; j < no_pha; j++)
            upd2[j] = (complex<double> *)fftw_malloc(params.N1 * params.N2 * sizeof(complex<double>));
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for (int i = 0; i < params.N1; i++)
        {
            for (int j = 0; j < params.N2; j++)
            {
                double k1, k2;
                if (i <= params.N1/2)
                {
                    k1 = (2 * pii) * i;
                }
                else
                {
                    k1 = (2 * pii) * (i - params.N1);
                }
                if (j <= params.N2/2)
                {
                    k2 = (2 * pii) * j;
                }
                else
                {
                    k2 = (2 * pii) * (j - params.N2);
                }
                
                complex<double> tempor = {SQR(SQR(temp1))*SQR(SQR(k1)+SQR(k2)) , 0 };

                

                upd1[i * params.N2 + j] = dt * SQR(surface_width) * tempor ;
                for(int k = 0  ; k < no_pha ; k++)
                    upd2[k][i * params.N2 + j] = 1./(1.+dt * Di[k] * SQR(surface_width) * tempor );
            }
        }

    }

    MultiPhaseKR_Update(const MultiPhaseKR_Update &a) : updateRules(a.par), Di(a.Di), Chi(a.Chi)
    {
        no_pha = a.no_pha;
        dt = a.dt;
        temp1 = a.temp1;
        A = a.A;
        surface_width = a.surface_width;

        upd1 = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));
        upd2 = new complex<double> *[no_pha];
        for (int j = 0; j < no_pha; j++)
            upd2[j] = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                for (int k = 0; k < no_pha; k++)
                    upd2[k][i * par.N2 + j] = a.upd2[k][i * par.N2 + j];
            }
        }

    }
    MultiPhaseKR_Update& operator=(const MultiPhaseKR_Update &a)
    {
        fftw_free(upd1);
        for (int i = 0; i < no_pha; i++)
        {
            fftw_free(upd2[i]);
        }
        delete upd2;
        par = a.par;
        Di = a.Di;
        Chi = a.Chi;
        no_pha = a.no_pha;
        surface_width = a.surface_width;
        dt = a.dt;
        temp1 = a.temp1;

        upd1 = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));
        upd2 = new complex<double> *[no_pha];
        for (int j = 0; j < no_pha; j++)
            upd2[j] = (complex<double> *)fftw_malloc(par.N1 * par.N2 * sizeof(complex<double>));

        for (int i = 0; i < par.N1; i++)
        {
            for (int j = 0; j < par.N2; j++)
            {
                upd1[i * par.N2 + j] = a.upd1[i * par.N2 + j];
                for (int k = 0; k < no_pha; k++)
                    upd2[k][i * par.N2 + j] = a.upd2[k][i * par.N2 + j];
            }
        }

        return *this;
    }

    ~MultiPhaseKR_Update()
    {
        fftw_free(upd1);
        for (int i = 0; i < no_pha; i++)
        {
            fftw_free(upd2[i]);
        }
        delete upd2;
        //delete upd2;
    }

    void operator()(complex<double> **res, complex<double> **f, complex<double> **w, complex<double> **r, int j, const CH_builder &p)
    {
        if(j>0) {

        }
        else{

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {

            for(int k = 0  ; k < no_pha ; k++) {

                    // complex<double> evo1 = upd2[k][i] * upd1[i] * (Di[k]) * f[k][i];
                    // complex<double> evo2 = upd2[k][i]*dt*w[k][i];
                    // complex<double> evo3 = dt * upd2[k][i] * r[k][i];

                    res[k][i] = upd2[k][i] * f[k][i] + upd2[k][i] * upd1[i] * (Di[k]) * f[k][i] + upd2[k][i]*dt*w[k][i] + dt * upd2[k][i] * r[k][i];

 
                }   
            }

        }
    }

    MultiPhaseKR_Update *clone() const
    {
        return new MultiPhaseKR_Update(*this);
    }
    void print()
    {
        cout << "MultiPhaseKRUpdate" << endl;
    }
};

//MultiPhaseKR(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, matrix<double> surface_widthh, matrix<double> Chii)
//MultiPhaseKR_Update(const CH_builder &params, int no_phaa, double dtt, vector1<double> Dii, double temp11, double surface_widthh, matrix<double> Chii)
  
struct MultiPhaseBundleReal {
    CoupledPhaseSeparatingSystem<double> A; 
    MultiPhaseKR B;
    MultiPhaseKR_Update C;

    MultiPhaseBundleReal(CoupledPhaseSeparatingSystem<double> &A2, const CH_builder &params, double dtt, vector1<double> Dii, double temp11, matrix<double> &surface_widthh)
        : A(A2),
          B(params, A2.get_number_phase(), dtt, Dii, temp11, surface_widthh, A2.getchi()),
          C(params, A2.get_number_phase(), dtt, Dii, temp11, 1., A2.getchi()) {

          }
};

struct MultiPhaseBundleComplex
{
    CoupledPhaseSeparatingSystem<complex<double> > A;
    MultiPhaseKR B;
    MultiPhaseKR_Update C;

    MultiPhaseBundleComplex(CoupledPhaseSeparatingSystem<complex<double> > &A2, const CH_builder &params, double dtt, vector1<double> Dii, double temp11, matrix<double> &surface_widthh, double A00)
        : A(A2),
          B(params, A2.get_number_phase(), dtt, Dii, temp11, surface_widthh, A2.getchi()),
          C(params, A2.get_number_phase(), dtt, Dii, temp11, A00, A2.getchi())
    {
    }
};

#endif /* UPDATERULES_CPP */
