#ifndef UPDATERULESFRAC_CPP
#define UPDATERULESFRAC_CPP

template <class T>
struct FractionalSubDiffusionWithInteraction : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double alpha;
    double *upd1;
    double *upd2;
    double temp1;

    FractionalSubDiffusionWithInteraction(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = dt * Di * temp1 * tempor/alpha;
                upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor/alpha);
            }
        }
    }

    FractionalSubDiffusionWithInteraction(const FractionalSubDiffusionWithInteraction &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
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

    FractionalSubDiffusionWithInteraction &operator=(const FractionalSubDiffusionWithInteraction &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
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

    ~FractionalSubDiffusionWithInteraction()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    { // f here must have the correct form for the correct implementation of the update rule

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + (dt/alpha) * upd2[i] * r[j][i];
        }
    }

    FractionalSubDiffusionWithInteraction *clone() const
    {
        return new FractionalSubDiffusionWithInteraction(*this);
    }
    void print()
    {
        cout << "FractionalWave Diffusion With Interaction" << endl;
    }
};


template <class T>
struct FractionalWaveDiffusionWithInteraction : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double alpha;
    double *upd1;
    double *upd2;
    double temp1;

    FractionalWaveDiffusionWithInteraction(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = dt * dt * Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (0.5 * (alpha) * (alpha - 1) + 0.5 * dt * alpha * (2 - alpha) + dt * dt * Di * temp1 * tempor);
            }
        }
    }

    FractionalWaveDiffusionWithInteraction(const FractionalWaveDiffusionWithInteraction &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
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

    FractionalWaveDiffusionWithInteraction &operator=(const FractionalWaveDiffusionWithInteraction &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
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

    ~FractionalWaveDiffusionWithInteraction()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    { // f here must have the correct form for the correct implementation of the update rule

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + SQR(dt) * upd2[i] * r[j][i];
        }
    }

    FractionalWaveDiffusionWithInteraction *clone() const
    {
        return new FractionalWaveDiffusionWithInteraction(*this);
    }
    void print()
    {
        cout << "FractionalWave Diffusion With Interaction" << endl;
    }
};

template <class T>
struct FractionalSubDiffusionWithSurfaceTension : updateRules<T, T, T, T>
{
    double dt;
    double Di;
    double alpha;
    double *upd1;
    double *upd2;
    double temp1;
    double eps;

    FractionalSubDiffusionWithSurfaceTension(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        eps = epss;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = (dt/alpha) * Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (1. + (dt/alpha) * Di * SQR(eps) * SQR(temp1) * SQR(tempor));
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

    FractionalSubDiffusionWithSurfaceTension(const FractionalSubDiffusionWithSurfaceTension &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        alpha = a.alpha;
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

    FractionalSubDiffusionWithSurfaceTension &operator=(const FractionalSubDiffusionWithSurfaceTension &a)
    {
        // cout << "called = " << endl;
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        alpha = a.alpha;
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

    ~FractionalSubDiffusionWithSurfaceTension()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + (dt / alpha) * upd2[i] * r[j][i];
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

    FractionalSubDiffusionWithSurfaceTension *clone() const
    {
        return new FractionalSubDiffusionWithSurfaceTension(*this);
    }
    void print()
    {
        cout << "Sub Diffusion With Surface Tension" << endl;
    }
};

template <class T>
struct FractionalWaveDiffusionWithSurfaceTension : updateRules<T, T, T, T>
{
    double dt;
    double Di;
    double alpha;
    double *upd1;
    double *upd2;
    double temp1;
    double eps;

    FractionalWaveDiffusionWithSurfaceTension(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        eps = epss;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = dt * dt * Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (0.5 * (alpha) * (alpha - 1) + 0.5 * dt * alpha * (2 - alpha) + SQR(dt) * Di * SQR(eps) * SQR(temp1) * SQR(tempor));
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

    FractionalWaveDiffusionWithSurfaceTension(const FractionalWaveDiffusionWithSurfaceTension &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        alpha = a.alpha;
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

    FractionalWaveDiffusionWithSurfaceTension &operator=(const FractionalWaveDiffusionWithSurfaceTension &a)
    {
        // cout << "called = " << endl;
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        eps = a.eps;
        alpha = a.alpha;
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

    ~FractionalWaveDiffusionWithSurfaceTension()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = upd2[i] * f[j][i] - upd2[i] * upd1[i] * w[j][i] + SQR(dt) * upd2[i] * r[j][i];
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

    FractionalWaveDiffusionWithSurfaceTension *clone() const
    {
        return new FractionalWaveDiffusionWithSurfaceTension(*this);
    }
    void print()
    {
        cout << "Diffusion With Surface Tension" << endl;
    }
};

template <class T>
struct FractionalSubDiffusionCalculateUpdateWeight : updateRules<T, T, T, T>
{
    double dt;
    double alpha;

    FractionalSubDiffusionCalculateUpdateWeight(const CH_builder &params, double dtt, double alphaa) : updateRules<T, T, T, T>(params)
    {
        dt = dtt;
        alpha = alphaa;
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

    FractionalSubDiffusionCalculateUpdateWeight(const FractionalSubDiffusionCalculateUpdateWeight &a) : updateRules<T, T, T, T>(a.par)
    {

        dt = a.dt;
        alpha = a.alpha;
    }

    FractionalSubDiffusionCalculateUpdateWeight &operator=(const FractionalSubDiffusionCalculateUpdateWeight &a)
    {
        // cout << "called = " << endl;

        this->setpar(a.par);
        dt = a.dt;
        alpha = a.alpha;

        return *this;
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            // w is the initial condition
            res[j][i] = f[j][i]-(1-alpha)*(dt/alpha)*(f[j][i]-w[j][i]);

            

            // cout << res[j][i] << endl;
            // cout << f[j][i] << endl;
            // cout << w[j][i] << endl;
            // cout << r[j][i] << endl;
            // pausel();
            //if(abs(res[j][i])>1000.) {
            //     cout << i << endl;
            //     cout << f[j][i] << endl;
            //     cout << w[j][i] << endl;
            //     cout << r[j][i] << endl;
            //     pausel();
            // }
        }
    }

    FractionalSubDiffusionCalculateUpdateWeight *clone() const
    {
        return new FractionalSubDiffusionCalculateUpdateWeight(*this);
    }
    void print()
    {
        cout << "sub diffusion weight" << endl;
    }
};

template <class T>
struct FractionalWaveDiffusionCalculateUpdateWeight : updateRules<T, T, T, T>
{
    double dt;
    double alpha;

    FractionalWaveDiffusionCalculateUpdateWeight(const CH_builder &params, double dtt, double alphaa) : updateRules<T, T, T, T>(params)
    {
        dt = dtt;
        alpha = alphaa;
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

    FractionalWaveDiffusionCalculateUpdateWeight(const FractionalWaveDiffusionCalculateUpdateWeight &a) : updateRules<T, T, T, T>(a.par)
    {

        dt = a.dt;
        alpha = a.alpha;
    }

    FractionalWaveDiffusionCalculateUpdateWeight &operator=(const FractionalWaveDiffusionCalculateUpdateWeight &a)
    {
        // cout << "called = " << endl;

        this->setpar(a.par);
        dt = a.dt;
        alpha = a.alpha;

        return *this;
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    {

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            // NB this assumes that the initial field has a zero time derivative.
            res[j][i] = -(0.5) * alpha * (alpha - 1.) * (-2. * f[j][i] + w[j][i]) - 0.5 * dt * alpha * (2 - alpha) * (-w[j][i]) - 0.5 * dt * dt * (alpha - 1) * (alpha - 2) * f[j][i] + dt * dt * 0.5 * (alpha - 1) * (alpha - 2) * r[j][i];

            // cout << res[j][i] << endl;
            // cout << f[j][i] << endl;
            // cout << w[j][i] << endl;
            // cout << r[j][i] << endl;
            // pausel();
            // if(abs(res[j][i])>1000.) {
            //     cout << i << endl;
            //     cout << f[j][i] << endl;
            //     cout << w[j][i] << endl;
            //     cout << r[j][i] << endl;
            //     pausel();
            // }
        }
    }

    FractionalWaveDiffusionCalculateUpdateWeight *clone() const
    {
        return new FractionalWaveDiffusionCalculateUpdateWeight(*this);
    }
    void print()
    {
        cout << "Diffusion With Surface Tension" << endl;
    }
};

#endif