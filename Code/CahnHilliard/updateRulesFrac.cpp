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
    { // f is cn-1
    // w is cn
    // r is the initial condition

        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            // NB this assumes that the initial field has a zero time derivative.
            res[j][i] = -(0.5) * alpha * (alpha - 1.) * (-2. * f[j][i] + w[j][i]) - 0.5 * dt * alpha * (2 - alpha) * (-w[j][i]) - 0.5 * dt * dt * (alpha - 1) * (alpha - 2) * w[j][i] + dt * dt * 0.5 * (alpha - 1) * (alpha - 2) * r[j][i];

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

template <class T>
struct FractionalDtDiffusionWithInteractionI : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    // double *upd2;
    double temp1;
    double alpha;

    FractionalDtDiffusionWithInteractionI(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        alpha = alphaa;
        // upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
            }
        }
    }

    FractionalDtDiffusionWithInteractionI(const FractionalDtDiffusionWithInteractionI &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        alpha = a.alpha;
        // upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                // upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }
    }

    FractionalDtDiffusionWithInteractionI &operator=(const FractionalDtDiffusionWithInteractionI &a)
    {
        fftw_free(upd1);
        // fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        // upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                // upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }

        return *this;
    }

    ~FractionalDtDiffusionWithInteractionI()
    {
        fftw_free(upd1);
        // fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) // f is cn, w is wn
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = -dt*upd1[i]*w[j][i] - dt*upd1[i]*f[j][i];//-dt * upd1[i] * f[j][i] - dt * upd1[i] * w[j][i] - SQR(dt) * upd1[i] * ((1 - alpha) / alpha) * w[j][i];
        }
    }

    FractionalDtDiffusionWithInteractionI *clone() const
    {
        return new FractionalDtDiffusionWithInteractionI(*this);
    }
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithInteraction1 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    //double *upd2;
    double temp1;
    double alpha;

    FractionalDtDiffusionWithInteraction1(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        alpha = alphaa;
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));

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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
            }
        }
    }

    FractionalDtDiffusionWithInteraction1(const FractionalDtDiffusionWithInteraction1 &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        alpha = a.alpha;
        //upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                //upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }
    }

    FractionalDtDiffusionWithInteraction1 &operator=(const FractionalDtDiffusionWithInteraction1 &a)
    {
        fftw_free(upd1);
        //fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        //upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                //upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }

        return *this;
    }

    ~FractionalDtDiffusionWithInteraction1()
    {
        fftw_free(upd1);
        //fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) //f is cn, w is wn
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = /*-dt * upd1[i] * f[j][i]*/ - dt * upd1[i] * w[j][i] -SQR(dt)*upd1[i]*((1-alpha)/alpha)*w[j][i];
        }
    }

    FractionalDtDiffusionWithInteraction1 *clone() const
    {
        return new FractionalDtDiffusionWithInteraction1(*this);
    }
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithInteraction2 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    // double *upd2;
    double temp1;
    double alpha;

    FractionalDtDiffusionWithInteraction2(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        upd1 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        //upd2 = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
            }
        }
    }

    FractionalDtDiffusionWithInteraction2(const FractionalDtDiffusionWithInteraction2 &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        // upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                //upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }
    }

    FractionalDtDiffusionWithInteraction2 &operator=(const FractionalDtDiffusionWithInteraction2 &a)
    {
        fftw_free(upd1);
        //fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        upd1 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        //upd2 = (double *)fftw_malloc(a.par.N1 * a.par.N2 * sizeof(double));
        for (int i = 0; i < a.par.N1; i++)
        {
            for (int j = 0; j < a.par.N2; j++)
            {
                upd1[i * a.par.N2 + j] = a.upd1[i * a.par.N2 + j];
                //upd2[i * a.par.N2 + j] = a.upd2[i * a.par.N2 + j];
            }
        }

        return *this;
    }

    ~FractionalDtDiffusionWithInteraction2()
    {
        fftw_free(upd1);
        //fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) // f is cn-1, w is cn, r is wn-1
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = -f[j][i]+2.*w[j][i]+(dt/2)*((1-alpha)/(alpha))*f[j][i]+ dt*upd1[i]*(r[j][i]+0.5*f[j][i]) ;// -dt * upd1[i] * f[j][i] - dt * upd1[i] * w[j][i] - SQR(dt) * upd1[i] * ((1 - alpha) / alpha) * w[j][i];
        }
    }

    FractionalDtDiffusionWithInteraction2 *clone() const
    {
        return new FractionalDtDiffusionWithInteraction2(*this);
    }
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithInteraction3 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double alpha;

    FractionalDtDiffusionWithInteraction3(const CH_builder &params, double dtt, double Dii, double temp11, double alphaa) : updateRules<T, T, T, T>(params)
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (1. + (dt / 2) * (1 - alpha) / (alpha) + (dt/2)*Di * temp1 * tempor + SQR(dt) * (1 - alpha) / (alpha)*Di * temp1 * tempor);
            }
        }
    }

    FractionalDtDiffusionWithInteraction3(const FractionalDtDiffusionWithInteraction3 &a) : updateRules<T, T, T, T>(a.par)
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

    FractionalDtDiffusionWithInteraction3 &operator=(const FractionalDtDiffusionWithInteraction3 &a)
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

    ~FractionalDtDiffusionWithInteraction3()
    {
        fftw_free(upd1);
        //fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    { //f is the zero point
        // w is the combined weight of previous results
        // r is the chemistry
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = /* -upd2[i]*SQR(dt)*(1-alpha)/(alpha) *f[j][i] */ + upd2[i] * w[j][i] +SQR(dt)/alpha * upd2[i] * r[j][i];
        }
    }

    FractionalDtDiffusionWithInteraction3 *clone() const
    {
        return new FractionalDtDiffusionWithInteraction3(*this);
    }
    void print()
    {
        cout << "Diffusion With Interaction" << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithSurfaceTensionI : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double alpha;
    double eps;

    FractionalDtDiffusionWithSurfaceTensionI(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                upd2[i * params.N2 + j] = Di * SQR(eps) * SQR(temp1) * SQR(tempor);
            }
        }
    }

    FractionalDtDiffusionWithSurfaceTensionI(const FractionalDtDiffusionWithSurfaceTensionI &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    FractionalDtDiffusionWithSurfaceTensionI &operator=(const FractionalDtDiffusionWithSurfaceTensionI &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    ~FractionalDtDiffusionWithSurfaceTensionI()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) // f is cn, w is wn
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = -dt*upd2[i]*f[j][i]-dt*upd1[i]*w[j][i];//-dt * upd2[i] * f[j][i] - dt * upd1[i] * w[j][i] - SQR(dt) * upd1[i] * ((1 - alpha) / alpha) * w[j][i];
        }
    }

    FractionalDtDiffusionWithSurfaceTensionI *clone() const
    {
        return new FractionalDtDiffusionWithSurfaceTensionI(*this);
    }
    void print()
    {
        cout << "Surface tension I: " << eps << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithSurfaceTension1 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double alpha;
    double eps;

    FractionalDtDiffusionWithSurfaceTension1(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                upd2[i * params.N2 + j] = Di * SQR(eps) * SQR(temp1) * SQR(tempor);
            }
        }
    }

    FractionalDtDiffusionWithSurfaceTension1(const FractionalDtDiffusionWithSurfaceTension1 &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    FractionalDtDiffusionWithSurfaceTension1 &operator=(const FractionalDtDiffusionWithSurfaceTension1 &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    ~FractionalDtDiffusionWithSurfaceTension1()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) // f is cn, w is wn
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] =/*  -dt * upd2[i] * f[j][i] */ - dt * upd1[i] * w[j][i] - SQR(dt) * upd1[i] * ((1 - alpha) / alpha) * w[j][i];
        }
    }

    FractionalDtDiffusionWithSurfaceTension1 *clone() const
    {
        return new FractionalDtDiffusionWithSurfaceTension1(*this);
    }
    void print()
    {
        cout << "Surface tension 1: " << eps << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithSurfaceTension2 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double alpha;
    double eps;

    FractionalDtDiffusionWithSurfaceTension2(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                upd2[i * params.N2 + j] = Di * SQR(eps) * SQR(temp1) * SQR(tempor);
            }
        }
    }

    FractionalDtDiffusionWithSurfaceTension2(const FractionalDtDiffusionWithSurfaceTension2 &a) : updateRules<T, T, T, T>(a.par)
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

    FractionalDtDiffusionWithSurfaceTension2 &operator=(const FractionalDtDiffusionWithSurfaceTension2 &a)
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

    ~FractionalDtDiffusionWithSurfaceTension2()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p) // f is cn-1, w is cn, r is wn-1
    {
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = -f[j][i] + 2. * w[j][i] + (dt / 2) * ((1 - alpha) / (alpha)) * f[j][i] + dt * upd1[i] * (r[j][i]) + dt * 0.5 * upd2[i] *(f[j][i]); // -dt * upd1[i] * f[j][i] - dt * upd1[i] * w[j][i] - SQR(dt) * upd1[i] * ((1 - alpha) / alpha) * w[j][i];
        }
    }

    FractionalDtDiffusionWithSurfaceTension2 *clone() const
    {
        return new FractionalDtDiffusionWithSurfaceTension2(*this);
    }
    void print()
    {
        cout << "Surface tension 2: " << eps << endl;
    }
};

template <class T>
struct FractionalDtDiffusionWithSurfaceTension3 : updateRules<T, T, T, T>
{
    // double **calculated;

    double dt;
    double Di;
    double *upd1;
    double *upd2;
    double temp1;
    double alpha;
    double eps;

    FractionalDtDiffusionWithSurfaceTension3(const CH_builder &params, double dtt, double Dii, double temp11, double epss, double alphaa) : updateRules<T, T, T, T>(params)
    {
        Di = Dii;
        dt = dtt;
        temp1 = temp11;
        alpha = alphaa;
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

                upd1[i * params.N2 + j] = Di * temp1 * tempor;
                upd2[i * params.N2 + j] = 1. / (1. + (dt / 2) * (1 - alpha) / (alpha) +0.5*dt *Di * SQR(eps) * SQR(temp1) * SQR(tempor) + SQR(dt) * (1 - alpha) / (alpha)*Di * SQR(eps) * SQR(temp1) * SQR(tempor));
            }
        }
    }

    FractionalDtDiffusionWithSurfaceTension3(const FractionalDtDiffusionWithSurfaceTension3 &a) : updateRules<T, T, T, T>(a.par)
    {
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    FractionalDtDiffusionWithSurfaceTension3 &operator=(const FractionalDtDiffusionWithSurfaceTension3 &a)
    {
        fftw_free(upd1);
        fftw_free(upd2);
        this->setpar(a.par);
        Di = a.Di;
        dt = a.dt;
        temp1 = a.temp1;
        alpha = a.alpha;
        eps = a.eps;
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

    ~FractionalDtDiffusionWithSurfaceTension3()
    {
        fftw_free(upd1);
        fftw_free(upd2);
    }

    void operator()(T **res, T **f, T **w, T **r, int j, const CH_builder &p)
    { // f is the zero point
        // w is the combined weight of previous results
        // r is the chemistry
        int end = p.get_total();
        for (int i = 0; i < end; i++)
        {
            res[j][i] = /*-upd2[i]*SQR(dt) * (1 - alpha) / (alpha)*f[j][i]*/ + upd2[i] * w[j][i] + SQR(dt) / alpha * upd2[i] * r[j][i];
        }
    }

    FractionalDtDiffusionWithSurfaceTension3 *clone() const
    {
        return new FractionalDtDiffusionWithSurfaceTension3(*this);
    }
    void print()
    {
        cout << "Surface tension 3: " << eps << endl;
    }
};

struct FracRuleBundle {
    updateRules<complex<double>, complex<double>, complex<double>, complex<double> > *A;
    updateRules<complex<double>, complex<double>, complex<double>, complex<double> > *B;
    updateRules<complex<double>, complex<double>, complex<double>, complex<double>> *C;
    updateRules<complex<double>, complex<double>, complex<double>, complex<double>> *I;
};

struct IntBundle : FracRuleBundle {

    IntBundle(
            CH_builder &params,
            double dt,
              double Di,
              // double *upd2;
              double temp1,
              double alpha) {
        I = new FractionalDtDiffusionWithInteractionI<complex<double> >(params, dt, Di, temp1, alpha);
        A = new FractionalDtDiffusionWithInteraction1<complex<double> >(params, dt, Di, temp1, alpha);
        B = new FractionalDtDiffusionWithInteraction2<complex<double> >(params, dt, Di, temp1, alpha);
        C = new FractionalDtDiffusionWithInteraction3<complex<double> >(params, dt, Di, temp1, alpha);
              }
    ~IntBundle() {
        delete I;
        delete A;
        delete B;
        delete C;
    }
};

struct SurfBundle : FracRuleBundle
{

    SurfBundle(
        CH_builder &params,
        double dt,
        double Di,
        // double *upd2;
        double temp1,
        double eps,
        double alpha) 
    {
        I = new FractionalDtDiffusionWithSurfaceTensionI<complex<double> >(params, dt, Di, temp1, eps, alpha);
        A = new FractionalDtDiffusionWithSurfaceTension1<complex<double> >(params, dt, Di, temp1, eps, alpha);
        B = new FractionalDtDiffusionWithSurfaceTension2<complex<double> >(params, dt, Di, temp1, eps, alpha);
        C = new FractionalDtDiffusionWithSurfaceTension3<complex<double> >(params, dt, Di, temp1, eps, alpha);
        cout << "created bundle" << endl;
    }
    ~SurfBundle()
    {
        delete I;
        delete A;
        delete B;
        delete C;
    }
};

#endif