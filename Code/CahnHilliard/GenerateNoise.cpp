#ifndef GENERATENOISE_CPP
#define GENERATENOISE_CPP
template <class T>
GenNoise<T>::GenNoise()
{
    params.number_of_fields = 1;
    params.N1 = 128;
    params.N2 = 128;
    random_field = new T *[params.number_of_fields];

    random_field[0] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
    for (int j = 0; j < params.N1 * params.N2; j++)
    {
        T b;
        random_field[0][j] = b; // uninitialized, you need to initiliaze the field in order to set the values
    }


}

template <class T>
GenNoise<T>::GenNoise(const CH_builder &a)
{
    params.number_of_fields = a.number_of_fields;
    params.N1 = a.N1;
    params.N2 = a.N2;
    random_field = new T *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        random_field[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            T b;
            random_field[i][j] = b;
        }
    }
}

template <class T>
GenNoise<T>::GenNoise(const GenNoise<T> &a)
{
    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;
    random_field = new T *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        random_field[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            random_field[i][j] = a.random_field[i][j];
        }
    }


}

template <class T>
GenNoise<T> &GenNoise<T>::operator=(const GenNoise<T> &a)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        fftw_free(random_field[i]);
    }

    delete random_field;

    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;

    random_field = new T *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        random_field[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            random_field[i][j] = a.random_field[i][j];
        }
    }


    return *this;
}

template <class T>
GenNoise<T>::~GenNoise()
{
    // cout << "field wrapper free" << endl;

    for (int i = 0; i < params.number_of_fields; i++)
    {

        fftw_free(random_field[i]);
    }
    delete random_field;
}

void create_k_values(double *val, int n, double samprate) {
    val[0] = 0.;
    double dn = (double)n;
    for(int i =  1 ; i < n/2 ; i++) {
        val[i] = (double)i/(dn*samprate);
    }
    val[n/2] = - 1/(2*samprate);
    for(int i =  1 ; i < n/2 ; i++) {
        val[n - i] = -(double)i / (dn * samprate);
    }
}

// typedef complex<double> cc;

void create_rands(complex<double> *val, int n) //ran is a random 
{
    double *ran = new double[(n/2)-1];
    for(int i = 0  ; i < n/2-1 ; i++) {
        ran[i] = 2 * pi * (double)rand() / (double)RAND_MAX;
    }
    val[0] = (2*(double)rand()/(double)RAND_MAX-1.0);
    for (int i = 1; i < n / 2; i++)
    {
        val[i] = complex<double>(cos(ran[i-1]),sin(ran[i-1]));
    }
    val[n / 2] = (2 * (double)rand() / (double)RAND_MAX - 1.0);
    for (int i = 1; i < n / 2; i++)
    {
        val[n - i] = complex<double>(cos(ran[i - 1]), -sin(ran[i - 1]));
    }
}

template <class T>
template <class Q>
void GenNoise<T>::GenFields(Q &func, double *strs, double samprat)
{

    if (params.N1 != params.N2)
        error("for now, generating noise requires square matrices");

    double *ks = new double [params.N1];

    create_k_values(ks, params.N1, samprat);
    int tot = params.N1 * params.N2;

    for (int i = 0; i < params.number_of_fields; i++)
    {
        //random_field[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1; j++) {
            for (int k = 0; k < params.N2; k++)
            {
            random_field[i][j*params.N1+k] = func(ks[j],ks[k]); //generate frequency data
            }
        }
}

complex<double> *rand1 = new complex<double>[params.N1];
complex<double> *rand2 = new complex<double>[params.N1];

for (int i = 0; i < params.number_of_fields; i++)
{
    //make the randoms
    create_rands(rand1, params.N1);
    create_rands(rand2, params.N1);
    for (int j = 0; j < params.N1; j++)
    {
        for (int k = 0; k < params.N2; k++)
        {
            random_field[i][j * params.N1 + k] = strs[i]*rand1[j] * rand2[k] * random_field[i][j * params.N1 + k]; // generate frequency data
        }
    }
}

delete ks;
delete rand1;
delete rand2;


}

#endif /* GENERATENOISE_CPP */
