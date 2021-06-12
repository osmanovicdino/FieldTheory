#ifndef RULE_WRAPPER_CPP
#define RULE_WRAPPER_CPP

Rule_Wrapper::Rule_Wrapper()
{
    params.number_of_fields = 1;
    params.N1 = 128;
    params.N2 = 128;
    calculated_reactions = new double *[params.number_of_fields];

    calculated_reactions[0] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
    for (int j = 0; j < params.N1 * params.N2; j++)
    {
        calculated_reactions[0][j] = 0.0;
    }

    chems = new updateRules *[params.number_of_fields];

    NormalDiffusion *chem1 = new NormalDiffusion(params, 1.0, 1.0, 1.0);
    chems[0] = chem1;
}

Rule_Wrapper::Rule_Wrapper(const CH_builder &a)
{
    params.number_of_fields = a.number_of_fields;
    params.N1 = a.N1;
    params.N2 = a.N2;
    calculated_reactions = new double *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = 0.0;
        }
    }

    chems = new updateRules *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        NormalDiffusion *chem1 = new NormalDiffusion(params,1.0,1.0,1.0);
        chems[i] = chem1;
    }
}

Rule_Wrapper::Rule_Wrapper(const Rule_Wrapper &a)
{
    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;
    calculated_reactions = new double *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = a.calculated_reactions[i][j];
        }
    }

    chems = new updateRules *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        updateRules *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }
}

Rule_Wrapper &Rule_Wrapper::operator=(const Rule_Wrapper &a)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        delete chems[i];
        fftw_free(calculated_reactions[i]);
    }
    delete chems;
    delete calculated_reactions;

    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;

    calculated_reactions = new double *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (double *)fftw_malloc(params.N1 * params.N2 * sizeof(double));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = a.calculated_reactions[i][j];
        }
    }

    chems = new updateRules *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        updateRules *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }

    return *this;
}

Rule_Wrapper::~Rule_Wrapper()
{
    //cout << "called rule wrapper destructor" << endl;
    for (int i = 0; i < params.number_of_fields; i++)
    {
        delete chems[i];
        fftw_free(calculated_reactions[i]);
    }
    delete chems;
    delete calculated_reactions;
}

void Rule_Wrapper::add_method(updateRules &a, int j)
{

    //cout << "trying to add method" << endl;
    updateRules *chem1 = a.clone();

    //cout << "method cloned" << endl;
   
    chems[j] = chem1;

    //cout << "method set" << endl;
}

void Rule_Wrapper::Calculate_Results(double **fields, double **weis, double **reacs)
{
    #pragma omp parallel for
    for (int i = 0; i < params.number_of_fields; i++)
    {
        chems[i]->operator()(calculated_reactions, fields, weis, reacs, i, params);
    }
}

void Rule_Wrapper::Check_fields()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        int tot = params.N1 * params.N2;
        for (int j = 0; j < tot; j++)
        {
            if (calculated_reactions[i][j] != calculated_reactions[i][j])
            {
                cout << "NaN found" << endl;
                cout << i << " " << j << endl;
                pausel();
            }
        }
    }
}

#endif /* RULE_WRAPPER_CPP */
