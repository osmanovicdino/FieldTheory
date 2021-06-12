#ifndef FIELD_WRAPPER_CPP
#define FIELD_WRAPPER_CPP


Field_Wrapper::Field_Wrapper()
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

    chems = new Weight*[params.number_of_fields];

    NoWeight *chem1 = new NoWeight;
    chems[0] = chem1;
}


Field_Wrapper::Field_Wrapper(const CH_builder &a)
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

    chems = new Weight *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        NoWeight *chem1 = new NoWeight;
        chems[i] = chem1;
    }
}


Field_Wrapper::Field_Wrapper(const Field_Wrapper &a)
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

    chems = new Weight *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        Weight *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }
}


Field_Wrapper &Field_Wrapper::operator=(const Field_Wrapper &a)
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

    chems = new Weight *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        Weight *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }

    return *this;
}


Field_Wrapper::~Field_Wrapper()
{
   // cout << "field wrapper free" << endl;
    
    for (int i = 0; i < params.number_of_fields; i++)
    {
        delete chems[i];
        fftw_free(calculated_reactions[i]);
    }
    delete chems;
    delete calculated_reactions;
}


void Field_Wrapper::add_method(Weight &a, int j)
{
    //Chemistry *chem1 = a.clone();
    chems[j] = a.clone();
}


void Field_Wrapper::Calculate_Results(double **fields)
{

    #pragma omp parallel for
    for (int i = 0; i < params.number_of_fields; i++)
    {
        chems[i]->operator()(calculated_reactions, fields, i, params);
    }
}



void Field_Wrapper::Check_fields()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        int tot = params.get_total();
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

void Field_Wrapper::GetMaximas() {
    for(int i = 0  ; i < params.number_of_fields ; i++) {
        double a = calculated_reactions[i][0];
        int tot =  params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] > a)
            {
                a = calculated_reactions[i][j];
            }
        }
        cout << a << " ";
    }
    cout << endl;
 
}

void Field_Wrapper::GetMinimas()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        double a = calculated_reactions[i][0];
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] < a)
            {
                a = calculated_reactions[i][j];
            }
        }
        cout << a << " ";
    }
    cout << endl;
}

#endif /* FIELD_WRAPPER_CPP */
