 #ifndef FIELD_WRAPPER_CPP
#define FIELD_WRAPPER_CPP



//the template class here must be of the correct form for fftw_malloc to work on it
template <class T, class Q>
Field_Wrapper<T, Q>::Field_Wrapper()
{
    params.number_of_fields = 1;
    params.N1 = 128;
    params.N2 = 128;
    calculated_reactions = new T *[params.number_of_fields];

    calculated_reactions[0] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
    for (int j = 0; j < params.N1 * params.N2; j++)
    {
        T b;
        calculated_reactions[0][j] = b; //uninitialized, you need to initiliaze the field in order to set the values
    }

    chems = new Weight<T,Q>*[params.number_of_fields];

    NoWeight<T,Q> *chem1 = new NoWeight<T,Q>;
    chems[0] = chem1->clone();
    delete chem1;
}

template <class T, class Q>
Field_Wrapper<T, Q>::Field_Wrapper(const CH_builder &a)
{
    params.number_of_fields = a.number_of_fields;
    params.N1 = a.N1;
    params.N2 = a.N2;
    calculated_reactions = new T *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            T b;
            calculated_reactions[i][j] = b;
        }
    }

    chems = new Weight<T, Q> *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        NoWeight<T, Q> *chem1 = new NoWeight<T, Q>;
        chems[i] = chem1->clone();
        delete chem1;
    }
}


template <class T, class Q>
Field_Wrapper<T, Q>::Field_Wrapper(const Field_Wrapper<T,Q> &a)
{
    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;
    calculated_reactions = new T *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = a.calculated_reactions[i][j];
        }
    }

    chems = new Weight<T, Q> *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        Weight<T, Q> *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }
}

template <class T, class Q>
Field_Wrapper<T, Q> &Field_Wrapper<T,Q>::operator=(const Field_Wrapper<T,Q> &a)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        delete chems[i];
        fftw_free(calculated_reactions[i]);
    }
    delete[] chems;
    delete[] calculated_reactions;

    params.number_of_fields = a.params.number_of_fields;
    params.N1 = a.params.N1;
    params.N2 = a.params.N2;

    calculated_reactions = new T *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = a.calculated_reactions[i][j];
        }
    }

    chems = new Weight<T, Q> *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        Weight<T, Q> *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }

    return *this;
}

template <class T, class Q>
Field_Wrapper<T, Q>::~Field_Wrapper()
{
   // cout << "field wrapper free" << endl;
    
    for (int i = 0; i < params.number_of_fields; i++)
    {
        delete chems[i];
        fftw_free(calculated_reactions[i]);
    }
    delete[] chems;
    delete[] calculated_reactions;
}

template <class T, class Q>
void Field_Wrapper<T, Q>::add_method(Weight<T, Q> &a, int j)
{
    // Weight<T,Q> *chem1 = a.clone();
    delete chems[j];
    chems[j] = a.clone();
    // Weight<T, Q> *chem1 = a.clone();
    // chems[j] = chem1;
}


template <class T, class Q>
void Field_Wrapper<T, Q>::Calculate_Results(Q **fields)
{

    for (int i = 0; i < params.number_of_fields; i++)
        {
        chems[i]->operator()(calculated_reactions, fields, i, params);
        }



}

template <class T, class Q>
void Field_Wrapper<T,Q>::Add_Noise(GenNoise<T> &a) {
    for (int i = 0; i < params.number_of_fields; i++)
    {
        int end = params.get_total();
        for (int j = 0; j < end; j++)
            calculated_reactions[i][j] += a.random_field[i][j];
    }
}

template <class T, class Q>
void Field_Wrapper<T, Q>::Check_fields()
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
                exit(1);
            }
        }
    }
}

template <class T, class Q>
void Field_Wrapper<T, Q>::means(vector1<T> &vv)
{

    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {

            a += calculated_reactions[i][j];
            
        }
        vv[i] = a/T(tot);
    }
}

template <class T, class Q>
void Field_Wrapper<T, Q>::means() {
    vector1<T> vv(params.number_of_fields);
    this->means(vv);
    cout << vv << endl;
}

template <class T, class Q>
void Field_Wrapper<T, Q>::rescale()
{
    vector1<T> means(params.number_of_fields);
    vector1<T> mins(params.number_of_fields);
    vector1<T> maxs(params.number_of_fields);

    this->means(means);
    this->GetMinimas(mins);
    this->GetMaximas(maxs);


    T fac1 = T(2.);
    T fac2 = T(0.1);



    for (int i = 0; i < params.number_of_fields; i++)
    {
        T max = fac1 * means[i];
        T min = fac2 * means[i];
        int tot = params.get_total();
        for (int j = 0; j < tot; j++)
        {
            calculated_reactions[i][j] = min+(max-min)*(calculated_reactions[i][j]-mins[i])/(maxs[i]-mins[i]);
        }
    }
}

template <class T, class Q>
void Field_Wrapper<T, Q>::GetMaximas() {
    for(int i = 0  ; i < params.number_of_fields ; i++) {
        T a = calculated_reactions[i][0];
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
template <class T, class Q>
void Field_Wrapper<T, Q>::GetMaximas(vector1<T> &vv)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] > a)
            {
                a = calculated_reactions[i][j];
            }
        }
        vv[i] = a;
    }
}

template <class T>
void GetMaximas(T **calculated_reactions, const CH_builder &params)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int tot = params.get_total();
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

template <class T>
void GetMinimas(T **calculated_reactions, const CH_builder &params)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
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


template <class T, class Q>
void Field_Wrapper<T,Q>::GetMaximasIndex()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int ind = 0;
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] > a)
            {
                a = calculated_reactions[i][j];
                ind = j;
            }
        }
        int myi = floor((ind/params.N2));
        int myj = ind % params.N2;

        cout << myi <<" " << myj << "\t";
    }
    cout << endl;
}

template <class T, class Q>
void Field_Wrapper<T, Q>::GetMinimas()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
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

template <class T, class Q>
void Field_Wrapper<T, Q>::GetMinimas(vector1<T> &vv)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] < a)
            {
                a = calculated_reactions[i][j];
            }
        }
        vv[i] = a;
    }
}

template <class T, class Q>
void Field_Wrapper<T,Q>::GetMinimasIndex()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int ind = 0;
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {
            if (calculated_reactions[i][j] < a)
            {
                a = calculated_reactions[i][j];
                ind = j;
            }
        }
        int myi = floor((ind / params.N2));
        int myj = ind % params.N2;

        cout << myi << " " << myj << "\t";
    }
    cout << endl;
}

template <class T>
void GetTotals(T **calculated_reactions, const CH_builder &params)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        T a = calculated_reactions[i][0];
        int tot = params.get_total();
        for (int j = 1; j < tot; j++)
        {

                a += calculated_reactions[i][j];
            
        }
        cout << a << " ";
    }
    cout << endl;
}

#endif /* FIELD_WRAPPER_CPP */
