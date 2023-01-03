#ifndef RULE_WRAPPER_CPP
#define RULE_WRAPPER_CPP

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3>::Rule_Wrapper()
{
    params.number_of_fields = 1;
    params.N1 = 128;
    params.N2 = 128;
    calculated_reactions = new T *[params.number_of_fields];

    calculated_reactions[0] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
    for (int j = 0; j < params.N1 * params.N2; j++)
    {
        T b = 0.0;
        calculated_reactions[0][j] = b;
    }

    chems = new updateRules<T,T1,T2,T3> *[params.number_of_fields];

    NoRule<T, T1, T2, T3> *chem1 = new NoRule<T, T1, T2, T3>(params);
    chems[0] = chem1;
}

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3>::Rule_Wrapper(const CH_builder &a)
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
            T b = 0.0;
            calculated_reactions[i][j] = b;
        }
    }
    chems = new updateRules<T, T1, T2, T3> *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        NoRule<T, T1, T2, T3> *chem1 = new NoRule<T, T1, T2, T3>(params);
        chems[i] = chem1;
    }
}

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3>::Rule_Wrapper(const Rule_Wrapper<T,T1,T2,T3> &a)
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

    chems = new updateRules<T, T1, T2, T3> *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        updateRules<T,T1,T2,T3> *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }
}

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3> &Rule_Wrapper<T, T1, T2, T3>::operator=(const Rule_Wrapper<T, T1, T2, T3> &a)
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

    calculated_reactions = new T *[params.number_of_fields];
    for (int i = 0; i < params.number_of_fields; i++)
    {
        calculated_reactions[i] = (T *)fftw_malloc(params.N1 * params.N2 * sizeof(T));
        for (int j = 0; j < params.N1 * params.N2; j++)
        {
            calculated_reactions[i][j] = a.calculated_reactions[i][j];
        }
    }


    chems = new updateRules<T, T1, T2, T3> *[params.number_of_fields];

    for (int i = 0; i < params.number_of_fields; i++)
    {
        updateRules<T, T1, T2, T3> *chem1 = a.chems[i]->clone();
        chems[i] = chem1;
    }

    return *this;
}

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3> &Rule_Wrapper<T, T1, T2, T3>::operator+=(const Rule_Wrapper<T, T1, T2, T3> &a) {
    int tot = params.N1 * params.N2;
    
        for (int i = 0; i < params.number_of_fields; i++)
        {
            for (int j = 0; j < tot; j++) {
                calculated_reactions[i][j] += a.calculated_reactions[i][j];
            }
        }

        return *this;
}

template <class T, class T1, class T2, class T3>
Rule_Wrapper<T, T1, T2, T3>::~Rule_Wrapper()
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

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::add_method(updateRules<T, T1, T2, T3> &a, int j)
{
    updateRules<T,T1,T2,T3> *chem1 = a.clone();
   
    chems[j] = chem1;

}

template <class T, class T1, class T2, class T3> //allow for the different fields to be of different types
void Rule_Wrapper<T, T1, T2, T3>::Calculate_Results(T1 **fields, T2 **weis, T3 **reacs)
{
    #pragma omp parallel for
    for (int i = 0; i < params.number_of_fields; i++)
    {
        //cout << "doing rule " << i << endl;
        chems[i]->operator()(calculated_reactions, fields, weis, reacs, i, params);
    }
}

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::Check_fields()
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
                exit(1);
            }
        }
    }
}

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::Check_negatives(bool &neg)
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        int tot = params.N1 * params.N2;
        for (int j = 0; j < tot; j++)
        {
            if (calculated_reactions[i][j] < 0 )
            {
                neg = true;
            }
        }
    }
}

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::GetMaximas()
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

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::GetMaximasIndex()
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
        int myi = floor((ind / params.N2));
        int myj = ind % params.N2;

        cout << myi << " " << myj << "\t";
    }
    cout << endl;
}

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::GetMinimas()
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

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::GetMinimasIndex()
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

template <class T, class T1, class T2, class T3>
void Rule_Wrapper<T, T1, T2, T3>::GetTotal()
{
    for (int i = 0; i < params.number_of_fields; i++)
    {
        int ind = 0;
        int tot = params.get_total();
        T ever = T(0.0);
        for (int j = 0; j < tot; j++)
        {

                ever += calculated_reactions[i][j];
 
            
        }


        cout << ever << "\t";
    }
    cout << endl;
}

#endif /* RULE_WRAPPER_CPP */
