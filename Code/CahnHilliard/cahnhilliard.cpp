#ifndef CAHNHILLIARD_CPP
#define CAHNHILLIARD_CPP

template<class T>
CH<T>::CH(const CH_builder &p) : myp(p), chems(Field_Wrapper<T, T>(p)),
                              weigs(Field_Wrapper<T,T>(p)), transformed1(Field_Wrapper<T, T>(p)),
                              transformed2(Field_Wrapper<T,T>(p)), transformed3(Field_Wrapper<T, T>(p)),
                              rules(Rule_Wrapper<T, T,T, T>(p)), reverse_transform(Field_Wrapper<T, T>(p))
{

    fields = new T *[myp.number_of_fields];
    int totp =  myp.get_total();
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        fields[i] = (T *)fftw_malloc(totp * sizeof(T));
        for (int j = 0; j < totp; j++)
        {
            fields[i][j] = 0.0;
        }
    }

    cout << "dim: " << myp.dimension << endl;

    if(myp.dimension == 2) {
        for(int i = 0  ; i < myp.number_of_fields ; i++) {
            if constexpr (std::is_same_v<T, double>) {
                CosineWeightForward fw;
                transformed1.add_method(fw, i);
                transformed2.add_method(fw, i);
                transformed3.add_method(fw, i);
                CosineWeightBackward fw2;
                reverse_transform.add_method(fw2,i);
            }
            else if constexpr (std::is_same_v<T, complex<double> >) {
                FourierWeightForward2D fw;
                transformed1.add_method(fw, i);
                transformed2.add_method(fw, i);
                transformed3.add_method(fw, i);
                FourierWeightBackward2D fw2;
                reverse_transform.add_method(fw2, i);
            }
            else{
                cout << "will neeed to specify your own functional rules for custom types" << endl;
            }
        }
    }
    else{
        for (int i = 0; i < myp.number_of_fields; i++)
        {
            if constexpr (std::is_same_v<T, double>)
            {
                CosineWeightForward fw;
                transformed1.add_method(fw, i);
                transformed2.add_method(fw, i);
                transformed3.add_method(fw, i);
                CosineWeightBackward fw2;
                reverse_transform.add_method(fw2, i);
            }
            else if constexpr (std::is_same_v<T, complex<double>>) {
                cout << "current 3D rules only supported for periodic boundary conditions" << endl;
                FourierWeightForward3D fw;
                transformed1.add_method(fw, i);
                transformed2.add_method(fw, i);
                transformed3.add_method(fw, i);
                FourierWeightBackward3D fw2;
                reverse_transform.add_method(fw2, i);
            }
            else{
                
            }

        }
    }

}

template <class T>
CH<T>::~CH() {
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        fftw_free(fields[i]);
    }
    delete fields;

}

template <class T>
void CH<T>::set_field(const matrix<T> &mymat, int k)
{
    for(int i =  0 ;  i < myp.N1 ; i++) {
        for(int j = 0  ; j < myp.N2 ; j++) {
            fields[k][i*myp.N2+j] =  mymat.gpcons(i,j);
        }
    }
}


template<class T>
void CH<T>::set_field(T **orig) { // only set the field to real values
    if(myp.dimension ==2) {
        for(int k = 0 ; k < myp.number_of_fields ; k++) {
            for (int i = 0; i < myp.N1; i++)
            {
                for (int j = 0; j < myp.N2; j++)
                {
                    //fields[k][i * myp.N2 + j] = {orig[k][i * myp.N2 + j].real(), 0.};
                    fields[k][i * myp.N2 + j] = T(orig[k][i * myp.N2 + j]);
                }
            }
        }
    }
    else{
        for (int k = 0; k < myp.number_of_fields; k++)
        {
            for (int i = 0; i < myp.N1; i++)
            {
                for (int j = 0; j < myp.N2; j++)
                {
                    for (int lk = 0; lk < myp.N3; lk++)
                    {
                        //fields[k][i * myp.N2 * myp.N3 + j * myp.N3 + lk] = {orig[k][i * myp.N2 * myp.N3 + j * myp.N3 + lk].real(), 0.};
                        fields[k][i * myp.N2 + j] = T(orig[k][i * myp.N2 + j]);
                    }
                }
            }
        }
    }
}


template <>
void CH<complex<double>>::set_field(complex<double> *orig, int k)
{ // only set the field to real values

int totp = myp.get_total();

    for(int i  = 0 ; i < totp ; i++)            
        fields[k][i] = orig[i].real();
            
        
    
}

template <class T>
void CH<T>::set_field(T *orig, int k)
{ // only set the field to real values

    int totp = myp.get_total();

    for (int i = 0; i < totp; i++)
        fields[k][i] = T(orig[i]);
}

template <class T>
void CH<T>::high_k_correct(int cut_offf) {
    transformed1.Calculate_Results(fields); // calculate FT of fields
    int totp = myp.get_total();
    int nof = myp.number_of_fields;

    // cut off form
    int cut_off = cut_offf;
    if(myp.dimension==2)
    for (int fn = 0; fn < nof; fn++)
    {
        for (int i = 0; i < myp.N1; i++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                int k1, k2;
                if (i <= myp.N1 / 2)
                {
                    k1 = i;
                }
                else
                {
                    k1 = (i - myp.N1);
                }
                if (j <= myp.N2 / 2)
                {
                    k2 = j;
                }
                else
                {
                    k2 = (j - myp.N2);
                }
                double tempor = SQR(k1) + SQR(k2);
                if (tempor > cut_off)
                {
                    transformed1.calculated_reactions[fn][i * myp.N2 + j] = 0.; // cut off the high frequency modes
                }
            }
        }
    }
    else if(myp.dimension==3) {
        for (int fn = 0; fn < nof; fn++)
        {
            for (int i = 0; i < myp.N1; i++)
            {
                for (int j = 0; j < myp.N2; j++)
                {
                    for(int k = 0 ; k < myp.N3 ; k++) {
                        int k1, k2, k3;
                        if (i <= myp.N1 / 2)
                        {
                            k1 = i;
                        }
                        else
                        {
                            k1 = (i - myp.N1);
                        }
                        if (j <= myp.N2 / 2)
                        {
                            k2 = j;
                        }
                        else
                        {
                            k2 = (j - myp.N2);
                        }
                        if (k <= myp.N3 / 2)
                        {
                            k3 = k;
                        }
                        else
                        {
                            k3 = (k - myp.N2);
                        }
                        double tempor = SQR(k1) + SQR(k2) + SQR(k3);
                        if (tempor > cut_off)
                        {
                            transformed1.calculated_reactions[fn][i * myp.N2 * myp.N3 + j * myp.N3 + k] = 0.; // cut off the high frequency modes
                        }
                    }
                }
            }
        }
    }


    reverse_transform.Calculate_Results(transformed1.calculated_reactions);

    set_field(reverse_transform.calculated_reactions);
}


template<class T>
void CH<T>::Update() {
    //cout << fields[0][0] << endl;



    cout << "start" << endl;

    GetMaximas(fields,myp);
    GetMinimas(fields,myp);
    cout << endl;

    weigs.Calculate_Results(fields);

    cout << "weigs done" << endl;
    // weigs.GetMaximas();
    // weigs.GetMinimas();

    chems.Calculate_Results(fields);

    cout << "chems done" << endl;
    // chems.GetMaximas();
    // chems.GetMinimas();

    transformed1.Calculate_Results(fields);

    cout << "FFT1 done" << endl;
    // transformed1.GetMaximas();
    // transformed1.GetMinimas();
    // cout << 3 << endl;

    transformed2.Calculate_Results(weigs.calculated_reactions);

    cout << "FFT2 done" << endl;
    // transformed2.GetMaximas();
    // transformed2.GetMinimas();
    // cout << 4 << endl;

    transformed3.Calculate_Results(chems.calculated_reactions);

    cout << "FFT3 done" << endl;
    // transformed3.GetMaximas();
    // transformed3.GetMinimas();
    // cout << 5 << endl;

    // cout << transformed1.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[1][0] << endl;
    // cout << transformed1.calculated_reactions[2][0] << endl;

    // cout << endl;
    rules.Calculate_Results(transformed1.calculated_reactions,transformed2.calculated_reactions,transformed3.calculated_reactions);


    rules.Check_fields();
    cout << "Rules done" << endl;
    //rules.GetMaximas();
    // cout << transformed1.calculated_reactions[0][0] << endl;
    // cout << rules.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[0][1] << endl;
    // cout << rules.calculated_reactions[0][1] << endl;

    // cout << transformed1.calculated_reactions[0][2] << endl;
    // cout << rules.calculated_reactions[0][2] << endl;
    string filename1 = "xfield";
    string filename2 = "xrealspaceweigh";
    string filename3 = "xfieldtransform";
    string filename4 = "xweighttransform";
    string filename5 = "xnewfield";
    string filename6 = "xrealspacenewfield";



    // cout << endl;
    //  cout << 6 << endl;
    reverse_transform.Calculate_Results(rules.calculated_reactions);
    // cout << "beginning renorm" << endl;
    // for (int i = 0; i < myp.number_of_fields; i++)
    //     for (int j = 0; j < myp.get_total(); j++)
    //         reverse_transform.calculated_reactions[i][j] = {reverse_transform.calculated_reactions[i][j].real(), 0.};
    // cout << "renorm" << endl;
    // outfunc(fields[5], filename1, myp);
    // outfunc(weigs.calculated_reactions[5], filename2, myp);
    // outfunc(transformed1.calculated_reactions[5], filename3, myp);
    // outfunc(transformed2.calculated_reactions[5], filename4, myp);
    // outfunc(rules.calculated_reactions[5], filename5, myp);
    // outfunc(reverse_transform.calculated_reactions[5], filename6, myp);

    cout << "Reverse Transform done" << endl;
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;

     //pausel();
    //cout << reverse_transform.calculated_reactions[0][0] << endl;

    set_field(reverse_transform.calculated_reactions);

    // cout << 8 << endl;
}

template<class T>
void CH<T>::print_all_results(string s1) {
    for(int i = 0  ; i < myp.number_of_fields ; i++) {
        stringstream ss;
        ss << i;
        string fiel =  ss.str();
        string fie =  "field";
        string filename = fie + fiel + s1;
        outfunc(fields[i],filename,myp);
    }

    //outfunc(fields[0],);
}

template <class T>
void CH<T>::print_some_results(string s1,vector1<bool> &ps)
{
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        if(ps[i]) {
        stringstream ss;
        ss << i;
        string fiel = ss.str();
        string fie = "field";
        string filename = fie + fiel + s1;
        outfunc(fields[i], filename, myp);
        }
    }

    // outfunc(fields[0],);
}

void CHN::Update() {



    weigs.Calculate_Results(fields);

    // cout << "weigs done" << endl;
    // weigs.GetMaximas();
    // weigs.GetMinimas();
    // pausel();

    cout << endl;
    // cout << 1 << endl;

    chems.Calculate_Results(fields);

    // cout << 2 << endl;

    transformed1.Calculate_Results(fields);


    // cout << 3 << endl;
    // cout << "FFT1 done" << endl;
    // transformed1.GetMaximas();
    // transformed1.GetMinimas();
    // cout << 3 << endl;

    // cout << transformed1.calculated_reactions[0][0] << endl;

    transformed2.Calculate_Results(weigs.calculated_reactions);

    // cout << transformed2.calculated_reactions[0][0] << endl;
    // cout << "FFT2 done" << endl;
    // transformed2.GetMaximas();
    // transformed2.GetMinimas();
    // cout << 4 << endl;

    transformed3.Calculate_Results(chems.calculated_reactions);


    // cout << transformed3.calculated_reactions[0][0] << endl; //total change in weights
    // cout << "FFT3 done" << endl;
    // transformed3.GetMaximas();
    // transformed3.GetMinimas();
    // cout << 5 << endl;
    // string filename10 = "0field";
    // string filename20 = "0realspaceweigh";
    // string filename30 = "0fieldtransform";
    // string filename40 = "0weighttransform";
    // string filename50 = "0newweighttransform";

    // string filename1 = "5field";
    // string filename2 = "5realspaceweigh";
    // string filename3 = "5fieldtransform";
    // string filename4 = "5weighttransform";
    // string filename5 = "5newweighttransform";
    //string filename6 = "xrealspacenewfield";
    // cout << transformed1.calculated_reactions[0][0] << endl

    // outfunc(fields[5], filename1, myp);
    // outfunc(weigs.calculated_reactions[5], filename2, myp);
    // outfunc(transformed1.calculated_reactions[5], filename3, myp);
    // outfunc(transformed2.calculated_reactions[5], filename4, myp);

    // outfunc(fields[0], filename10, myp);

    // outfunc(weigs.calculated_reactions[0], filename20, myp);
    // outfunc(transformed1.calculated_reactions[0], filename30, myp);
    // outfunc(transformed2.calculated_reactions[0], filename40, myp);


    newrules.Calculate_Results(fields,transformed2.calculated_reactions, transformed1.calculated_reactions); // don't need chemistry here
    //new rules destroys the second one here


    // outfunc(newrules.calculated_reactions[5], filename5, myp);
    // outfunc(newrules.calculated_reactions[0], filename50, myp);
    // cout << newrules.calculated_reactions[0][0] << endl;
    // cout << 6 << endl;
    // cout << endl;
    rules.Calculate_Results(transformed1.calculated_reactions, newrules.calculated_reactions, transformed3.calculated_reactions);

    // cout << rules.calculated_reactions[0][0] << endl;
    // cout << 7 << endl;

    rules.Check_fields();
    // cout << "Rules done" << endl;
    //rules.GetMaximas();
    // cout << transformed1.calculated_reactions[0][0] << endl;
    // cout << rules.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[0][1] << endl;
    // cout << rules.calculated_reactions[0][1] << endl;

    // cout << transformed1.calculated_reactions[0][2] << endl;
    // cout << rules.calculated_reactions[0][2] << endl;


    // cout << endl;
    //  cout << 6 << endl;
    reverse_transform.Calculate_Results(rules.calculated_reactions);

    // outfunc(fields[5], filename1, myp);
    // outfunc(weigs.calculated_reactions[5], filename2, myp);
    // outfunc(transformed1.calculated_reactions[5], filename3, myp);
    // outfunc(transformed2.calculated_reactions[5], filename4, myp);
    // outfunc(rules.calculated_reactions[5], filename5, myp);
    // outfunc(reverse_transform.calculated_reactions[5], filename6, myp);
    cout << 8 << endl;

    
    cout << "Reverse Transform done" << endl;
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    // cout << endl;
    // cout << 7 << endl;

    // pausel();
    //pausel();
    //cout << reverse_transform.calculated_reactions[0][0] << endl;

    set_field(reverse_transform.calculated_reactions); 
}

template <class Q>
void CHWithNoise<Q>::Update() {

    weigs.Calculate_Results(fields);

    // cout << "weigs done" << endl;
    // weigs.GetMaximas();
    // weigs.GetMinimas();
    // string chem22 = "fields2";
    // outfunc(fields[0], chem22, myp);
    

    chems.Calculate_Results(fields);

    // string chem2= "chem2";
    // outfunc(chems.calculated_reactions[0],chem2,myp);
    // pausel();
    
    // cout << "chemical dynamics: ";
    // for (int i = 0; i < myp.number_of_fields; i++)
    //     cout << chems.calculated_reactions[i][0] << ",";
    // cout << endl;
    // cout << "chems done" << endl;
    // chems.GetMaximas();
    // chems.GetMinimas();

    transformed1.Calculate_Results(fields);

    // cout << "FFT1 done" << endl;
    // transformed1.GetMaximas();
    // transformed1.GetMinimas();
    // cout << 3 << endl;

    transformed2.Calculate_Results(weigs.calculated_reactions);


    // cout << "FFT2 done" << endl;
    // transformed2.GetMaximas();
    // transformed2.GetMinimas();
    // cout << 4 << endl;

    transformed3.Calculate_Results(chems.calculated_reactions);

    // cout << "FFT3 done" << endl;
    // transformed3.GetMaximas();
    // transformed3.GetMinimas();
    // cout << 5 << endl;



    mynoise.GenFields(func,str,1.0/2.0);

    transformed3.Add_Noise(mynoise);


    // cout << transformed1.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[1][0] << endl;
    // cout << transformed1.calculated_reactions[2][0] << endl;

    // cout << endl;
    rules.Calculate_Results(transformed1.calculated_reactions, transformed2.calculated_reactions, transformed3.calculated_reactions);

    
    // cout << "Rules done" << endl;
    // rules.GetMaximas();
    // cout << transformed1.calculated_reactions[0][0] << endl;
    // cout << rules.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[0][1] << endl;
    // cout << rules.calculated_reactions[0][1] << endl;

    // cout << transformed1.calculated_reactions[0][2] << endl;
    // cout << rules.calculated_reactions[0][2] << endl;


    // cout << endl;
    //  cout << 6 << endl;
    reverse_transform.Calculate_Results(rules.calculated_reactions);



    //cout << "beginning renorm" << endl;
    for(int i = 0  ; i < myp.number_of_fields ; i++)
    for(int j = 0 ; j < myp.get_total() ; j++ )
        reverse_transform.calculated_reactions[i][j] = {reverse_transform.calculated_reactions[i][j].real(),0.};
    //cout << "renorm" << endl;
    // outfunc(fields[5], filename1, myp);
    // outfunc(weigs.calculated_reactions[5], filename2, myp);
    // outfunc(transformed1.calculated_reactions[5], filename3, myp);
    // outfunc(transformed2.calculated_reactions[5], filename4, myp);
    // outfunc(rules.calculated_reactions[5], filename5, myp);
    // outfunc(reverse_transform.calculated_reactions[5], filename6, myp);
    // cout << "fields: ";
    // for (int i = 0; i < myp.number_of_fields; i++)
    //     cout << fields[i][0] << ",";
    // cout << endl;
    // cout << "updated: ";
    // for(int i = 0  ; i < myp.number_of_fields ; i++)
    // cout << reverse_transform.calculated_reactions[i][0] <<",";
    // cout << endl;
    // pausel();

    cout << "Reverse Transform done" << endl;
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;
    cout << 7 << endl;
    // pausel();
    // cout << reverse_transform.calculated_reactions[0][0] << endl;

    set_field(reverse_transform.calculated_reactions);
}

void CHSubFrac::Update()
{


    fractional_weight.Calculate_Results(fields, initial_cond.calculated_reactions, initial_cond.calculated_reactions);


    weigs.Calculate_Results(fields);


    chems.Calculate_Results(fields);


    transformed1.Calculate_Results(fractional_weight.calculated_reactions);


    transformed2.Calculate_Results(weigs.calculated_reactions);


    transformed3.Calculate_Results(chems.calculated_reactions);

    rules.Calculate_Results(transformed1.calculated_reactions, transformed2.calculated_reactions, transformed3.calculated_reactions);

    reverse_transform.Calculate_Results(rules.calculated_reactions);

    cout << "Reverse Transform done" << endl;
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;
    cout << 7 << endl;


    set_field(reverse_transform.calculated_reactions);
}

void CHFrac::Update() {
    // cout << "calculating weight" << endl;


    fractional_weight.Calculate_Results(fields,old_fields.calculated_reactions,initial_cond.calculated_reactions);

    // cout << "weight? calculated" << endl;

    weigs.Calculate_Results(fields);

    // cout << "fractionalweight: " << endl;
    // cout << fractional_weight.calculated_reactions[0][0] << endl;
    // cout << fractional_weight.calculated_reactions[1][0] << endl;
    // cout << fractional_weight.calculated_reactions[2][0] << endl;
    // cout << fractional_weight.calculated_reactions[3][0] << endl;
    // pausel();
    // cout << "weigs done" << endl;
    // weigs.GetMaximas();
    // weigs.GetMinimas();

    chems.Calculate_Results(fields);

    // cout << "chems: " << endl;
    // cout << chems.calculated_reactions[0][0] << endl;
    // cout << chems.calculated_reactions[1][0] << endl;
    // cout << chems.calculated_reactions[2][0] << endl;
    // cout << chems.calculated_reactions[3][0] << endl;

    // pausel();
    // cout << "chemical dynamics: ";
    // for (int i = 0; i < myp.number_of_fields; i++)
    //     cout << chems.calculated_reactions[i][0] << ",";
    // cout << endl;
    // cout << "chems done" << endl;
    // chems.GetMaximas();
    // chems.GetMinimas();

    transformed1.Calculate_Results(fractional_weight.calculated_reactions);

    // cout << "FFT1 done" << endl;
    // transformed1.GetMaximas();
    // transformed1.GetMinimas();
    // cout << 3 << endl;

    transformed2.Calculate_Results(weigs.calculated_reactions);

    // cout << "FFT2 done" << endl;
    // transformed2.GetMaximas();
    // transformed2.GetMinimas();
    // cout << 4 << endl;

    transformed3.Calculate_Results(chems.calculated_reactions);

    rules.Calculate_Results(transformed1.calculated_reactions, transformed2.calculated_reactions, transformed3.calculated_reactions);

    // cout << transformed1.calculated_reactions[0][0] << endl;
    // cout << transformed2.calculated_reactions[0][0] << endl;
    // cout << transformed3.calculated_reactions[0][0] << endl;
    cout << rules.calculated_reactions[0][0] << " " << rules.calculated_reactions[1][0] << " " << rules.calculated_reactions[2][0] << " " << rules.calculated_reactions[3][0] << endl;
    
    reverse_transform.Calculate_Results(rules.calculated_reactions);

    cout << "Reverse Transform done" << endl;
    reverse_transform.GetMaximas();
    reverse_transform.GetMaximasIndex();
    reverse_transform.GetMinimas();
    reverse_transform.GetMinimasIndex();
    cout << endl;
    cout << 7 << endl;

    old_fields.set_field(fields);

    set_field(reverse_transform.calculated_reactions);
}

#endif /* CAHNHILLIARD_CPP */
