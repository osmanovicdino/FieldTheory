#ifndef CAHNHILLIARD_CPP
#define CAHNHILLIARD_CPP

template<class T>
CH<T>::CH(const CH_builder &p) : myp(p), chems(Field_Wrapper<T, T>(p)),
                              weigs(Field_Wrapper<T,T>(p)), transformed1(Field_Wrapper<T, T>(p)),
                              transformed2(Field_Wrapper<T,T>(p)), transformed3(Field_Wrapper<T, T>(p)),
                              rules(Rule_Wrapper<T, T,T, T>(p)), reverse_transform(Field_Wrapper<T, T>(p))
{

    fields = new T *[myp.number_of_fields];

    for (int i = 0; i < myp.number_of_fields; i++)
    {
        fields[i] = (T *)fftw_malloc(myp.N1 * myp.N2 * sizeof(T));
        for (int j = 0; j < myp.N1 * myp.N2; j++)
        {
            fields[i][j] = 0.0;
        }
    }

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
void CH<T>::set_field(T **orig) {
    for(int k = 0 ; k < myp.number_of_fields ; k++) {
        for (int i = 0; i < myp.N1; i++)
        {
            for (int j = 0; j < myp.N2; j++)
            {
                fields[k][i * myp.N2 + j] = orig[k][i * myp.N2 + j];
            }
        }
    }
}




template<class T>
void CH<T>::Update() {
    //cout << fields[0][0] << endl;
    //cout << "trying" << endl;

    weigs.Calculate_Results(fields);

    // cout << "weigs done" << endl;
    // weigs.GetMaximas();
    // weigs.GetMinimas();

    chems.Calculate_Results(fields);

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

    // cout << transformed1.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[1][0] << endl;
    // cout << transformed1.calculated_reactions[2][0] << endl;

    // cout << endl;
    rules.Calculate_Results(transformed1.calculated_reactions,transformed2.calculated_reactions,transformed3.calculated_reactions);


    rules.Check_fields();
    // cout << "Rules done" << endl;
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
    cout << "beginning renorm" << endl;
    for (int i = 0; i < myp.number_of_fields; i++)
        for (int j = 0; j < myp.get_total(); j++)
            reverse_transform.calculated_reactions[i][j] = {reverse_transform.calculated_reactions[i][j].real(), 0.};
    cout << "renorm" << endl;
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
     cout << 7 << endl;
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

    chems.Calculate_Results(fields);

    
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

#endif /* CAHNHILLIARD_CPP */
