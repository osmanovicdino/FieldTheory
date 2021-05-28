#ifndef CAHNHILLIARD_CPP
#define CAHNHILLIARD_CPP

void outfunc(double *a, string s, const CH_builder &p)
{
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if (!myfileg.is_open())
        error("failed to open file");
    //myfileg <<= a;
    int nr = p.N1;
    int nc = p.N2;

    for (int i = 0; i < nr; ++i)
    {

        for (int j = 0; j < nc; ++j)
        {
            j == nc - 1 ? myfileg << a[i * nc + j] << endl : myfileg << a[i * nc + j] << ",";
        }
    }

    myfileg.close();
}


CH::CH(const CH_builder &p) : myp(p), chems(Field_Wrapper(p)) ,  
weigs(Field_Wrapper(p)), transformed1(Field_Wrapper(p)) ,
transformed2(Field_Wrapper(p)), transformed3(Field_Wrapper(p)) ,
 rules(Rule_Wrapper(p)), reverse_transform(Field_Wrapper(p))
{

    fields = new double *[myp.number_of_fields];

    for (int i = 0; i < myp.number_of_fields; i++)
    {
        fields[i] = (double *)fftw_malloc(myp.N1 * myp.N2 * sizeof(double));
        for (int j = 0; j < myp.N1 * myp.N2; j++)
        {
            fields[i][j] = 0.0;
        }
    }

    for(int i = 0  ; i < myp.number_of_fields ; i++) {
        FourierWeightForward fw;
        transformed1.add_method(fw, i);
        transformed2.add_method(fw, i);
        transformed3.add_method(fw, i);
        FourierWeightBackward fw2;
        reverse_transform.add_method(fw2,i);
    }

}

CH::~CH() {
    for (int i = 0; i < myp.number_of_fields; i++)
    {
        fftw_free(fields[i]);
    }
    delete fields;

}

void CH::set_field(const matrix<double> &mymat, int k) {
    for(int i =  0 ;  i < myp.N1 ; i++) {
        for(int j = 0  ; j < myp.N2 ; j++) {
            fields[k][i*myp.N2+j] =  mymat.gpcons(i,j);
        }
    }

}

void CH::set_field(double **orig) {
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





void CH::Update() {
    //cout << fields[0][0] << endl;
    //cout << "trying" << endl;

    weigs.Calculate_Results(fields);

    // cout << 1 << endl;

    chems.Calculate_Results(fields);

    // cout << 2 << endl;

    transformed1.Calculate_Results(fields);


    // cout << 3 << endl;

    transformed2.Calculate_Results(weigs.calculated_reactions);

    // cout << 4 << endl;

    transformed3.Calculate_Results(chems.calculated_reactions);

    // cout << 5 << endl;

    // cout << transformed1.calculated_reactions[0][0] << endl;
    // string filename2 = "weig_transform";
    // string fn2 = "chem_transform";
    // outfunc(transformed2.calculated_reactions[0], filename2, myp);
    // outfunc(transformed3.calculated_reactions[0], fn2, myp);
    // cout << transformed1.calculated_reactions[1][0] << endl;
    // cout << transformed1.calculated_reactions[2][0] << endl;

    // cout << endl;
    rules.Calculate_Results(transformed1.calculated_reactions,transformed2.calculated_reactions,transformed3.calculated_reactions);

    // cout << transformed1.calculated_reactions[0][0] << endl;
    // cout << rules.calculated_reactions[0][0] << endl;

    // cout << transformed1.calculated_reactions[0][1] << endl;
    // cout << rules.calculated_reactions[0][1] << endl;

    // cout << transformed1.calculated_reactions[0][2] << endl;
    // cout << rules.calculated_reactions[0][2] << endl;
    // string filename = "rule_transform";
    // string fn = "field_transform";
    // outfunc(rules.calculated_reactions[0], filename, myp);
    // outfunc(transformed1.calculated_reactions[0], fn, myp);

    // pausel();

    // cout << endl;
    // cout << 6 << endl;
    reverse_transform.Calculate_Results(rules.calculated_reactions);

    // cout << 7 << endl;
    //cout << reverse_transform.calculated_reactions[0][0] << endl;

    set_field(reverse_transform.calculated_reactions);

    // cout << 8 << endl;
}


void CH::print_all_results(string s1) {
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

void Setup_Invasion(const CH &a) {

}


#endif /* CAHNHILLIARD_CPP */
