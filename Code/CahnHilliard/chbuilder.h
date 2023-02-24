#ifndef CHBUILDER_H
#define CHBUILDER_H

struct CH_builder
{

    int number_of_fields;
    int dimension = 2; //baseline
    int N1;
    int N2;
    int N3;

    int get_total() const { if(dimension==2) return N1*N2; else return N1*N2*N3; }
};

void outfunc(double *a, string s, const CH_builder &p)
{
    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if (!myfileg.is_open())
        error("failed to open file");
    //myfileg <<= a;
    if( p.dimension == 2) {
    int nr = p.N1;
    int nc = p.N2;

    for (int i = 0; i < nr; ++i)
    {

        for (int j = 0; j < nc; ++j)
        {
            j == nc - 1 ? myfileg << a[i * nc + j] << endl : myfileg << a[i * nc + j] << ",";
        }
    }
    }
    else{
    int nr = p.N1;
    int nc = p.N2;
    int nd = p.N3;

    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            for (int k = 0; k < nd; ++j)
            {
                k == nd - 1 ? myfileg << a[i * nc * nd + j * nd + k] << endl : myfileg << a[i * nc * nd + j * nd + k] << ",";
            }
        }
    }
    }

    myfileg.close();
}

void outfunc(complex<double> *a, string s, const CH_builder &p)
{
    string s1 = s + "_real.csv";
    ofstream myfileg;
    myfileg.open(s1.c_str());
    if (!myfileg.is_open())
        error("failed to open file");
    //myfileg <<= a;
    int nr = p.N1;
    int nc = p.N2;
    if (p.dimension == 2)
    {
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                j == nc - 1 ? myfileg << std::setprecision(3) << a[i * nc + j].real() << endl : myfileg << std::setprecision(3) << a[i * nc + j].real() << ",";
            }
        }
    }
    else{
        int nr = p.N1;
        int nc = p.N2;
        int nd = p.N3;

        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                for (int k = 0; k < nd; ++k)
                {
                k == nd - 1 ? myfileg << std::setprecision(3) << a[i * nc * nd + j * nd + k].real() << endl : myfileg << std::setprecision(3) << a[i * nc * nd + j * nd + k].real() << ",";
                }
            }
        }
    }

    // string s2 = s + "_imag.csv";
    // ofstream myfileg2;
    // myfileg2.open(s2.c_str());
    // for (int i = 0; i < nr; ++i)
    // {

    //     for (int j = 0; j < nc; ++j)
    //     {
    //         j == nc - 1 ? myfileg2 << a[i * nc + j].imag() << endl : myfileg2 << a[i * nc + j].imag() << ",";
    //     }
    // }

    myfileg.close();
    // myfileg2.close();
    }

#endif /* CHBUILDER_H */
