#ifndef CHBUILDER_H
#define CHBUILDER_H

struct CH_builder
{

    int number_of_fields;
    int N1;
    int N2;

    int get_total() const { return N1*N2; }
};

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

    for (int i = 0; i < nr; ++i)
    {

        for (int j = 0; j < nc; ++j)
        {
            j == nc - 1 ? myfileg << std::setprecision(3) << a[i * nc + j].real() << endl : myfileg << std::setprecision(3) << a[i * nc + j].real() << ",";
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
