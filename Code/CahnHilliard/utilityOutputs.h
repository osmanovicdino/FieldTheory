
#ifndef UTILITYOUTPUTS_H
#define UTILITYOUTPUTS_H

#include "../DataStructures/basic.h"
#include "../DataStructures/vector1.h"
#include "../DataStructures/vector1.cpp"
#include "../DataStructures/matrix2.h"

void outfunc2D(double *in, int N0, int N1, string s) {
    //for a filename s and dimensions N0 and N1

    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if (!myfileg.is_open())
        error("failed to open file");

    for (int i = 0; i < N0; ++i)
    {

        for (int j = 0; j < N1; ++j)
        {
            j == N1 - 1 ? myfileg << in[i * N0 + j] << endl : myfileg << in[i * N0 + j] << ",";
        }
    }

    //myfileg <<= a;
    myfileg.close();
}

void outfunc2D(complex<double>  *in, int N0, int N1, string s)
{
    // for a filename s and dimensions N0 and N1

    s = s + ".csv";
    ofstream myfileg;
    myfileg.open(s.c_str());
    if (!myfileg.is_open())
        error("failed to open file");

    for (int i = 0; i < N0; ++i)
    {

        for (int j = 0; j < N1; ++j)
        {
            j == N1 - 1 ? myfileg << in[i * N0 + j].real() << endl : myfileg << in[i * N0 + j].real() << ",";
        }
    }

    // myfileg <<= a;
    myfileg.close();
}

#endif