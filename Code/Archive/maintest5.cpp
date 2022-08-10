#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/ioctl.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <random>
//#include <mutex>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <chrono>
#include <type_traits>
//#include <stdafx>
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_threads() { return 1; }
#endif

#include "DataStructures/basic.h"
#include "DataStructures/vector1.h"
#include "DataStructures/matrix2.h"
#include "DataStructures/matrix2.cpp"
#include "CahnHilliard/cahnhilliard.h"

#include "fftw3.h"

using namespace std;

int main(int argc, char **argv)
{
    srand(time(NULL));
    /*
    CH_builder p;
    p.number_of_fields = 4;
    p.N1 = 1024;
    p.N2 = 1024;
    Field_Wrapper a(p);
    Field_Wrapper a2(p);

    double T;
    bool err1, err2, err3, err4;
    matrix<double> field1 = importcsv("ci.csv", T, err1);
    matrix<double> field2 = importcsv("Ii.csv", T, err2);
    matrix<double> field3 = importcsv("ρi.csv", T, err3);
    matrix<double> field4 = importcsv("Ai.csv", T, err4);

    CH b(p);

    b.set_field(field1, 0);
    b.set_field(field2, 1);
    b.set_field(field3, 2);
    b.set_field(field4, 3);

    for (int i = 0; i < p.number_of_fields; i++)
    {
        FourierWeightForward fw;
        a.add_method(fw, i);
        FourierWeightBackward fw2;
        a2.add_method(fw2, i);
    }


    cout << b.fields[0][0] << endl;
    a.Calculate_Results(b.fields);

    cout << a.calculated_reactions[0][0] << endl;
    a2.Calculate_Results(a.calculated_reactions);

    cout << a2.calculated_reactions[0][0] << endl;
    */

    // double p0 = atof(argv[1]);
    // double rate1 = atof(argv[2]);

    // CH_builder p;
    // int nof = 1;
    // p.number_of_fields = nof;
    // p.N1 = 1024;
    // p.N2 = 1024;

    // Field_Wrapper<complex<double>, complex<double>> revk1(p);
    // Field_Wrapper<complex<double>, complex<double>> revk2(p);

    // FourierWeightBackwardGradient2D a1;
    // a1.set_direction(0);
    // FourierWeightBackwardGradient2D a2;
    // a2.set_direction(1);

    // revk1.add_method(a1, 0);
    // revk2.add_method(a2, 0);

    // complex<double> **field = new complex<double> *[p.number_of_fields];

    // field[0] = (complex<double> *)fftw_malloc(p.N1 * p.N2 * sizeof(complex<double>));

    // Field_Wrapper<complex<double>, complex<double>> fourier1(p);
    // Field_Wrapper<complex<double>, complex<double>> fourier2(p);
    // FourierWeightForward2D a3;
    // fourier1.add_method(a3, 0);
    // fourier2.add_method(a3, 0);

    // for (int i = 0; i < p.N1; i++)
    // {
    //     for (int j = 0; j < p.N2; j++)
    //     {
    //         field[0][i * p.N2 + j] = sin(2. * pii * i / 1024.) * sin(2 * pii * (j / 1024.)); // realspace
    //     }
    // }

    // cout << field[0][300000] << endl;
    // Field_Wrapper<complex<double>, complex<double>> invf(p);
    // FourierWeightBackward2D a4;

    // invf.add_method(a4, 0);


    // fourier1.Calculate_Results(field);

    // invf.Calculate_Results(fourier1.calculated_reactions);

    // cout << invf.calculated_reactions[0][300000] << endl;

/* 

    
    CH_builder p;
    int nof = 1;
    p.number_of_fields = nof;
    p.N1 = 1024;
    p.N2 = 1024;

    Field_Wrapper<complex<double>, complex<double> > revk1(p);
    Field_Wrapper<complex<double>, complex<double> > revk2(p);

    FourierWeightBackwardGradient2D a1;
    a1.set_direction(0);
    FourierWeightBackwardGradient2D a2;
    a2.set_direction(1);


    revk1.add_method(a1, 0);
    revk2.add_method(a2, 0);

    complex<double> **field = new complex<double> *[p.number_of_fields];

    field[0] = (complex<double> *)fftw_malloc(p.N1 * p.N2 * sizeof(complex<double>));

    Field_Wrapper<complex<double>, complex<double>> fourier1(p);
    Field_Wrapper<complex<double>, complex<double>> fourier2(p);
    Field_Wrapper<complex<double>, complex<double>> fourier3(p);
    FourierWeightForward2D a3;
    fourier1.add_method(a3, 0);
    fourier2.add_method(a3, 0);
    fourier3.add_method(a3, 0);

    for(int i = 0  ; i < p.N1 ; i++) {
        for(int j = 0  ; j < p.N2 ; j++) {
            field[0][i * p.N2 + j] = 1.+log((double)rand() / (double)RAND_MAX); // realspace
        }
    }

   
    
    fourier1.Calculate_Results(field); //k space

    string file_fourier = "file_fourier";
    outfunc(fourier1.calculated_reactions[0], file_fourier, p);


    string file1 = "file1";
    string file2 = "file2";
    string file3 = "file3";

 

    revk1.Calculate_Results(fourier1.calculated_reactions); // real space

    // outfunc(revk1.calculated_reactions[0], file2,p);



    revk2.Calculate_Results(fourier1.calculated_reactions); // real space

    // outfunc(revk2.calculated_reactions[0], file3, p);

    fourier2.Calculate_Results(revk1.calculated_reactions); //k space
    fourier3.Calculate_Results(revk2.calculated_reactions);

    outfunc(fourier2.calculated_reactions[0], file2, p);
    outfunc(fourier3.calculated_reactions[0], file3, p);

    complex<double> **res = new complex<double> *[p.number_of_fields];
    res[0] = (complex<double> *)fftw_malloc(p.N1 * p.N2 * sizeof(complex<double>));


        for (int i1 = 0; i1 < p.N1; i1++)
        {
            for (int j1 = 0; j1 < p.N2; j1++)
            {

                double k1, k2;
                if (i1 <= p.N1/2)
                {
                    k1 = (2 * pii) * i1;
                }
                else
                {
                    k1 = (2 * pii) * (i1 - 1.0*p.N1);
                }
                if (j1 <= p.N2/2)
                {
                    k2 = (2 * pii) * j1;
                }
                else
                {
                    k2 = (2 * pii) * (j1 - 1.0*p.N2);
                }
                complex<double> myfac1 = {0.0, k1};
                complex<double> myfac2 = {0.0, k2};

                //we need a negative sign here for the two differentials we are doing


                res[0][i1 * p.N2 + j1] = myfac1 * fourier2.calculated_reactions[0][i1 * p.N2 + j1] + myfac2 * fourier3.calculated_reactions[0][i1 * p.N2 + j1];
            }
        }

        int ss = 0;
        int ds = 0;
        for (int i1 = 0; i1 < p.N1; i1++)
        {
            for (int j1 = 0; j1 < p.N2; j1++)
            {

                double k1, k2;
                if (i1 <= p.N1 / 2)
                {
                    k1 = (2 * pii) * i1;
                }
                else
                {
                    k1 = (2 * pii) * (i1 - 1.0 * p.N1);
                }
                if (j1 <= p.N2 / 2)
                {
                    k2 = (2 * pii) * j1;
                }
                else
                {
                    k2 = (2 * pii) * (j1 - 1.0 * p.N2);
                }
                complex<double> myfac1 = {0.0, k1};
                complex<double> myfac2 = {0.0, k2};


                if (sign(fourier1.calculated_reactions[0][i1 * p.N2 + j1].real()) == sign(res[0][i1 * p.N2 + j1].real()))
                {
                    // cout << myfac1 << endl;
                    // cout << myfac2 << endl;
                    // cout << fourier1.calculated_reactions[0][i1 * p.N2 + j1] << endl;
                    // cout << res[0][i1 * p.N2 + j1] << endl;

                    ss += 1;
                }
                else
                    ds += 1;
            }
        }
        cout << ss << endl;
        cout << ds << endl;
        pausel();
        string fileres = "fileres";
        outfunc(res[0], fileres, p);



        Field_Wrapper<complex<double>, complex<double>> invf(p);
        FourierWeightBackward2D a4;

        invf.add_method(a4, 0);

        invf.Calculate_Results(res);

        outfunc(invf.calculated_reactions[0], file1, p);

        pausel(); */




    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc,myc> FWCC; 
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    CH_builder p;
    int nof = 7;
    p.number_of_fields = nof;
    p.N1 = 128;
    p.N2 = 128;

    CHN a(p);

    double ratef1 = 0.0;
    double rateb1 = 0.0;

    ILB c0(ratef1, 0.0 * ratef1 /*no growth from below as this is the lowest field*/, 0.0 * rateb1 /*no decay down*/, rateb1, 0, 0, 1, 5);
    ILB c1(ratef1, ratef1, rateb1, rateb1, 1, 0, 2, 5);
    ILB c2(ratef1, ratef1, rateb1, rateb1, 2, 1, 3, 5);
    ILB c3(ratef1, ratef1, rateb1, rateb1, 3, 2, 4, 5);
    ILB c4(0.0 * ratef1 /*no loss to above*/, ratef1, rateb1, 0.0 * rateb1 /*no decay above*/, 4, 3, 4, 5);

    MultipleReactions<myc> c5(4);
    ILA d0(ratef1, rateb1, 0, 5, 1);
    ILA d1(ratef1, rateb1, 1, 5, 2);
    ILA d2(ratef1, rateb1, 2, 5, 3);
    ILA d3(ratef1, rateb1, 3, 5, 4);

    c5.add_chemical_reaction(d0, 0);
    c5.add_chemical_reaction(d1, 1);
    c5.add_chemical_reaction(d2, 2);
    c5.add_chemical_reaction(d3, 3);

    FWCC my_chemsitry(p);

    NoWeight<myc,myc> c6; // the solvent;

    my_chemsitry.add_method(c0, 0);
    my_chemsitry.add_method(c1, 1);
    my_chemsitry.add_method(c2, 2);
    my_chemsitry.add_method(c3, 3);
    my_chemsitry.add_method(c4, 4);
    my_chemsitry.add_method(c5, 5);
    my_chemsitry.add_method(c6, 6);

    cout << "added chemistries" << endl;

    // myc ** calculated_reactions = new myc *[p.number_of_fields];

    // for(int k =0 ; k < p.number_of_fields ; k++) {
    // calculated_reactions[k] = (myc *)fftw_malloc(p.N1 * p.N2 * sizeof(myc));
    // for (int j = 0; j < p.N1 * p.N2; j++)
    // {
    //     calculated_reactions[k][j] = 1.0;
    // }
    // }
    // my_chemsitry.Calculate_Results(calculated_reactions);

    // cout << "up to here" << endl;

    FWCC my_weights(p);

    matrix<double> myint(nof,nof);

    double energy_scale = -0.25;
    for(int i = 0  ; i < nof ; i++) {
        for(int j = 0  ; j < nof ; j++) {
            if(i > 4 || j > 4) myint(i,j) = 0.0;
            else myint(i,j) = energy_scale*(4-i)*(4-j);
        }
    }

    cout << "int done" << endl;

    vector1<double> specvolumes(nof,1.);

    specvolumes[nof-1] = 1.;


    

    CPSS w0(nof, myint, specvolumes);
    CPSS w1(nof, myint, specvolumes);
    CPSS w2(nof, myint, specvolumes);
    CPSS w3(nof, myint, specvolumes);
    CPSS w4(nof, myint, specvolumes);
    CPSS w5(nof, myint, specvolumes);
    CPSS w6(nof, myint, specvolumes);

    my_weights.add_method(w0, 0);
    my_weights.add_method(w1, 1);
    my_weights.add_method(w2, 2);
    my_weights.add_method(w3, 3);
    my_weights.add_method(w4, 4);
    my_weights.add_method(w5, 5);
    my_weights.add_method(w6, 6);

    cout << "weights added" << endl;

    a.set_chems(my_chemsitry);
    a.set_weights(my_weights);

    cout << "everything set" << endl;
    matrix<double> ChiMatrix = w0.getchi();


    double dt = 0.5;
    //double dx = 0.05;

    //double temp1 = (1. / (dx * p.N1)) ;

    double surf = 100.0;
    double temp1 = 0.01*(1./(2*pii*pow(surf,0.25)*pow(dt,0.25)*p.N1));

    //cout << (1. / (dx * p.N1)) << endl;
    cout << temp1 << endl;
    pausel();

    vector1<double> Di(7,1.);

    double diffusion_constant_invader =  1.;
    Di[5] = diffusion_constant_invader;

    //double surface_width = 2.0;

    RWC my_rules(p);
    matrix<double> surface_width(nof,nof,1.);
    // for(int k = 0 ; k < nof ; k++) {
    //     surface_width(k,k) = 0.6;
    // }

    MultiPhaseBundleComplex Bund(w0, p, dt, Di, temp1, surface_width, surf);

    // MultiPhase r1(p, 6, dt, Di, temp1, surface_width, ChiMatrix);
    // MultiPhase r2(p);
    // MultiPhase r3(p);
    // MultiPhase r4(p);
    // MultiPhase r5(p);
    // MultiPhase r6(p); //all empty

    // cout << "created diffusion" << endl;

    //for(int k = 0  ; k < nof ; k++)
    my_rules.add_method(Bund.C, 0); //as this is a collective one, only need to define the first rule;
    

    NoRule<myc,myc,myc,myc> noru(p);
    for(int k = 1; k < nof ; k++)
        my_rules.add_method(noru,k);

    cout << "all rules added" << endl;
    // cout << "all methods added" << endl;
    cout << ChiMatrix << endl;
    pausel();

    a.set_rules(my_rules);

    // cout << "Setting rules" << endl;

    cout << "old rules set" << endl;

    RWC new_rules(p);

    new_rules.add_method(Bund.B, 0); //collective rule

    for(int k = 1  ; k < nof ; k++)
    new_rules.add_method(noru, k);

    cout << "all new methods added" << endl;

    a.set_new_rules(new_rules);

    cout << "all rules set" << endl;

    // for(int k = 0  ; k < nof ; k++)
    //     a.rules.chems[k]->print();

    // for (int k = 0; k < nof; k++)
    //     a.newrules.chems[k]->print();

    // cout << "added methods" << endl;

    matrix<myc> field1(p.N1, p.N2);
    matrix<myc> field2(p.N1, p.N2);
    matrix<myc> field3(p.N1, p.N2);

    double dens = 0.4;
    double densi = 0.01;
    for (int i = 0; i < p.N1; i++)
    {
        for (int j = 0; j < p.N2; j++)
        {
            double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            double r2 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
            field1(i, j) = dens + 0.1 * r1;
            field2(i, j) = densi;// + 0.001 * r2;
            field3(i, j) = 1. - field1(i,j) - 5.*field2(i,j);
        }
    }

    // double val;
    // (field1+(5.*field2)+field3).maxima(val);
    // if(val>1.) error("initial packing fracs not correct");

    a.set_field(field1, 0);
    a.set_field(field2, 1);
    a.set_field(field2, 2);
    a.set_field(field2, 3);
    a.set_field(field2, 4);
    a.set_field(field2, 5);
    a.set_field(field3, 6);

    cout << "all fields set" << endl;

    // auto start = std::chrono::high_resolution_clock::now();
    int runtime = 100000;
    int every = 1000;

    int tf = ceil((double)runtime / (double)every);
    int number_of_digits = 0;
    do
    {
        ++number_of_digits;
        tf /= 10;
    } while (tf);

    // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runtime; i++)
    {

        if (i % every == 0)
        {
            stringstream strep1;
            stringstream strep2;
            stringstream strep3;
            stringstream strep4;

            strep1 << dens;
            strep4 << densi;
            strep2 << ratef1;
            strep3 << rateb1;

            string s1 = "denp=" + strep1.str() + "denpi=" + strep4.str() + "_ratef1=" + strep2.str() + "_rateb1=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();

            a.print_all_results(s1 + s2);
        }
        cout << i << endl;
        a.Update();
    }

    cout << "did an update" << endl;
    pausel();

    // auto stop = std::chrono::high_resolution_clock::now();

    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // cout << duration.count() << endl;

    // cout << "he's done ya again" << endl;

    /*
    int N1 = 512;
    int N2 = 512;
    double **an_in_array = new double *[4];
    double **an_out_array = new double *[4];
    an_in_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));
    an_out_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));

    // fftw_plan p;

    auto start = std::chrono::high_resolution_clock::now();

    // p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

    for(int i = 0 ; i < 1000 ; i++) {
        fftw_plan p;

        p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

        fftw_execute(p);

        fftw_destroy_plan(p);
    }

    auto stop =  std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

    cout << duration.count() << endl;
    */

    // fftw_destroy_plan(p);

    /*
CH_builder p;
p.number_of_fields = 4;
p.N1=1024;
p.N2=1024;



Field_Wrapper a(p);


InvasionAntiInvasionLinear b(1.0,0.5);
NoWeight c;

a.add_method(b,0);
a.add_method(c, 1);
a.add_method(c, 2);
a.add_method(c, 3);

double **fields = new double *[4];
fields[0] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[1] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[2] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));
fields[3] = (double *)fftw_malloc(1024 * 1024 * sizeof(double));

for(int i = 0 ; i < 4 ; i++) {
    for(int j = 0 ; j < SQR(1024) ; j++) {
        fields[i][j] = 0.2;
    }
}

cout << "fine to here" << endl;

a.Calculate_Results(fields);

cout << a.calculated_reactions[0][0] << endl;
cout << "done" << endl;
pausel();
*/
    //Reaction_Wrapper()

    // CH a(p);

    // a.chems[0]->what_am_I();
    // a.update_rules[0]->what_am_I();

    // InvasionAntiInvasion *b = new InvasionAntiInvasion;

    // cout << a.chems.size() << endl;

    // a.chems[0] =  b;

    // a.chems[0]->what_am_I();

    /*
    int N1 = 512;
    int N2 = 512;
    double **an_in_array = new double*[4];
    double **an_out_array = new double*[4];
    an_in_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));
    an_out_array[0] = (double *)fftw_malloc(N1 * N2 * sizeof(double));

    fftw_plan p;

    p = fftw_plan_r2r_2d(N1, N2, an_in_array[0], an_out_array[0], FFTW_REDFT10, FFTW_REDFT10, 1);

    fftw_execute(p);
    */

    /*
int N1 = 512;
int N2 = 512;



double *an_array;
an_array = (double *)fftw_malloc(N1 * N2 * sizeof(double));

double *out;
out = (double *)fftw_malloc(N1 * N2 * sizeof(double));
// in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
// out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

for (int i = 0; i < N1; i++)
{
    for (int j = 0; j < N2; j++)
    {
        an_array[i * N1 + j] = cos(pi*i*j/10.);
    }
}

fftw_plan p;

p = fftw_plan_r2r_2d(N1, N2, an_array, out, FFTW_REDFT10, FFTW_REDFT10, 1);

fftw_execute(p);
cout << an_array[0] << endl;
cout << out[0] / (4 * 512) << endl;




outfunc2D(out,N1,N2,"cosine");

fftw_plan p2;
p2 = fftw_plan_r2r_2d(N1, N2, out, an_array, FFTW_REDFT01, FFTW_REDFT01, 1);

fftw_execute(p2);

cout << an_array[0]/SQR(2*512) << endl;
cout << an_array[1] << endl;
cout << out[0] / (4*512) << endl;
*/

    // cout << a <<endl;
}