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
    // double A0 = atof(argv[2]);

    // double eps1 = -0.8;
    // double eps2 = 0.8;

    CH_builder p;
    p.number_of_fields = 3;
    p.N1 = 1024;
    p.N2 = 1024;

    CH a(p);

    double rate1 = 0.15;
    // InvasionSelfRelease c1(rate1, rate2, 0, 1, 2);
    // InvasionSelfRelease c2(-rate1, -rate2, 0, 1, 2);
    // InvasionSelfRelease c3(rate1, rate2, 0, 1, 2);
    //NoWeight c4;

    InvasionLinear c1(rate1, 0, 2);
    InvasionLinear c2(-rate1, 0, 2);
    InvasionLinear c3(rate1, 0, 2);

    Field_Wrapper my_chemsitry(p);
    my_chemsitry.add_method(c1, 0);
    my_chemsitry.add_method(c2, 1);
    my_chemsitry.add_method(c3, 2);

    Field_Wrapper my_weights(p);

    double c00 = 0.2;
    double c11 = 0.8;
    double nu = 0.871687587;

    double eps1 = 2.00;
    double eps2 = 0.0;

    CahnHilliardWithCouplingWeightSQR d1(c00, c11, nu, eps1, eps2, 1, 0);
    DiffDiffusiveWeightSQR d2(eps1, 0);
    DiffusiveWeight d3(1.0);

        // DiffusiveWeight d2(1.0);
        // DiffDiffusiveWeightSQR d3(eps1, 0);

    my_weights.add_method(d1, 0);
    my_weights.add_method(d2, 1);
    my_weights.add_method(d3, 2);

    a.set_chems(my_chemsitry);
    a.set_weights(my_weights);

    cout << "weights added" << endl;

    Rule_Wrapper my_rules(p);
    cout << "rule wrapper created" << endl;

    double dt = 0.05;
    double dx = 0.05;
    double D = 0.1;
    double temp1 = SQR(pi / (dx * p.N1));
    double eps = 0.3;

    DiffusionWithSurfaceTension e1(p, dt, D, temp1, eps);
    DiffusionWithInteraction e2(p, dt, D, temp1);
    NormalDiffusion e3(p, dt, D, temp1);
    

    cout << "created diffusion" << endl;

    my_rules.add_method(e1, 0);
    cout << "method1" << endl;
    my_rules.add_method(e2, 1);
    cout << "method2" << endl;
    my_rules.add_method(e3, 2);
    cout << "method3" << endl;

    cout << "added methods" << endl;

    a.set_rules(my_rules);

    cout << "done" << endl;

    // matrix<double> field1 = importcsv("ci.csv", T, err1);
    // matrix<double> field2 = importcsv("Ii.csv", T, err2);
    // matrix<double> field3 = importcsv("ρi.csv", T, err3);
    // matrix<double> field4 = importcsv("Ai.csv", T, err4);

    matrix<double> field1(p.N1, p.N2);
    matrix<double> field2(p.N1, p.N2);
    matrix<double> field3(p.N1, p.N2);

    double bc = 0.166;
    double p0 = 0.596;

    for (int i = 0; i < p.N1; i++)
    {
        for (int j = 0; j < p.N2; j++)
        {
            field1(i, j) = bc + (0.1 * ((double)rand() / (double)RAND_MAX) - 0.1);
            field2(i, j) = 0.0;
            field3(i, j) = p0 + (0.2 * ((double)rand() / (double)RAND_MAX) - 0.1);
            //field4(i, j) = A0 + (0.2 * ((double)rand() / (double)RAND_MAX) - 0.1);
        }
    }

    a.set_field(field1, 0);
    a.set_field(field2, 1);
    a.set_field(field3, 2);

    cout << "set fields" << endl;

    int runtime = 10000;
    int every = 100;

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

            strep1 << eps;
            strep2 << eps1;
            strep3 << eps2;

            string s1 = "eps=" + strep1.str() + "_eps1=" + strep2.str() + "_eps2=" + strep3.str();
            stringstream ss;
            ss << setw(number_of_digits) << setfill('0') << i / every;
            string s2 = "_i=" + ss.str();

            a.print_all_results(s1 + s2);
        }
        cout << i << endl;
        a.Update();
    }

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