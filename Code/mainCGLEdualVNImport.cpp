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
    string importstring;
    double c1,c3,c4;
    // int dir_no;
    if (argc == 3)
    {
        c1 = atof(argv[1]);
        c3 = atof(argv[2]);
    }
    else
    {
        error("no");
    }

    string ps = "";
    stringstream ss;
    ss << c1;
    string s1  =ss.str();

    stringstream ss2;
    ss2 << c3;
    string s2 = ss2.str();

    ps += "c1=";
    ps += s1;
    ps += "_c3=";
    ps += s2;



    // double T;
    // bool err1;
    // matrix<double> mat1 = importcsv(importstring, T, err1);


    int N1 = 512;
    int N2 = 512;
    complex<double> *upd1;
    complex<double> *upd2;
    // double c1,c2,c3,c4;
    // double sigma;
    // double L;
    // double deltat;

    // string strbase = "./";

    // for(int ih = dir_no  ; ih < dir_no+50 ; ih++) {

    // string dirs="dir";
    // stringstream ssg;
    // ssg << setw(2) << setfill('0') << ih;

    // dirs += ssg.str();

    // string dir_to_make = strbase + dirs;

    // int status = mkdir(dir_to_make.c_str(),0777);

    double r1 = (double)rand() / (double)RAND_MAX;
    double r2 = (double)rand() / (double)RAND_MAX;
    double r3 = (double)rand() / (double)RAND_MAX;
    double r4 = (double)rand() / (double)RAND_MAX;
    double r5 = (double)rand() / (double)RAND_MAX;

    // double c1 = 2 * (2 * r1 - 1);
    // double c3 = 2 * (2 * r3 - 1);
    double sigma  = r5;
    // double c1 = 0.8829490258;
    // double c2 = -0.543766673;
    // double c3 = -0.1703678305;
    // double c4 = 0.8733192491;
    // double sigma = 0.7543335486;

    double L = 200.;
    double deltat = 0.0005;

    // string params = "param";
    // vector1<double> pa(2);

    // pa[0] = c1;
    // pa[1] = c3;

    // outfunc(pa,params);

    upd1 = (complex<double> *)fftw_malloc(N1 * N2 * sizeof(complex<double>));
    upd2 = (complex<double> *)fftw_malloc(N1 * N2 * sizeof(complex<double>));

    c4 = 1;

    // -fac3*(|B|^2)*A
    // -fac7*(|A|^2)*B

    complex<double> fac1(1., c1);

    complex<double> fac2(1., c3);

    complex<double> fac3(-1., c4);


    double r;
    
    double gamma5,delta5;

    double gamma6,delta6;

    double gamma7,delta7;


    r = 1;
    gamma5=1.;
    delta5=0.;

    gamma6=1.;
    delta6= -0.588235;

    gamma7=-1.;
    delta7=1.;

    complex<double> fac4(r, 0.);

    complex<double> fac5(gamma5,delta5);

    complex<double> fac6(gamma6, delta6);

    complex<double> fac7(gamma7, delta7);

    for (int i1 = 0; i1 < N1; i1++)
    {
        for (int j = 0; j < N2; j++)
        {
            double k1, k2;
            k1 = i1;
            k2 = j;
            // if (i1 <= N1 / 2)
            // {
            //     k1 = i1;
            // }
            // else
            // {
            //     k1 = (i1 - N1);
            // }
            // if (j <= N2 / 2)
            // {
            //     k2 = j;
            // }
            // else
            // {
            //     k2 = (j - N2);
            // }

            upd1[i1 * N2 + j] = (1. / (1. + fac1*deltat * (2 * pii / L) * (2 * pii / L)*(SQR(k1)+SQR(k2))));

            upd2[i1 * N2 + j] = (1. / (1. + fac5*deltat * (2 * pii / L) * (2 * pii / L) * (SQR(k1) + SQR(k2))));

            // double tempor = SQR(k1) + SQR(k2);

            // upd2[i * params.N2 + j] = 1. / (1. + dt * Di * temp1 * tempor);
        }
    }

    CH_builder myp;
    myp.number_of_fields=12;
    myp.dimension=2;
    myp.N1 = N1;
    myp.N2 = N2;
    myp.N3 = 0;


    CH_builder myp2;
    myp2.number_of_fields = 4;
    myp2.dimension = 2;
    myp2.N1 = N1;
    myp2.N2 = N2;
    myp2.N3 = 0;

    Field_Wrapper<double,double> store(myp);
    Field_Wrapper<double,double > transformed1(myp);

    // Field_Wrapper<complex<double>, complex<double>> transformed2(myp);

    // Field_Wrapper<complex<double>, complex<double>> transformed3(myp);
    Field_Wrapper<double,double> store2(myp2);
    Field_Wrapper<double,double> reverse_transform(myp2);

    // complex<double>* initial_field;

    // initial_field = (complex<double> *)fftw_malloc(N1 * N2 * sizeof(complex<double>));
    // matrix<complex<double> > initial_field(N1,N2);
    // for(int i = 0  ; i < N1 ; i++) {
    //     for(int j  = 0 ; j < N2 ; j++) {
    //         double a=5*(2.*(double)rand()/(double)RAND_MAX-1);
    //         double b = 5*(2. * (double)rand() / (double)RAND_MAX - 1);
    //         initial_field(i,j) = complex<double>(a,b);
    //     }
    // }


    double r0 = 20.;
    double theta0 = pii/4;
    double A0=1;
    double rc=10.;

    // matrix<complex<double>> initial_field(N1, N2);
    // matrix<complex<double>> initial_field2(N1, N2);

    bool err1,err2,err3,err4;
    double T;
    matrix<double> r1r = importcsv("f1r.csv", T, err1);
    matrix<double> r1i = importcsv("f1i.csv", T, err2);
    matrix<double> r2r = importcsv("f2r.csv", T, err3);
    matrix<double> r2i = importcsv("f2i.csv", T, err4);
    /*
    for (int i = 0; i < N1; i++)
    {
        for (int j = 0; j < N2; j++)
        {
            double x = (i - 0.25 * N1);
            double y = (j - 0.25 * N2);
            double r =  sqrt(SQR(x)+SQR(y));
            double theta = atan2(y,x);
            if(r<r0)
                initial_field(i, j) = { A0 * tanh(r / rc) * cos(theta),
                                        A0 * tanh(r / rc) * sin(theta)};
            else {
                    initial_field(i, j) = { A0 * cos(theta0),
                                            A0 * sin(theta0) };
            }
                                        // double a = 5 * (2. * (double)rand() / (double)RAND_MAX - 1);
                                        // double b = 5 * (2. * (double)rand() / (double)RAND_MAX - 1);
                                        // initial_field(i, j) = complex<double>(a, b);
        }
    }
    */
    matrix<complex<double>> initial_field = combinematrix(r1r, r1i);
    matrix<complex<double>> initial_field2 = combinematrix(r2r, r2i);

    CosineWeightForward fw;
    transformed1.add_method(fw, 0);
    transformed1.add_method(fw, 1);
    transformed1.add_method(fw, 2);
    transformed1.add_method(fw, 3);
    transformed1.add_method(fw, 4);
    transformed1.add_method(fw, 5);
    transformed1.add_method(fw, 6);
    transformed1.add_method(fw, 7);
    transformed1.add_method(fw, 8);
    transformed1.add_method(fw, 9);
    transformed1.add_method(fw, 10);
    transformed1.add_method(fw, 11);
    CosineWeightBackward fw2;
    reverse_transform.add_method(fw2, 0);
    reverse_transform.add_method(fw2, 1);
    reverse_transform.add_method(fw2, 2);
    reverse_transform.add_method(fw2, 3);

    int iter= 0;
    int every =  100;
    int number_of_digits=5;

    // using clock = std::chrono::high_resolution_clock;
    // using duration_t = clock::duration; // native tick type

    // duration_t total1{0};
    // duration_t total2{0};
    // duration_t total3{0};
    // duration_t total4{0};

    matrix<complex<double>> initial_field3(N1, N2);
    matrix<complex<double>> initial_field4(N1, N2);
    matrix<complex<double>> initial_field5(N1, N2);
    matrix<complex<double>> initial_field6(N1, N2);
    matrix<complex<double>> res(N1, N2);
    matrix<complex<double>> res2(N1, N2);

    for(;;) {
        cout << iter << endl;
        if (iter % every == 0 && iter > 000000 )
        {
            matrix<double> if_real(N1, N2);
            matrix<double> if_imag(N1, N2);
            matrix<double> if_real2(N1, N2);
            matrix<double> if_imag2(N1, N2);

            for (int i = 0; i < N1; i++)
            {
                for (int j = 0; j < N2; j++)
                {
                    if_real(i, j) = sqrt(SQR(initial_field(i, j).real()) + SQR(initial_field(i, j).imag()));
                    if_imag(i, j) = atan2(initial_field(i,j).imag(),initial_field(i,j).real());

                    if_real2(i, j) = sqrt(SQR(initial_field2(i, j).real()) + SQR(initial_field2(i, j).imag()));
                    if_imag2(i, j) = atan2(initial_field2(i, j).imag(), initial_field2(i, j).real());
                }
            }

            // double max1;
            // double max2;
            // if_real.maxima(max1);
            // if_imag.maxima(max2);
            // cout << iter << endl;

            // cout << max1 << " " << max2 << endl;

            stringstream ii;
            ii << setw(number_of_digits) << setfill('0') << iter / every;
            string inx = "_i=";
            outfunc(if_real,"res1r" + ps + inx + ii.str());
            outfunc(if_imag, "res1i" + ps + inx + ii.str());
            outfunc(if_real2, "res2r" + ps + inx + ii.str());
            outfunc(if_imag2, "res2i" + ps + inx + ii.str());
        }
    
    // auto start1=clock::now();

    for(int i = 0  ; i < N1 ; i++) {
        for(int j  = 0 ; j < N2 ; j++) {
            double aif = abs(initial_field(i, j));
            double aif2 = abs(initial_field2(i, j));
            initial_field3(i, j) = SQR(aif) * initial_field(i, j);
            initial_field4(i, j) = SQR(aif2) * initial_field(i, j);
            initial_field5(i, j) = SQR(aif2) * initial_field2(i, j);
            initial_field6(i, j) = SQR(aif) * initial_field2(i, j);
        }
    }

    store.set_field(realmatrix(initial_field), 0);
    store.set_field(imaginarymatrix(initial_field), 1);
    store.set_field(realmatrix(initial_field3), 2);
    store.set_field(imaginarymatrix(initial_field3), 3);
    store.set_field(realmatrix(initial_field4), 4);
    store.set_field(imaginarymatrix(initial_field4), 5);
    store.set_field(realmatrix(initial_field2), 6);
    store.set_field(imaginarymatrix(initial_field2), 7);
    store.set_field(realmatrix(initial_field5), 8);
    store.set_field(imaginarymatrix(initial_field5), 9);
    store.set_field(realmatrix(initial_field6), 10);
    store.set_field(imaginarymatrix(initial_field6), 11);
    // total1 += clock::now() - start1;

    // auto start2=clock::now();

    transformed1.Calculate_Results(store.calculated_reactions);
    // total2 += clock::now() - start2;

    // auto start3 = clock::now();

    //#pragma omp simd
    for (int i1 = 0; i1 < N1; i1++)
    {
        for (int j = 0; j < N2; j++)
        {
            int indx =i1 * N2 + j;
            res(i1, j) = upd1[indx] * (1. + deltat) * (complex<double>(transformed1.calculated_reactions[0][indx], transformed1.calculated_reactions[1][indx])) - (upd1[indx]) * deltat * fac2 * (complex<double>(transformed1.calculated_reactions[2][indx], transformed1.calculated_reactions[3][indx])) - (upd1[indx]) * deltat * fac3 * (complex<double>(transformed1.calculated_reactions[4][indx], transformed1.calculated_reactions[5][indx]));

            res2(i1, j) = upd2[indx] * (fac4 + deltat) * (complex<double>(transformed1.calculated_reactions[6][indx], transformed1.calculated_reactions[7][indx])) - (upd2[indx]) * deltat * fac6 * (complex<double>(transformed1.calculated_reactions[8][indx], transformed1.calculated_reactions[9][indx])) - (upd2[indx]) * deltat * fac7 * (complex<double>(transformed1.calculated_reactions[10][indx], transformed1.calculated_reactions[11][indx]));

            // cout << i1 << " " << j << endl;
            // cout << upd1[indx] << endl;
            // cout << transformed1.calculated_reactions[0][indx] << endl;
            // cout << (upd2[indx]) * transformed1.calculated_reactions[1][indx] << endl;
            // cout << deltat * fac4 * transformed1.calculated_reactions[2][indx] << endl;
            // cout << res(i1,j) << endl;
            // pausel();
        }
    }

    //
    store2.set_field(realmatrix(res), 0);
    store2.set_field(imaginarymatrix(res), 1);
    store2.set_field(realmatrix(res2), 2);
    store2.set_field(imaginarymatrix(res2), 3);
    // total3 += clock::now() - start3;
    

    // auto start4 = clock::now();
    reverse_transform.Calculate_Results(store2.calculated_reactions);


    for(int i = 0  ; i < N1 ; i++) {
        for(int j  = 0 ; j < N2 ; j++) {
            initial_field(i, j) = complex<double>(reverse_transform.calculated_reactions[0][i * N2 + j],reverse_transform.calculated_reactions[1][i * N2 + j]);
            initial_field2(i, j) = complex<double>(reverse_transform.calculated_reactions[2][i * N2 + j], reverse_transform.calculated_reactions[3][i * N2 + j]);
        }
    }
    // total4 += clock::now() - start4;

    iter++;

    if(iter>1000000) break;
   // pausel();
    //transformed1.Calculate_Results()
       
    //    std::cout << "total = " << total1.count() / 10E6 << " µs\n";
    //    std::cout << "total = " << total2.count() / 10E6 << " µs\n";
    //    std::cout << "total = " << total3.count() / 10E6 << " µs\n";
    //    std::cout << "total = " << total4.count() / 10E6 << " µs\n";
    }
    // auto total_us1 = std::chrono::duration_cast<std::chrono::microseconds>(total1);
    // auto total_us2 = std::chrono::duration_cast<std::chrono::microseconds>(total2);
    // auto total_us3 = std::chrono::duration_cast<std::chrono::microseconds>(total3);
    // auto total_us4 = std::chrono::duration_cast<std::chrono::microseconds>(total4);


    // ofstream myfile;
    // ofstream myfile2;
    // string fil1 = string("field0") + importstring;
    // string fil2 = string("field1") + importstring;
    // myfile.open(fil1.c_str());
    // myfile2.open(fil2.c_str());
    // for (int i = 0; i < runtime; i++)
    // {

    //     if (i % every == 0 && i > 00000)
    //     {

    //         for(int j = 0 ; j < p.N1-1 ; j++) {
    //             // cout << a.fields[0][j * p.N1] << endl;
    //             myfile << a.fields[0][j*p.N1] <<",";
    //         }
    //         myfile << a.fields[0][p.N1*p.N1-1] << endl;

    //         for (int j = 0; j < p.N1 - 1; j++) {
    //             myfile2 << a.fields[1][j * p.N1] << ",";
    //         }
    //         myfile2 << a.fields[1][p.N1 * p.N1 - 1] << endl;

    //         // pausel();
    //         // stringstream strep1;
    //         // stringstream strep2;
    //         // stringstream strep3;
    //         // stringstream strep4;

    //         // strep1 << dens;
    //         // strep4 << c0;
    //         // strep2 << c1;
    //         // strep3 <<  surf;
    //         // string s1 = importstring;
    //         // // // string s1 = "denp=" + strep1.str() + "c0=" + strep4.str() + "_c1=" + strep2.str() + "_surf=" + strep3.str();
    //         // stringstream ss;
    //         // ss << setw(number_of_digits) << setfill('0') << i / every;
    //         // string s2 = "_i=" + ss.str();
    //         // string su = s1.substr(0, s1.size() - 4);
    //         // cout << su + s2 << endl;
    //         // a.print_some_results(su + s2,ps2);
    //     }
    //     cout << i << endl;
    //     cout << "begin" << endl;
    //     a.Update();
    //     bool chck = true;
    //     a.check_field(chck);
    //     if (!chck)
    //         break;
    // }
    
}