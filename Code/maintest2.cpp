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
#include "ParticleReactionDiffusion/cellarray.h"


int main(int argc, char **argv)
{

    srand(time(NULL));




    int Nt = 1;
    int Nr = 2;

    cell c(Nt,Nr);

    Subdiffusion sd1;

    sd1.alpha=0.5;
    sd1.tau=1;
    sd1.p = 100.;
    sd1.normal = false;

    c.set_subdiffusion(sd1, 0);
    // c.set_subdiffusion(sd1, 1);
    // c.set_subdiffusion(sd1, 2);

    // chemical_reaction cr1,cr2;

    // cr1.orderin = 2;
    // cr2.orderout = 1;
    // cr1.i1 = 0;
    // cr1.i2 = 1;
    // cr1.o1 = 2;
    // cr1.rate = 1.0;

    // cr2.orderin = 1;
    // cr2.orderout = 2;
    // cr2.i1 = 2;
    // cr2.o1 = 0;
    // cr2.o2 = 1;
    // cr2.rate = 1.0;

    // c.set_chemical_reaction(cr1, 0);
    // c.set_chemical_reaction(cr2, 1);

    // vector1<int> nu(3);
    // nu[0] = 10;
    // nu[1] = 10;
    // nu[2] = 10;
    // c.initialize(nu);

    int Nx;
    double eps1;
    int neq1;
    double eps2;
    int nav;
    if(argc == 6) {
        Nx = atof(argv[1]);
        eps1 = atof(argv[2]);
        neq1 = atof(argv[3]);
        eps2 = atof(argv[4]);
        nav = atof(argv[5]);
    }
    string str1 = "density";
    stringstream ss1,ss2,ss3,ss4,ss5;
    ss1 << Nx;
    ss2 << eps1;
    ss3 << neq1;
    ss4 << eps2;
    ss5 << nav;

    string s1 = "_N=" + ss1.str();
    string s2 = "_eps1=" + ss2.str();
    string s3 = "_neq=" + ss3.str();
    string s4 = "_eps2=" + ss4.str();
    string s5 = "_nav=" + ss5.str();

    str1 += s1;
    str1 += s2;
    str1 += s3;
    str1 += s4;
    str1 += s5;
    cell_array full(Nx,1,c);

    attractive_self_interaction attpot;
    attpot.eps=eps1;
    attpot.Neq = neq1;

    linear_surface_interaction surfpot;
    surfpot.eps = eps2;

    full.set_si_interaction(attpot,0);
    full.set_cnni_interaction(surfpot, 0, 0);

    for(int i  = 0 ; i < Nx ; i++) {
        for(int j = 0 ; j < Nx ; j++) {
            vector1<int> nu(1);

            //int r1 = -5 + rand() % 11;

            nu[0] = nav;

            full.reactors(i,j).initialize(nu);
            // full.reactors(i,j).generate_diffusion_time();
            // full.reactors(i,j).generate_time(0.0);
        }
    }

    
    full.initialize_next_diffusion_events(); //construct the sorted list
    cout << full.next_diffusion_events[0].a << endl;
    
    //full.iterate();
    // full.find_next_reaction_event();
    // full.find_next_diffusion_event();
    int nmax = nav*Nx*Nx*10000;
    int iter=0;
    int sch = (nmax/1000);
    vector1<int> counter(3);
    ofstream file1;
    file1.open("times.csv");
    for (;;)
    {
        full.do_a_diffusion_event(counter);
        iter++;
        if(iter > 0&& iter % sch == 0) {
            cout << full.current_time << endl;
            matrix<int> A(Nx,Nx);

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Nx; j++)
                {
                    A(i, j) = full.reactors(i, j).typenumber[0];
                    // B(i, j) = full.reactors(i, j).typenumber[1];
                    // C(i, j) = full.reactors(i, j).typenumber[2];
                }
            }
            stringstream ss;
            ss << setw(4);
            ss << std::setfill('0');
            ss << iter/sch;
            string p = str1+"_i="+ss.str(); 
            outfunc(A, p);
            file1 << full.current_time << endl;

            full.reset_all_times();
        }
        if (iter > nmax)
            break;
    }
    file1.close();

    // matrix<int> A(100, 100);
    // matrix<int> B(100, 100);
    // matrix<int> C(100, 100);

    // for (int i = 0; i < 100; i++)
    // {
    //     for (int j = 0; j < 100; j++)
    //     {
    //         A(i, j) = full.reactors(i, j).typenumber[0];
    //         // B(i, j) = full.reactors(i, j).typenumber[1];
    //         // C(i, j) = full.reactors(i, j).typenumber[2];
    //     }
    // }
    // outfunc(A, "A");
    // outfunc(B, "B");
    // outfunc(C, "C");

    /*
    int iter=0;
    for(;;) {
     full.iterate();
     iter++;
     cout << iter << endl;
     if(iter>10000) break;
    }
    
    */

    // c.generate_time(0.0);
    // c.generate_diffusion_time();
    // for(int j = 0 ; j < 1000 ; j++) {
    // for(int i = 0 ; i < Nt ; i++)
    // cout << c.diffusion_times[i].size() << ",";
    // cout << endl;
    // c.do_reaction(0.0);
    // }

    // string importstring;
    // CH_builder p;
    // p.dimension = 1;
    // p.number_of_fields = 4;
    // p.N1 = 256;

    // Field_Wrapper<double,double> a(p);
    // Field_Wrapper<double, double> a2(p);

    // double T;
    // bool err1, err2, err3, err4;
    // double *f1 = new double[p.N1];
    // double *f2 = new double[p.N1];
    // double *f3 = new double[p.N1];
    // double *f4 = new double[p.N1];

    // for (int i = 0; i < p.N1; i++)
    // {
    //     f1[i] = rand() / (double)(RAND_MAX);
    //     f2[i] = rand() / (double)(RAND_MAX);
    //     f3[i] = rand() / (double)(RAND_MAX);
    //     f4[i] = rand() / (double)(RAND_MAX);
    // }

    // CH<double> b(p);

    // b.set_field(f1, 0);
    // b.set_field(f2, 1);
    // b.set_field(f3, 2);
    // b.set_field(f4, 3);

    // for (int i = 0; i < p.number_of_fields; i++)
    // {
    //     CosineWeightForward1D fw;
    //     a.add_method(fw, i);
    //     CosineWeightBackward1D fw2;
    //     a2.add_method(fw2, i);
    // }

    // cout << b.fields[0][0] << endl;
    // a.Calculate_Results(b.fields);

    // cout << a.calculated_reactions[0][0] << endl;
    // a2.Calculate_Results(a.calculated_reactions);

    // cout << a2.calculated_reactions[0][0] << endl;

    // pausel();
    // vector<int> OLD;
    // for(int i = 0  ; i < 10 ; i++)
    // OLD.push_back(i);

    // OLD.pop_back();
    // OLD.push_back(10);
    // std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());

    // OLD.pop_back();
    // OLD.push_back(11);
    // std::rotate(OLD.rbegin(), OLD.rbegin() + 1, OLD.rend());

    // for(int i =0 ; i < OLD.size() ; i++)
    // cout << OLD[i] <<",";

    // cout << endl;

    // int n = 128;   // size of the image
    // int m = n * n; // number of unknowns (=number of pixels)

    // // Assembly:
    // std::vector<T> coefficients; // list of non-zeros coefficients
    // Eigen::VectorXd b(m);        // the right hand side-vector resulting from the constraints
    // buildProblem(coefficients, b, n);

    // SpMat A(m, m);
    // A.setFromTriplets(coefficients.begin(), coefficients.end());

    // // Solving:
    // // auto start = std::chrono::high_resolution_clock::now();
    // // Eigen::SimplicialCholesky<SpMat> chol(A); // performs a Cholesky factorization of A
    // // auto stop = std::chrono::high_resolution_clock::now();
    // // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // // auto start2 = std::chrono::high_resolution_clock::now();
    // // Eigen::VectorXd x = chol.solve(b);        // use the factorization to solve for the given right hand side
    // // auto stop2 = std::chrono::high_resolution_clock::now();
    // // auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);

    // // cout << "Time taken by function: "
    // //      << duration.count() << " microseconds" << endl;

    // // cout << "Time taken by function: "
    // //      << duration2.count() << " microseconds" << endl;

    // Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
    // auto start = std::chrono::high_resolution_clock::now();
    // cg.compute(A);
    // Eigen::VectorXd x = cg.solve(b);
    // auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // cout << "Time taken by function: "
    //      << duration.count() << " microseconds" << endl;
    // std::cout << "#iterations:     " << cg.iterations() << std::endl;
    // std::cout << "estimated error: " << cg.error() << std::endl;

    // update b, and solve again
    // x = cg.solve(b);

    // Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver;
    // // fill A and b;
    // // Compute the ordering permutation vector from the structural pattern of A
    // solver.analyzePattern(A);
    // // Compute the numerical factorization
    // solver.factorize(A);
    // // Use the factors to solve the linear system
    // Eigen::VectorXd x = solver.solve(b);
    // cout << x << endl;
    // Export the result to a file:
    // saveAsBitmap(x, n, argv[1]);

    return 0;
        }