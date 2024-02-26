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


int main(int argc, char **argv)
{

    srand(time(NULL));
    string importstring;
    CH_builder p;
    p.dimension = 1;
    p.number_of_fields = 4;
    p.N1 = 256;

    Field_Wrapper<double,double> a(p);
    Field_Wrapper<double, double> a2(p);

    double T;
    bool err1, err2, err3, err4;
    double *f1 = new double[p.N1];
    double *f2 = new double[p.N1];
    double *f3 = new double[p.N1];
    double *f4 = new double[p.N1];

    for (int i = 0; i < p.N1; i++)
    {
        f1[i] = rand() / (double)(RAND_MAX);
        f2[i] = rand() / (double)(RAND_MAX);
        f3[i] = rand() / (double)(RAND_MAX);
        f4[i] = rand() / (double)(RAND_MAX);
    }

    CH<double> b(p);

    b.set_field(f1, 0);
    b.set_field(f2, 1);
    b.set_field(f3, 2);
    b.set_field(f4, 3);

    for (int i = 0; i < p.number_of_fields; i++)
    {
        CosineWeightForward1D fw;
        a.add_method(fw, i);
        CosineWeightBackward1D fw2;
        a2.add_method(fw2, i);
    }

    cout << b.fields[0][0] << endl;
    a.Calculate_Results(b.fields);

    cout << a.calculated_reactions[0][0] << endl;
    a2.Calculate_Results(a.calculated_reactions);

    cout << a2.calculated_reactions[0][0] << endl;

    pausel();
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