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

#include </home/dino/Documents/IsingPolymer/Eigen/eigen-3.4.0/Eigen/Sparse>


typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

// void buildProblem(std::vector<T> &coefficients, Eigen::VectorXd &b, int n);
// void saveAsBitmap(const Eigen::VectorXd &x, int n, const char *filename);

void insertCoefficient(int id, int i, int j, double w, std::vector<T> &coeffs,
                       Eigen::VectorXd &b, const Eigen::VectorXd &boundary)
{
    int n = int(boundary.size());
    int id1 = i + j * n;

    if (i == -1 || i == n)
        b(id) -= w * boundary(j); // constrained coefficient
    else if (j == -1 || j == n)
        b(id) -= w * boundary(i); // constrained coefficient
    else
        coeffs.push_back(T(id, id1, w)); // unknown coefficient
}

void buildProblem(std::vector<T> &coefficients, Eigen::VectorXd &b, int n)
{
    b.setZero();
    Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0, M_PI).sin().pow(2);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            int id = i + j * n;
            insertCoefficient(id, i - 1, j, -1, coefficients, b, boundary);
            insertCoefficient(id, i + 1, j, -1, coefficients, b, boundary);
            insertCoefficient(id, i, j - 1, -1, coefficients, b, boundary);
            insertCoefficient(id, i, j + 1, -1, coefficients, b, boundary);
            insertCoefficient(id, i, j, 4, coefficients, b, boundary);
        }
    }
}



using namespace std;
typedef complex<double> cdf;
void test_function(complex<double> **fields) {
    CH_builder p;
    int nof = 5;
    p.number_of_fields = nof;
    p.N1 = 64;
    p.N2 = 64;

    // Field_Wrapper<complex<double>, complex<double>> tempf(p);
    cdf** calculated_reactions = new cdf *[p.number_of_fields];

    for (int i = 0; i < p.number_of_fields; i++)
    {
        calculated_reactions[i] = (cdf *)fftw_malloc(p.N1 * p.N2 * sizeof(cdf));
        for (int j = 0; j < p.N1 * p.N2; j++)
        {
            cdf b = 0.0;
            calculated_reactions[i][j] = b;
        }
    }

    cdf **fieldscopy = new cdf *[p.number_of_fields];

    for (int i = 0; i < p.number_of_fields; i++)
    {
        calculated_reactions[i] = (cdf *)fftw_malloc(p.N1 * p.N2 * sizeof(cdf));
        for (int j = 0; j < p.N1 * p.N2; j++)
        {
            cdf b = 0.0;
            fieldscopy[i][j] = fields[i][j];
        }
    }

    for(int j = 0  ; j < p.number_of_fields ; j++) {
    fftw_plan p2(fftw_plan_dft_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *>(fieldscopy[j]), reinterpret_cast<fftw_complex *>(calculated_reactions[j]), FFTW_FORWARD, FFTW_ESTIMATE));
    fftw_execute(p2);
    double corr = 1. / (p.N1);
    int end = p.get_total();
    for (int i = 0; i < end; i++)
    {
        calculated_reactions[j][i] *= corr;
    }
    fftw_destroy_plan(p2);

    }
    cout << "done" << endl;

    for (int i = 0; i < p.number_of_fields; i++)
    {
        fftw_free(calculated_reactions[i]);
        fftw_free(fieldscopy[i]);
    }
    delete[] calculated_reactions;
    delete[] fieldscopy;

    //pausel();


}

int main(int argc, char **argv)
{

    srand(time(NULL));
    int fftw_init_threads(void);
    string importstring;

    if (argc == 2)
    {
        stringstream ss;
        ss << argv[1];
        importstring = ss.str();
    }
    else
    {
        error("no");
    }

    double T;
    bool err1;
    matrix<double> mat1 = importcsv(importstring, T, err1);

    int n = mat1(0, 0);

    vector1<double> phasesepsparams(8);

    for (int i = 0; i < 8; i++)
        phasesepsparams[i] = mat1(1, i);

    vector1<double> epsi((n) * (n - 1) / 2);

    for (int i = 0; i < (n) * (n - 1) / 2; i++)
    {
        epsi[i] = -mat1(2, i);
    }

    int k = 0;
    matrix<double> epsa(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            epsa(i, j) = epsi[k];
            epsa(j, i) = epsi[k];
            k++;
        }
    }

    vector1<double> init(n);

    for (int i = 0; i < n; i++)
        init[i] = mat1(3, i);

    // if (argc == 5)
    // {
    //     c0 = atof(argv[1]);
    //     c1 = atof(argv[2]);
    //     eps = atof(argv[3]);
    //     dens = atof(argv[4]);
    // }
    // else {
    //     c0 = 0.2;
    //     c1 = 0.8;
    //     eps = 0.2;
    //     dens = 0.5;
    // }

    typedef complex<double> myc;
    typedef InvasionLinearReversibleB<myc> ILB;
    typedef InvasionLinearReversibleA<myc> ILA;
    typedef Field_Wrapper<myc, myc> FWCC;
    typedef CoupledPhaseSeparatingSystem<myc> CPSS;
    typedef Rule_Wrapper<myc, myc, myc, myc> RWC;

    CH_builder p;
    int nof = n;
    p.number_of_fields = nof;
    p.N1 = 64;
    p.N2 = 64;

    CHC a(p);

    for (int i = 0; i < nof; i++)
    {
        for (int j = i + 1; j < nof; j++)
        {
            // cout << i << " " << j << endl;
            a.set_interaction(epsa(i, j), i, j);
        }
    }

    a.set_diffusion(phasesepsparams[0]);
    a.set_epsilon(phasesepsparams[1]);
    a.set_c0_c1(phasesepsparams[2], phasesepsparams[3]);
    double L = phasesepsparams[4];
    double temp1 = SQR(2. * pii / L);
    a.set_temp1(temp1);

    a.set_alpha(phasesepsparams[5]);
    a.set_dt(phasesepsparams[6]);

    double rate_multiplier = phasesepsparams[7];

    FWCC my_chemsitry(p);
    NoWeight<myc, myc> nw;
    for (int j = 4; j < mat1.getnrows(); j++)
    {
        int no_chem = mat1(j, 0);
        // cout << no_chem << endl;

        if (no_chem == 0)
        {
            my_chemsitry.add_method(nw, j - 4);
        }
        else
        {
            MultipleReactions<myc> c6(no_chem);

            double tot = 0.0;

            for (int i = 0; i < no_chem; i++)
            {
                // cout << i << endl;
                vector1<int> jpow(nof);
                for (int k = i * (nof + 1) + 2; k < i * (nof + 1) + 2 + nof; k++)
                {

                    jpow[k - (i * (nof + 1) + 2)] = (int)mat1(j, k);
                }
                // cout << mat1(j,i*(nof+1)+1) << endl;
                // cout << jpow << endl;
                // pausel();
                GenericChemistry<myc> c6_0(rate_multiplier * mat1(j, i * (nof + 1) + 1), jpow);
                c6.add_chemical_reaction(c6_0, i);
                double tot1 = 1.0;
                for (int k = 0; k < nof; k++)
                {
                    tot1 *= Power(init[k], jpow[k]);
                    // cout << tot1 << ",";
                }
                // cout << endl;

                tot += mat1(j, i * (nof + 1) + 1) * tot1;

                // cout << endl;
            }
            cout << tot << endl;

            my_chemsitry.add_method(c6, j - 4);
        }
    }
    a.set_chems(my_chemsitry);

    // pausel();

    a.setup_matrices();

    vector<matrix<myc>> v;

    for (int j = 0; j < nof; j++)
    {
        matrix<myc> field1(p.N1, p.N2);
        v.push_back(field1);
    }

    double gt = 0.6;

    for (int lk = 0; lk < nof; lk++)
    {
        double x1 = init[lk];
        for (int i = 0; i < p.N1; i++)
        {
            for (int j = 0; j < p.N2; j++)
            {
                double r1 = (2. * ((double)rand() / (double)RAND_MAX) - 1.);
                v[lk](i, j) = x1 + gt * x1 * r1;
            }
        }
    }

    for (int lk = 0; lk < nof; lk++)
    {
        a.set_field(v[lk], lk);
    }
    GetMinimas(a.fields, p);
    GetMaximas(a.fields, p);

    a.calculate_initial_weight();

    cout << "calc" << endl;

    // Eigen::SparseMatrix<complex<double>> a2= a.CalculateJ();

    cout << "Sparse constructed" << endl;

    cout << "all fields set" << endl;

    a.SetupFracScheme(100);

    fftw_plan p2;
    complex<double> *in;
    complex<double> *out;
    p2 = fftw_plan_dft_2d(p.N1, p.N2, reinterpret_cast<fftw_complex *>(in), reinterpret_cast<fftw_complex *>(out), FFTW_FORWARD, FFTW_ESTIMATE);

    cout << "system set up" << endl;
    typedef complex<double> cd;
    int N = p.N1*p.N2;
    int totN = nof*N;
    #pragma omp parallel
    {
        vector<Trip> coeffs_private;
        coeffs_private.reserve(N);
        #pragma omp for nowait schedule(static)
        for (int index = 0; index < totN; index++)
        {
            int fieldno = floor(index / N);
            int jk = index % N;

            cd **fields1;
            cd **fields2;

            fields1 = new cd *[nof];
            fields2 = new cd *[nof];

            for (int i = 0; i < nof; i++)
            {
                fields1[i] = (cd *)fftw_malloc(N * sizeof(cd));
                fields2[i] = (cd *)fftw_malloc(N * sizeof(cd));
                for (int j = 0; j < N; j++)
                {
                    fields1[i][j] = a.fields[i][j];
                    fields2[i][j] = a.fields[i][j];
                }
            }

            double dc = 0.001;
            fields1[fieldno][jk] += dc;

            fields2[fieldno][jk] -= dc;
            // Eigen::VectorXcd w0 = CalculateRHS_real(fields1, transformed1.calculated_reactions, transformed3.calculated_reactions);
            // Eigen::VectorXcd w1 = CalculateRHS_real(fields2, transformed1.calculated_reactions, transformed3.calculated_reactions);

            //Eigen::VectorXcd w0 = a.CalculateWeight(fields1);
            test_function(fields1);

            for (int i = 0; i < nof; i++)
            {
                fftw_free(fields1[i]);
                fftw_free(fields2[i]);
            }

            delete[] fields1;
            delete[] fields2;

        }

    }
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