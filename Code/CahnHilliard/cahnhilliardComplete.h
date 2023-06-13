#ifndef CAHNHILLIARDCOMPLETE_H
#define CAHNHILLIARDCOMPLETE_H

#include "cahnhilliard.h"
#include </home/dino/Documents/IsingPolymer/Eigen/eigen-3.4.0/Eigen/Sparse>

struct CHCF : public CHC
{
    CHCF(const CH_builder &p);

    Eigen::SparseMatrix<complex<double>> CalculateJ(complex<double> **input, bool); // based on the current state calculate J
    Eigen::SparseMatrix<complex<double>> CalculateJ_2PS(complex<double> **input, bool); // based on the current state calculate J

    int M;                        // memory length
    vector<Eigen::VectorXcd> OLD; // old weights

    Eigen::VectorXcd CalculateWeight(complex<double> **input);
    Eigen::VectorXcd CalculateWeight_2PS(complex<double> **input);
    Eigen::VectorXcd CalculateWeightPartial(complex<double> *input);
    Eigen::VectorXcd CalculateRHS(complex<double> **, complex<double> **, complex<double> **);
    Eigen::VectorXcd CalculateRHS_2PS(complex<double> **, complex<double> **, complex<double> **);
    Eigen::VectorXd CalculateRHS_real(complex<double> **, complex<double> **, complex<double> **);
    Eigen::VectorXcd SolveLinearProblem(Eigen::SparseMatrix<complex<double>> &unchangedJ, complex<double> **, complex<double> **);

    void SetupFracScheme(int MM);
    void SetupFracScheme_2PS(int MM);
    void UpdateWithNewton(bool);
    void UpdateWithNewton_2PS(bool);
    void UpdateWithNewtonCalcJ();
    Eigen::SparseMatrix<double> PartialJCalculation(Rule_Wrapper<cd, cd, cd, cd> &);
    Eigen::SparseMatrix<double> CalculateInitialJ();
    void UpdateWithNewtonGivenJ(Eigen::SparseMatrix<double> &);
};

#include "cahnhilliardComplete.cpp"

#endif