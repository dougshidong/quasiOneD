#ifndef adjoint_h
#define adjoint_h

#include<vector>
#include<Eigen/Eigen>

using namespace Eigen;


VectorXd adjoint(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> S,
    std::vector<double> W,
    std::vector<double> designVar,
    VectorXd &psi);

VectorXd evalpsidRdS(
    VectorXd psiV,
    std::vector<double> Flux,
    std::vector<double> p);

VectorXd buildbMatrix(std::vector<double> dIcdW);

VectorXd evaldIcdW(
    std::vector<double> W,
    std::vector<double> S);

MatrixXd solveSparseAXB(
    SparseMatrix<double> A,
    MatrixXd b,
    int eig_solv);

VectorXd itSolve(
    SparseMatrix<double> A,
    VectorXd b);

#endif
