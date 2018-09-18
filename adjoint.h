#ifndef adjoint_h
#define adjoint_h

#include<vector>
#include<Eigen/Eigen>

using namespace Eigen;

VectorXd adjoint(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> W,
    std::vector<double> designVar,
    VectorXd &psi);

MatrixXd solveSparseAXB(
    SparseMatrix<double> A,
    MatrixXd b,
    int eig_solv);

#endif
