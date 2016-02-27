#ifndef adjoint_h
#define adjoint_h

#include<vector>
#include<Eigen/Eigen>

using namespace Eigen;


VectorXd adjoint(
    std::vector <double> x, 
    std::vector <double> dx, 
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> &psi,
    std::vector <double> designVar);

VectorXd evalpsidRdS(
    VectorXd psiV,
    std::vector <double> Flux,
    std::vector <double> p);

MatrixXd evaldRdS(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W);

MatrixXd evaldRdS_FD(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W);

void JacobianCenter(
    std::vector <double> &J,
    double u, double c);

SparseMatrix<double> evaldRdW(
    std::vector <double> Ap,
    std::vector <double> An,
    std::vector <double> W,
    std::vector <double> dQdW,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> S,
    double Min);

SparseMatrix<double> evaldRdW_FD(
    std::vector <double> W,
    std::vector <double> S,
    double Min);

VectorXd buildbMatrix(std::vector <double> dIcdW);

void StegerJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list,
    std::vector <double> &Flux);

void ScalarJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list);

void BCJac(
    std::vector <double> W,
    std::vector <double> dt,
    std::vector <double> dx,
    std::vector <double> &dBidWi,
    std::vector <double> &dBidWd,
    std::vector <double> &dBodWd,
    std::vector <double> &dBodWo);

VectorXd evaldIcdW(std::vector <double> W, std::vector <double> S);

MatrixXd evaldSdDesign(
    std::vector <double> x, 
    std::vector <double> dx, 
    std::vector <double> designVar);

void evaldQdW(
    std::vector <double> &dQdW,
    std::vector <double> W,
    std::vector <double> S);

MatrixXd solveSparseAXB(SparseMatrix <double> A, MatrixXd b, int eig_solv);
VectorXd itSolve(SparseMatrix <double> A, VectorXd b);

#endif
