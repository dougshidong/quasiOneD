#ifndef residuald1_h
#define residuald1_h

#include<vector>
#include<Eigen/Sparse>

using namespace Eigen;

SparseMatrix<double> evaldRdW(
    std::vector <double> W,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> S);

SparseMatrix<double> evaldRdW_FD(
    std::vector <double> W,
    std::vector <double> S);

void StegerJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list,
    std::vector <double> &Flux);

std::vector <double> evaldlambdadW(
    std::vector <double> W,
    int i);

void ScalarJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list);

std::vector <double> evaldpdW(
    std::vector <double> W,
    std::vector <double> S);

MatrixXd evaldRdS(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W);

MatrixXd evaldRdS_FD(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W);

#endif
