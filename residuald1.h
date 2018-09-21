#ifndef residuald1_h
#define residuald1_h

#include<vector>
#include<Eigen/Sparse>
#include"structures.h"

using namespace Eigen;

SparseMatrix<double> evaldRdW(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);
MatrixXd evaldRdS(
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);

SparseMatrix<double> evaldRdW_FD(
    std::vector<double> W,
    std::vector<double> area);

void StegerJac(
    std::vector<double> W,
    std::vector<double> &Ap_list,
    std::vector<double> &An_list,
    std::vector<double> &Flux);

std::vector<double> evaldlambdadW(
    std::vector<double> W,
    int i);

void ScalarJac(
    std::vector<double> W,
    std::vector<double> &Ap_list,
    std::vector<double> &An_list);

std::vector<double> evaldpdW(
    std::vector<double> W,
    std::vector<double> area);

MatrixXd evaldRdS_FD(
    std::vector<double> Flux,
    std::vector<double> area,
    std::vector<double> W);

#endif
