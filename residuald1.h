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
MatrixXd evaldRdArea(
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
	const double gam,
    const double rho,
	const double rho_u,
	const double e);

void dFluxdW_scalard(
	const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector<double> &Ap_list,
    std::vector<double> &An_list);

std::vector<double> evaldpdW(
	const double gamma,
    const std::vector<double> &W);

MatrixXd evaldRdArea_FD(
    std::vector<double> Flux,
    std::vector<double> area,
    std::vector<double> W);

#endif
