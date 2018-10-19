#ifndef RESIDUALD2_H
#define RESIDUALD2_H
#include<Eigen/SparseCore>
#include<vector>
using namespace Eigen;
std::vector < SparseMatrix<double> > evalddRdWdW(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data);

std::vector <SparseMatrix<double> > evalddRdWdW_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data);

std::vector < Eigen::SparseMatrix<double> > evalddRdWdArea(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data);

std::vector < Eigen::SparseMatrix<double> > evalddRdWdArea_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data);

#endif
