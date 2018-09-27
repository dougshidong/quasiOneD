#ifndef residuald2_h
#define residuald2_h
#include<Eigen/SparseCore>
#include<vector>
using namespace Eigen;
std::vector < SparseMatrix<double> > evalddRdWdW(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);

std::vector <SparseMatrix<double> > evalddRdWdW_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);

std::vector < Eigen::SparseMatrix<double> > evalddRdWdArea(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);

std::vector < Eigen::SparseMatrix<double> > evalddRdWdArea_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data);

#endif
