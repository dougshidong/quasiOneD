#ifndef residuald2_h
#define residuald2_h
#include<Eigen/Sparse>
#include<vector>
std::vector < Eigen::SparseMatrix <double> > evalddRdWdW(
    std::vector <double> W,
    std::vector <double> S);

std::vector < Eigen::SparseMatrix <double> > evalddRdWdW_FD(
    std::vector <double> W,
    std::vector <double> S);

std::vector < Eigen::SparseMatrix <double> > evalddRdWdS(
    std::vector <double> W,
    std::vector <double> S);

std::vector < Eigen::SparseMatrix <double> > evalddRdWdS_FD(
    std::vector <double> W,
    std::vector <double> S);

std::vector < Eigen::Matrix3d> ddWpdWdWp(
    std::vector <double> W,
    int i);
#endif
