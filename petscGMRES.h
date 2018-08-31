#ifndef petscgmres_h
#define petscgmres_h
#include<Eigen/Core>
#include<Eigen/Sparse>
using namespace Eigen;
MatrixXd solveGMRES(SparseMatrix <double> A, MatrixXd B);
#endif
