#ifndef objectiveDerivatives_h
#define obejctiveDerivatives_h
#include<Eigen/Core>
#include<vector>

Eigen::VectorXd evaldIcdW(
    std::vector <double> W,
    std::vector <double> dx);

SparseMatrix <double> evaldIcdWdW(
    std::vector <double> W,
    std::vector <double> dx);

Eigen::VectorXd evaldIcdS();
Eigen::MatrixXd evalddIcdSdS();
Eigen::MatrixXd evalddIcdWdS();
#endif
