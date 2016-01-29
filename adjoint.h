#ifndef adjoint_h
#define adjoint_h

#include<vector>

std::vector <double> adjoint(std::vector <double> x, 
                             std::vector <double> dx, 
                             std::vector <double> S,
                             std::vector <double> W,
                             std::vector <double> &psi,
                             std::vector <double> designVar);
#endif
