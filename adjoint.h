#ifndef adjoint_h
#define adjoint_h

#include<vector>

std::vector <long double> adjoint(std::vector <long double> x, 
                             std::vector <long double> dx, 
                             std::vector <long double> S,
                             std::vector <long double> W,
                             std::vector <long double> &psi,
                             std::vector <long double> designVar);
#endif
