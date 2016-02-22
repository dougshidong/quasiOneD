#ifndef adjoint_h
#define adjoint_h

#include<vector>
#include <complex>

std::vector <std::complex<double> > adjoint(std::vector <std::complex<double> > x, 
                             std::vector <std::complex<double> > dx, 
                             std::vector <std::complex<double> > S,
                             std::vector <std::complex<double> > W,
                             std::vector <std::complex<double> > &psi,
                             std::vector <std::complex<double> > designVar);
#endif
