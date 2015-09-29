#ifndef adjoint_h
#define adjoint_h

#include<vector>

void adjointBC(int nx,
        std::vector <double> &psi, 
        std::vector <double> W,
        std::vector <double> dx,
        std::vector <double> S);

std::vector <double> adjointSteger(int nx,
                                   std::vector <double> S,
                                   std::vector <double> dx,
                                   std::vector <double> W,
                                   std::vector <double> psi,
                                   std::vector <double> &pFlux);
#endif
