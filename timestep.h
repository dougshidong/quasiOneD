#ifndef TIMESTEP_H
#define TIMESTEP_H

#include<vector>

void EulerExplicitStep(std::vector <double> S,
                std::vector <double> V,
                std::vector <double> dt,
                std::vector <double> Flux,
                std::vector <double> Q,
                std::vector <double> &Resi,
                std::vector <double> &W);

void rk4(std::vector <double> dx, std::vector <double> S, 
                std::vector <double> dt, 
                std::vector <double> &W,
                std::vector <double> Q, 
                std::vector <double> &Resi,
                std::vector <double> &Flux);
#endif
