#ifndef TIMESTEP_H
#define TIMESTEP_H

#include<vector>

void EulerExplicitStep(std::vector <long double> S,
                std::vector <long double> V,
                std::vector <long double> dt,
                std::vector <long double> Q,
                std::vector <long double> &Resi,
                std::vector <long double> &W,
                std::vector <long double> F);

void rk4(std::vector <long double> dx, std::vector <long double> S, 
                std::vector <long double> dt, 
                std::vector <long double> &W,
                std::vector <long double> F,
                std::vector <long double> Q, 
                std::vector <long double> &Resi);

void jamesonrk(std::vector <long double> dx,
         std::vector <long double> S, 
         std::vector <long double> V, 
         std::vector <long double> dt, 
         std::vector <long double> &W,
         std::vector <long double> F,
         std::vector <long double> &Resi);
#endif
