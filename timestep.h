#ifndef TIMESTEP_H
#define TIMESTEP_H

#include<vector>
#include<complex>

void EulerExplicitStep(std::vector <std::complex<double> > S,
                std::vector <std::complex<double> > V,
                std::vector <std::complex<double> > dt,
                std::vector <std::complex<double> > Q,
                std::vector <std::complex<double> > &Resi,
                std::vector <std::complex<double> > &W,
                std::vector <std::complex<double> > F);

void rk4(std::vector <std::complex<double> > dx, std::vector <std::complex<double> > S, 
                std::vector <std::complex<double> > dt, 
                std::vector <std::complex<double> > &W,
                std::vector <std::complex<double> > F,
                std::vector <std::complex<double> > Q, 
                std::vector <std::complex<double> > &Resi);

void jamesonrk(std::vector <std::complex<double> > dx,
         std::vector <std::complex<double> > S, 
         std::vector <std::complex<double> > V, 
         std::vector <std::complex<double> > dt, 
         std::vector <std::complex<double> > &W,
         std::vector <std::complex<double> > F,
         std::vector <std::complex<double> > &Resi);
#endif
