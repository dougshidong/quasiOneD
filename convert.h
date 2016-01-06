#ifndef CONVERT_H
#define CONVERT_H

#include<vector>

void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &u,
          std::vector <double> &e);

void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &u,
          std::vector <double> &p,
          std::vector <double> &c);

void WtoF(std::vector <double> W,
          std::vector <double> &F);

#endif
