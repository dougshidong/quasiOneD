#ifndef CONVERT_H
#define CONVERT_H

#include<vector>

void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &u,
          std::vector <double> &e);

void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &e,
          std::vector <double> &u,
          std::vector <double> &p,
          std::vector <double> &c,
          std::vector <double> &T);

void dWpdW(std::vector <double> &M,
           std::vector <double> W,
           int k);

void WtoF(std::vector <double> W,
          std::vector <double> &F);

void WtoQ(std::vector <double> W,
          std::vector <double> &Q,
          std::vector <double> S);
#endif
