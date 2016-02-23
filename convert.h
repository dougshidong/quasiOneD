#ifndef CONVERT_H
#define CONVERT_H

#include<vector>

void WtoP(std::vector <long double> W,
          std::vector <long double> &rho,
          std::vector <long double> &u,
          std::vector <long double> &e);

void WtoP(std::vector <long double> W,
          std::vector <long double> &rho,
          std::vector <long double> &e,
          std::vector <long double> &u,
          std::vector <long double> &p,
          std::vector <long double> &c,
          std::vector <long double> &T);

void dWpdW(std::vector <long double> &M,
           std::vector <long double> W,
           int k);

void dWdWp(std::vector <long double> &M,
           std::vector <long double> W,
           int k);

void WtoF(std::vector <long double> W,
          std::vector <long double> &F);

void WtoQ(std::vector <long double> W,
          std::vector <long double> &Q,
          std::vector <long double> S);
#endif
