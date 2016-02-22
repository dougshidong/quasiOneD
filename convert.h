#ifndef CONVERT_H
#define CONVERT_H

#include<vector>

void WtoP(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &rho,
          std::vector <std::complex<double> > &u,
          std::vector <std::complex<double> > &e);

void WtoP(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &rho,
          std::vector <std::complex<double> > &e,
          std::vector <std::complex<double> > &u,
          std::vector <std::complex<double> > &p,
          std::vector <std::complex<double> > &c,
          std::vector <std::complex<double> > &T);

void dWpdW(std::vector <std::complex<double> > &M,
           std::vector <std::complex<double> > W,
           int k);

void dWdWp(std::vector <std::complex<double> > &M,
           std::vector <std::complex<double> > W,
           int k);

void WtoF(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &F);

void WtoQ(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &Q,
          std::vector <std::complex<double> > S);
#endif
