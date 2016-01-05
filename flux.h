#ifndef FLUX_H
#define FLUX_H

#include<vector>

void Flux_StegerWarmingV(int nx,
	       		std::vector <double> &Flux,
		       	std::vector <double> W,
		       	std::vector <double> u,
		       	std::vector <double> c,
		       	std::vector <double> rho);

void Flux_Scalar(int nx,
	       		std::vector <double> &Flux,
		       	std::vector <double> W,
                std::vector <double> F,
		       	std::vector <double> u,
		       	std::vector <double> c);
#endif
