#ifndef FITNESS_H
#define FITNESS_H
#include <vector>
#include "structures.hpp"
double evalFitness(
    const std::vector<double> &dx,
	const struct Flow_options &flow_options,
    const std::vector<double> &W,
	const struct Optimization_options<double> &opt_options);
#endif
