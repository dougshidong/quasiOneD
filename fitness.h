#ifndef fitness_h
#define fitness_h
#include <vector>
#include "structures.h"
double evalFitness(
    const std::vector<double> &dx,
	const struct Flow_options &flow_options,
    const std::vector<double> &W,
	const struct Optimization_options &opt_options);
#endif
