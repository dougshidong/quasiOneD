#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <math.h>
#include <string>
#include <vector>

struct Flow_data {
	std::vector<double> W;
	std::vector<double> W_stage;
	std::vector<double> fluxes;
	std::vector<double> residual;
};

struct Constants {
	const double PI = atan(1.0) * 4.0;
	std::string case_name;
};
struct Flow_options {
	std::string case_name;
	int n_elem;
	double grid_xstart, grid_xend;

	int flow_maxit;
	double CFL;
	double flow_tol;

	int time_scheme, flux_scheme;
	double scalar_d_eps;

	int print_freq, print_conv, print_solution;

	double gam, R, Cv;
	double a2;
	double inlet_mach, inlet_total_T, inlet_total_p, outlet_p;
};

struct Design {
	double h, t1, t2;

	int spline_degree;

	int parametrization;
	int n_design_variables;

	std::vector<double> design_variables;

};

struct Optimization_options {
	int perform_design;
	int opt_maxit;
	int cost_function;
	std::vector<double> target_pressure;
	int descent_type, gradient_type, hessian_type, exact_hessian;
	double opt_tol, adj_tol, dwdx_tol;

	struct Design               *initial_design;
	struct Design               *current_design;
	struct Design               *target_design;

	int n_design_variables, spline_degree;
};

struct Input_data{
	struct Constants            *constants;
	struct Flow_options         *flow_options;
	struct Optimization_options *optimization_options;
};

#endif
