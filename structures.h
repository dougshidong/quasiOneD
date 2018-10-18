#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <math.h>
#include <string>
#include <vector>

template<class T> void UNUSED( const T& ) { }

template<typename dreal>
struct Flow_data {
	std::vector<dreal> dt;
	std::vector<dreal> W;
	std::vector<dreal> W_stage;
	std::vector<dreal> fluxes;
	std::vector<dreal> residual;
};

struct Constants {
	const double PI = atan(1.0) * 4.0;
	std::string case_name;
};
template<typename dreal>
struct Flow_options {
	std::string case_name;
	int n_elem;
	dreal grid_xstart, grid_xend;

	int flow_maxit;
	dreal CFL;
	dreal flow_tol;

	int time_scheme, flux_scheme;
	dreal scalar_d_eps;

	int print_freq, print_conv, print_solution;

	dreal gam, R, Cv;
	dreal a2;
	dreal inlet_mach, inlet_total_T, inlet_total_p, outlet_p;
};

template<typename dreal>
struct Design {
	dreal h, t1, t2;

	int spline_degree;

	int parametrization;
	int n_design_variables;

	std::vector<dreal> design_variables;
};

template<typename dreal>
struct Optimization_options {
	int perform_design;
	int opt_maxit;
	int cost_function;
	std::vector<dreal> target_pressure;
	int descent_type, gradient_type, hessian_type, exact_hessian;
	dreal opt_tol, adj_tol, dwdx_tol;
    dreal default_pert=1e-4;

	struct Design<dreal> *initial_design;
	struct Design<dreal> *current_design;
	struct Design<dreal> *target_design;

	int n_design_variables, spline_degree;
};

struct Input_data{
	struct Constants                   *constants;
	struct Flow_options<double>         *flow_options;
	struct Optimization_options<double> *optimization_options;

};

#endif
