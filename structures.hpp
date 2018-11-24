#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <math.h>
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Core>

template<class T> void UNUSED( const T& ) { }

template<typename dreal> using Matrix3 = Eigen::Matrix<dreal, 3, 3>;
template<typename dreal> using Vector3 = Eigen::Matrix<dreal, 3, 1>;

template<typename dreal>
class Flow_data {

	public:
		int n_elem, n_face, n_resi, n_flux;
		std::vector<dreal> dt;
		std::vector<dreal> W;
		std::vector<dreal> W_stage1;
		std::vector<dreal> W_stage2;
		std::vector<dreal> fluxes;
		std::vector<dreal> residual;

		dreal current_CFL;
		dreal old_residual_norm;
		dreal current_residual_norm;

		Flow_data() {}
		Flow_data(const int n_elem_input) {
			const int n_elem_ghost = n_elem_input + 2;
			const int n_resi_alloc = 3*n_elem_ghost;
			n_elem = n_elem_input;
			n_resi = n_elem*3;
			n_face = n_elem+1;
			n_flux = n_face*3;

			dt.resize(n_elem_ghost);
			W.resize(n_resi_alloc);
			W_stage1.resize(n_resi_alloc);
			W_stage2.resize(n_resi_alloc);
			fluxes.resize(n_flux);
			residual.resize(n_resi_alloc);
		};
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
	double CFL_min, CFL_max, CFL_ramp;
	double flow_tol;

	int time_scheme, flux_scheme;
	double scalar_d_eps;

	int print_freq, print_conv, print_solution;

	double gam, R, Cv;
	double a2;
	double inlet_mach, inlet_total_T, inlet_total_p, outlet_p;
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
	struct Flow_options         *flow_options;
	struct Optimization_options<double> *optimization_options;

};

#endif
