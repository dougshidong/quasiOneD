#include "input.hpp"
#include "structures.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <string>
#include <algorithm> // remove_if

#define MAX_STRLEN 256

void inputfile(std::string filename, struct Input_data* const input_data)
{
	struct Constants* const constants            = input_data->constants;
	struct Flow_options* const flo_opts          = input_data->flow_options;
	struct Optimization_options<double>* const opt_opts  = input_data->optimization_options;

    FILE *file = fopen("input.in", "r");

    char buf[MAX_STRLEN];
	// Name of test case. Will affect output names
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
	std::string case_name = buf;
    case_name.erase(std::remove_if (case_name.begin(), case_name.end(), isspace), case_name.end());
	constants->case_name = case_name;
	flo_opts->case_name = case_name;


    // Grid size
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%d %lf %lf\n", &flo_opts->n_elem, &flo_opts->grid_xstart, &flo_opts->grid_xend);

    // Input Geometry
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf %lf", &opt_opts->initial_design->h, &opt_opts->initial_design->t1, &opt_opts->initial_design->t2);

    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    // Flow Solver Paramameter
    // Stepping Scheme
    // 0   -   Euler Explicit
    // 1   -   Runge - Kutta 4th order
    // 2   -   Jameson's Runge-Kutta 4th order
    // Flux Scheme
    // 0   -   Scalar Dissipation (SD)
    // 1   -   Steger Warming (SW)
    sscanf(buf, "%d %d %lf", &flo_opts->time_scheme, &flo_opts->flux_scheme, &flo_opts->scalar_d_eps);

    // Flow Convergence
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf %d", &flo_opts->CFL, &flo_opts->flow_tol, &flo_opts->flow_maxit);

    // Printing Flow Stuff for Debugging
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%d %d %d", &flo_opts->print_freq, &flo_opts->print_conv, &flo_opts->print_solution);

    // Flow Inputs
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf", &flo_opts->gam, &flo_opts->R);
    flo_opts->Cv = flo_opts->R / (flo_opts->gam - 1.0);

    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    double outlet_ratio;
    sscanf(buf, "%lf %lf %lf %lf", &flo_opts->inlet_mach, &flo_opts->inlet_total_T, &flo_opts->inlet_total_p, &outlet_ratio);
    flo_opts->outlet_p = outlet_ratio * flo_opts->inlet_total_p;
    flo_opts->a2 = 2.0 * flo_opts->gam * flo_opts->Cv * flo_opts->inlet_total_T * ((flo_opts->gam - 1.0) / (flo_opts->gam + 1.0));

    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    // Design Optimization Parameters
    // perform_design = 0 or 1 to Turn ON/OFF Optimization
    // Design Variables
    // 0  -  Individual Areas
    // 1  -  Sin Parametrization (Final Project MECH 539)
    // Fitness Function
    // 0  -  Total Pressure Loss
    // 1  -  Pressure Target
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
	int design_parametrization;
    sscanf(buf, "%d %d %d", &opt_opts->perform_design, &design_parametrization, &opt_opts->cost_function);

    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    // Descent Type for Optimization
    // 1  -  Steepest Descent
    // 2  -  Quasi-Newton (BFGS)
    // 3  -  Newton
    // 4  -  Truncated-Newton
    // Gradient Type
    //-3  -  FD Centered
    //-2  -  FD Centered
    //-1  -  FD Centered
    // 1  -  Adjoint Variable
    // 2  -  Direct Differentiation
    // Hessian Type
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%d %d %d %d %lf",
		&opt_opts->descent_type,
		&opt_opts->gradient_type,
		&opt_opts->hessian_type,
		&opt_opts->exact_hessian,
		&opt_opts->dwdx_tol);

    // Design Convergence
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%lf %d", &opt_opts->opt_tol, &opt_opts->opt_maxit);

    // Target Geometry
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf %lf", &opt_opts->target_design->h, &opt_opts->target_design->t1, &opt_opts->target_design->t2);

    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), file) == NULL) {abort();} // Read
	int n_design_variables, spline_degree;
    sscanf(buf, "%d %d", &n_design_variables, &spline_degree);

    opt_opts->n_design_variables                 = n_design_variables;
    opt_opts->spline_degree                      = spline_degree;

    opt_opts->initial_design->n_design_variables = n_design_variables;
    opt_opts->target_design->n_design_variables  = n_design_variables;
    opt_opts->current_design->n_design_variables = n_design_variables;

    opt_opts->initial_design->parametrization    = 1;//design_parametrization;
    opt_opts->initial_design->spline_degree      = spline_degree;

    opt_opts->target_design->parametrization     = 1;//opt_opts->initial_design->parametrization;
    opt_opts->target_design->spline_degree       = opt_opts->initial_design->spline_degree;

    opt_opts->current_design->parametrization    = opt_opts->initial_design->parametrization;
    opt_opts->current_design->spline_degree      = opt_opts->initial_design->spline_degree;

    fclose(file);
}

