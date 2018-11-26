#include "structures.hpp"
#include <iostream>
#include <vector>
#include <fenv.h>
#include "input.hpp"
#include "grid.hpp"
#include "spline.hpp"
#include "quasiOneD.hpp"
#include "second_order_flow.hpp"
#include "convert.hpp"
#include "optimizer.hpp"
#include "oneshot_adjoint.hpp"
#include "oneshot_dwdx.hpp"
#include"petsc.h"
#include"petscsys.h"

static char help[] = "QuasiOneD\n\n";
int main(int argc,char **argv)
{
    feraiseexcept(FE_INVALID | FE_OVERFLOW); // Will crash on Nan or Overflow
	struct Constants  *const constants			        = new Constants; // Returned
	struct Flow_options *const flo_opts         = new Flow_options; // Returned
	struct Optimization_options<double> *const opt_opts = new Optimization_options<double>; // Returned
	opt_opts->initial_design = new Design<double>;
	opt_opts->target_design = new Design<double>;
	opt_opts->current_design = new Design<double>;

	struct Input_data *const input_data         = new Input_data; // Returned
	input_data->constants			 = constants;
	input_data->flow_options         = flo_opts;
	input_data->optimization_options = opt_opts;

	std::string input_filename = "input.in";
	inputfile(input_filename, input_data);

	const int n_elem = input_data->flow_options->n_elem;


    std::vector<double> x = uniform_x(input_data->flow_options->grid_xstart, 
				  input_data->flow_options->grid_xend, 
				  input_data->flow_options->n_elem);
    std::vector<double> dx = eval_dx(x);

    PetscInitialize(&argc, &argv, (char*)0,help);
    if (input_data->optimization_options->perform_design == 0) {
		const struct Design<double> initial_design = *(input_data->optimization_options->initial_design); // Make a copy
        const std::vector<double> area = evalS(initial_design, x, dx);

		struct Flow_options flow_options = *(input_data->flow_options); // Make a copy
		class Flow_data<double> flow_data(n_elem);
        quasiOneD(x, area, flow_options, &flow_data);
    }
    else if (abs(input_data->optimization_options->perform_design) == 1) {
		class Flow_data<double> flow_data(n_elem);

		const struct Flow_options flow_options = *(input_data->flow_options); // Make a copy

		// Target design with sine parametrization
		printf("Creating target pressure...\n");
        std::vector<double> area = evalS(*(input_data->optimization_options->target_design), x, dx);
        quasiOneD(x, area, flow_options, &flow_data);

		input_data->optimization_options->target_pressure.resize(n_elem+2);
		get_all_p(flow_options.gam, flow_data.W, input_data->optimization_options->target_pressure);

		const struct Optimization_options<double> opt_options = *(input_data->optimization_options); // Make a copy

		// Initial design with sine parametrization
		struct Design<double> *initial_design = opt_options.initial_design;
		initial_design->design_variables.resize(initial_design->n_design_variables);
        area = evalS(*initial_design, x, dx);

		// Fit a B-Spline through the sine-parametrized area
		// Note the control points contain end-points which are not part of the design variables
		int n_control_pts = initial_design->n_design_variables+2;
		std::vector<double> control_points =
			fit_bspline(x, dx, area, n_control_pts, initial_design->spline_degree);
		for (int i=0; i<initial_design->n_design_variables; i++) {
			initial_design->design_variables[i] = control_points[i+1];
		}
		// Re-evaluate the are based on the B-spline-parametrization
		initial_design->parametrization = 2;

		optimizer(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);

        //area = evalS(*initial_design, x, dx);
        //quasiOneD(x, area, flow_options, &flow_data);
    } else if (input_data->optimization_options->perform_design >= 2) {
		class Flow_data<double> flow_data(n_elem);

		const struct Flow_options flow_options = *(input_data->flow_options); // Make a copy

		// Target design with sine parametrization
		printf("Creating target pressure...\n");
        std::vector<double> area = evalS(*(input_data->optimization_options->target_design), x, dx);
        quasiOneD(x, area, flow_options, &flow_data);

		input_data->optimization_options->target_pressure.resize(n_elem+2);
		get_all_p(flow_options.gam, flow_data.W, input_data->optimization_options->target_pressure);

		const struct Optimization_options<double> opt_options = *(input_data->optimization_options); // Make a copy

		// Initial design with sine parametrization
		struct Design<double> *initial_design = opt_options.initial_design;
		initial_design->design_variables.resize(initial_design->n_design_variables);
        area = evalS(*initial_design, x, dx);

		// Fit a B-Spline through the sine-parametrized area
		// Note the control points contain end-points which are not part of the design variables
		int n_control_pts = initial_design->n_design_variables+2;
		std::vector<double> control_points =
			fit_bspline(x, dx, area, n_control_pts, initial_design->spline_degree);
		for (int i=0; i<initial_design->n_design_variables; i++) {
			initial_design->design_variables[i] = control_points[i+1];
		}
		// Re-evaluate the are based on the B-spline-parametrization
		initial_design->parametrization = 2;

		if (input_data->optimization_options->perform_design == 2) {
			oneshot_adjoint(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);
		} else if (input_data->optimization_options->perform_design >= 3) {
			oneshot_dwdx(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);
		}

        //area = evalS(*initial_design, x, dx);
        //quasiOneD(x, area, flow_options, &flow_data);
    }
    PetscFinalize();

	delete opt_opts->initial_design;
	delete opt_opts->target_design;
	delete opt_opts->current_design;
	delete constants;
	delete flo_opts;
	delete opt_opts;
	delete input_data;

    return 0;
}
