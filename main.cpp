#include "structures.h"
#include <iostream>
#include <vector>
#include <fenv.h>
#include "input.h"
#include "grid.h"
#include "spline.h"
#include "quasiOneD.h"
#include "second_order_flow.h"
#include "convert.h"
#include "optimizer.h"
#include "oneshot.h"
#include"petsc.h"
#include"petscsys.h"

static char help[] = "QuasiOneD\n\n";
int main(int argc,char **argv)
{
    feraiseexcept(FE_INVALID | FE_OVERFLOW); // Will crash on Nan or Overflow
	struct Constants  *const constants			= new Constants; // Returned
	struct Flow_options *const flo_opts         = new Flow_options; // Returned
	struct Optimization_options *const opt_opts = new Optimization_options; // Returned
	opt_opts->initial_design = new Design;
	opt_opts->target_design = new Design;
	opt_opts->current_design = new Design;

	struct Input_data *const input_data         = new Input_data; // Returned
	input_data->constants			 = constants;
	input_data->flow_options         = flo_opts;
	input_data->optimization_options = opt_opts;

	std::string input_filename = "input.in";
	inputfile(input_filename, input_data);

	int n_elem = input_data->flow_options->n_elem;

    std::vector<double> x(n_elem), area(n_elem + 1);
    std::vector<double> dx(n_elem);

    x = uniform_x(input_data->flow_options->grid_xstart, 
				  input_data->flow_options->grid_xend, 
				  input_data->flow_options->n_elem);
    dx = eval_dx(x);

    if (input_data->optimization_options->perform_design == -1) {
		const struct Design initial_design = *(input_data->optimization_options->initial_design); // Make a copy
        area = evalS(initial_design, x, dx);
        second_order_flow(*constants, x, area, *(input_data->flow_options));
    }
    if (input_data->optimization_options->perform_design == 0) {
		const struct Design initial_design = *(input_data->optimization_options->initial_design); // Make a copy
        area = evalS(initial_design, x, dx);

		struct Flow_options flow_options = *(input_data->flow_options); // Make a copy
		struct Flow_data flow_data;
		flow_data.dt.resize(n_elem);
		flow_data.W.resize(3*n_elem);
		flow_data.W_stage.resize(3*n_elem);
		flow_data.fluxes.resize(3*n_elem);
		flow_data.residual.resize(3*n_elem);
        quasiOneD(x, area, flow_options, &flow_data);
    }
    else if (input_data->optimization_options->perform_design == 1) {
		struct Flow_data flow_data;
		flow_data.dt.resize(n_elem);
		flow_data.W.resize(3*n_elem);
		flow_data.W_stage.resize(3*n_elem);
		flow_data.fluxes.resize(3*n_elem);
		flow_data.residual.resize(3*n_elem);


		const struct Flow_options flow_options = *(input_data->flow_options); // Make a copy

		// Target design with sine parametrization
		printf("Creating target pressure...\n");
        area = evalS(*(input_data->optimization_options->target_design), x, dx);
        quasiOneD(x, area, flow_options, &flow_data);

		input_data->optimization_options->target_pressure.resize(n_elem);
		get_all_p(flow_options.gam, flow_data.W, input_data->optimization_options->target_pressure);

		const struct Optimization_options opt_options = *(input_data->optimization_options); // Make a copy

		// Initial design with sine parametrization
		struct Design *initial_design = opt_options.initial_design;
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

        PetscInitialize(&argc, &argv, (char*)0,help);
		optimizer(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);
        PetscFinalize();

        //area = evalS(*initial_design, x, dx);
        //quasiOneD(x, area, flow_options, &flow_data);
    } else if (input_data->optimization_options->perform_design == 2) {
		struct Flow_data flow_data;
		flow_data.dt.resize(n_elem);
		flow_data.W.resize(3*n_elem);
		flow_data.W_stage.resize(3*n_elem);
		flow_data.fluxes.resize(3*n_elem);
		flow_data.residual.resize(3*n_elem);

		const struct Flow_options flow_options = *(input_data->flow_options); // Make a copy

		// Target design with sine parametrization
		printf("Creating target pressure...\n");
        area = evalS(*(input_data->optimization_options->target_design), x, dx);
        quasiOneD(x, area, flow_options, &flow_data);

		input_data->optimization_options->target_pressure.resize(n_elem);
		get_all_p(flow_options.gam, flow_data.W, input_data->optimization_options->target_pressure);

		const struct Optimization_options opt_options = *(input_data->optimization_options); // Make a copy

		// Initial design with sine parametrization
		struct Design *initial_design = opt_options.initial_design;
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

		oneshot_dwdx(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);
		//oneshot(*constants, x, dx, *flo_opts, *opt_opts, *initial_design);

        //area = evalS(*initial_design, x, dx);
        //quasiOneD(x, area, flow_options, &flow_data);
    }

	delete opt_opts->initial_design;
	delete opt_opts->target_design;
	delete opt_opts->current_design;
	delete constants;
	delete flo_opts;
	delete opt_opts;
	delete input_data;

    return 0;
}
