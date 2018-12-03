#include "quasiOneD.hpp"
#include "structures.hpp"
#include "convert.hpp"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <vector>
#include "grid.hpp"
#include "flux.hpp"
#include "timestep.hpp"
#include "output.hpp"

template<typename dreal>
int quasiOneD(
	const std::vector<dreal> &x,
	const std::vector<dreal> &area,
	const Flow_options &flo_opts,
	class Flow_data<dreal>* const flow_data)
{
	const dreal gam = flo_opts.gam;
	const int n_elem = flo_opts.n_elem;

	if(n_elem+1 != (int)x.size()) abort();

    std::vector<dreal> dx = eval_dx(x);

    std::vector<int> itV;
    std::vector<dreal> normV;
    std::vector<dreal> timeVec;


    // Inlet flow properties
    dreal inlet_T = isenT(gam, flo_opts.inlet_total_T, flo_opts.inlet_mach);
    const dreal p_inlet = isenP(gam, flo_opts.inlet_total_p, flo_opts.inlet_mach);
    dreal p = p_inlet;
    dreal rho = p / (flo_opts.R * inlet_T);
    dreal c = sqrt(gam * p / rho);
    dreal u = flo_opts.inlet_mach * c;
    dreal e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

	flow_data->W[0*3+0] = rho;
	flow_data->W[0*3+1] = rho * u;
	flow_data->W[0*3+2] = e;
    // Flow properties initialization with outlet
    // State Vectors Initialization
    for (int i = 1; i < n_elem+2; i++) {
        p = p_inlet + (i/(n_elem+1)) * (flo_opts.outlet_p - p_inlet);
		const dreal T = flo_opts.inlet_total_T * pow(p/flo_opts.inlet_total_p,(gam-1.0)/gam);
        //p = flo_opts.outlet_p;
        rho = p / (flo_opts.R * T);
        c = sqrt(gam * p / rho);
        u = c * flo_opts.inlet_mach;
        e = rho * (flo_opts.Cv * T + 0.5 * pow(u, 2));

        flow_data->W[i*3+0] = rho;
        flow_data->W[i*3+1] = rho * u;
        flow_data->W[i*3+2] = e;
    }


    //clock_t tic = clock();
    //clock_t toc;
    //dreal elapsed;

    getDomainResi(flo_opts, area, flow_data);
    flow_data->current_residual_norm = 1.0;
    int iterations = 0;
    while(flow_data->current_residual_norm > flo_opts.flow_tol && iterations < flo_opts.flow_maxit) {
        iterations++;

        // Step in Time
        stepInTime(flo_opts, area, dx, flow_data);

        // Monitor Convergence
        if (iterations%flo_opts.print_freq == 0) {
            if (flo_opts.print_conv == 1) {
                std::cout<<"Iteration "<<iterations<<"    CFL "<<flow_data->current_CFL<<"   NormR "<<std::setprecision(15)<<flow_data->current_residual_norm<<std::endl;
            }
            //itV[iterations / flo_opts.print_freq - 1] = iterations;
            //normV[iterations / flo_opts.print_freq - 1] = residual_norm;

			double max_res = 0;
			int max_loc = 0;
			for (int i = 1; i < n_elem+1; i++) {
				for (int istate = 1; istate < 3; istate++) {
					if (abs(flow_data->residual[i*3+istate]) > abs(max_res)) {
						max_res = flow_data->residual[i*3+istate];
						max_loc = i;
					}
				}
			}
            //if (flo_opts.print_conv == 1) {
			//	std::cout<<"Max loc and residual "<<max_loc<<" "<<max_res<<std::endl;
            //}

            //toc = clock();
            //elapsed = (dreal)(toc-tic) / CLOCKS_PER_SEC;
            //timeVec[iterations / flo_opts.print_freq - 1] = elapsed;

			if (flo_opts.print_solution == 1) {
				for (int k = 0; k < 3; k++) {
					std::cout<<"W"<<k + 1<<std::endl;
					for (int i = 0; i < n_elem; i++) {
						printf("%d %25.15e\n", i, flow_data->W[i*3+k]);
					}
				}
			}
        }
    }

	std::vector<dreal> p_vec(n_elem+2);
	for (int i = 0; i < n_elem+2; i++) {
		p_vec[i] = get_p(flo_opts.gam, flow_data->W[i*3+0], flow_data->W[i*3+1], flow_data->W[i*3+2]);
	}
    outVec(flo_opts.case_name, "current_pressure.dat", "w", p_vec);
    std::cout<<"Flow iterations = "<<iterations<<"   Density Residual = "<<flow_data->current_residual_norm<<std::endl;
    return 0;
}
template int quasiOneD( const std::vector<double> &x, const std::vector<double> &area, const Flow_options &flo_opts, class Flow_data<double>* const flow_data);

