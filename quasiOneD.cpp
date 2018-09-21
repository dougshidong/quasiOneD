#include "structures.h"
#include "convert.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <vector>
#include "flux.h"
#include "timestep.h"
#include "output.h"

double isenP(const double gam, const double pt, const double M) {
	return pt * pow((1 + (gam - 1) / 2 * pow(M, 2)), ( - gam / (gam - 1)));
}

double isenT(const double gam,double Tt, const double M) {
	return Tt * pow((1 + (gam - 1) / 2 * pow(M, 2)), - 1);
}

void inletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);
void outletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);

double quasiOneD(
	const std::vector<double> &x,
	const std::vector<double> &area,
	const Flow_options &flo_opts,
	struct Flow_data* const flow_data)
{
	const double gam = flo_opts.gam;
	const int n_elem = flo_opts.n_elem;

	if(n_elem!=x.size()) abort();

    std::vector<double> dx(n_elem);
    for (int i = 1; i < n_elem-1; i++) {
        dx[i] =  (x[i] - x[i-1])/2  +  (x[i+1] - x[i])/2 ;
    }
    dx[0] = x[1] - x[0];
    dx[n_elem-1] = x[n_elem-1] - x[n_elem-2];

    std::vector<double> dt(n_elem);

    std::vector <int> itV(flo_opts.flow_maxit/flo_opts.print_freq);
    std::vector<double> normV(flo_opts.flow_maxit/flo_opts.print_freq);
    std::vector<double> timeVec(flo_opts.flow_maxit/flo_opts.print_freq);


    // Inlet flow properties
    double inlet_T = isenT(gam, flo_opts.inlet_total_T, flo_opts.inlet_mach);
    double p = isenP(gam, flo_opts.inlet_total_p, flo_opts.inlet_mach);
    double rho = p / (flo_opts.R * inlet_T);
    double c = sqrt(gam * p / rho);
    double u = flo_opts.inlet_mach * c;
    double e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

	flow_data->W[0*3+0] = rho;
	flow_data->W[0*3+1] = rho * u;
	flow_data->W[0*3+2] = e;
    // Flow properties initialization with outlet
    // State Vectors Initialization
    for (int i = 1; i < n_elem; i++) {
        p = flo_opts.outlet_p;
        rho = p / (flo_opts.R * inlet_T);
        c = sqrt(gam * p / rho);
        u = c * flo_opts.inlet_mach;
        e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

        flow_data->W[i*3+0] = rho;
        flow_data->W[i*3+1] = rho * u;
        flow_data->W[i*3+2] = e;
    }


    clock_t tic = clock();
    clock_t toc;
    double elapsed;

    double normR = 1.0;
    int iterations = 0;
    while(normR > flo_opts.flow_tol && iterations < flo_opts.flow_maxit) {
        iterations++;

        // Calculate Time Step
        for (int i = 0; i < n_elem; i++) {
			double u = flow_data->W[i*3+1] / flow_data->W[i*3+0];
			double c = get_c(gam, flow_data->W[i*3+0], flow_data->W[i*3+1], flow_data->W[i*3+2]);
            dt[i] = (flo_opts.CFL * dx[i]) / fabs(u + c);
		}

        // Step in Time
        stepInTime(flo_opts, area, dx, dt, flow_data);

		const int first_cell = 0;
		const int last_cell = n_elem-1;
        inletBC(flo_opts, dt[first_cell], dx[first_cell], flow_data);

        outletBC(flo_opts, dt[last_cell], dx[last_cell], flow_data);

        // Calculating the norm of the density residual
        normR = 0;
        for (int i = 0; i < n_elem; i++)
            normR = normR + pow(flow_data->residual[i*3+0], 2);
        normR = sqrt(normR);

        // Monitor Convergence
        if (iterations%flo_opts.print_freq == 0) {
            if (flo_opts.print_conv == 1) {
                std::cout<<"Iteration "<<iterations <<"   NormR "<<std::setprecision(15)<<normR<<std::endl;
            }
            itV[iterations / flo_opts.print_freq - 1] = iterations;
            normV[iterations / flo_opts.print_freq - 1] = normR;

            toc = clock();
            elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
            timeVec[iterations / flo_opts.print_freq - 1] = elapsed;

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

	std::vector<double> p_vec(n_elem);
	for (int i = 0; i < n_elem; i++) {
		p_vec[i] = get_p(flo_opts.gam, flow_data->W[i*3+0], flow_data->W[i*3+1], flow_data->W[i*3+2]);
	}
    outVec(flo_opts.case_name, "current_pressure.dat", "w", p_vec);
    std::cout<<"Flow iterations = "<<iterations<<"   Density Residual = "<<normR<<std::endl;
    return 1;
}


void inletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data)
{
	const double a2  = flo_opts.a2;
	const double gam = flo_opts.gam;
    double rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[i*3+0];
        u[i] = flow_data->W[i*3+1] / rho[i];
        e[i] = flow_data->W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
    if (u[0] < c[0]) {
        const double dpdu = flo_opts.inlet_total_p * (gam / (gam - 1.0))
             * pow(1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2,
                   1.0 / (gam - 1.0))
             * ( - 2.0 * ((gam - 1.0) / (gam + 1.0)) * u[0] / a2);
        const double dtdx = dt0 / dx0;
        const double eigenvalue = ((u[1] + u[0] - c[1] - c[0]) / 2.0) * dtdx;

        const double dpdx = p[1] - p[0];
        const double dudx = u[1] - u[0];
        const double du = -eigenvalue * (dpdx - rho[0] * c[0] * dudx)
                / (dpdu - rho[0] * c[0]);

//      flow_data->residual[0 * 3 + 1] = -((u[0] + du) - u[0]) / dtdx;
        u[0] = u[0] + du;

        const double T0 = flo_opts.inlet_total_T * (1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2);
//      flow_data->residual[0 * 3 + 2] = -(flo_opts.inlet_total_p * pow(T0 / flo_opts.inlet_total_T, gam / (gam - 1.0)) - p[0]) / dtdx;
        p[0] = flo_opts.inlet_total_p * pow(T0 / flo_opts.inlet_total_T, gam / (gam - 1.0));
        rho[0] = p[0] / (flo_opts.R * T0);
        e[0] = rho[0] * (flo_opts.Cv * T0 + 0.5 * u[0] * u[0]);

        flow_data->residual[0 * 3 + 0] = -(rho[0] - flow_data->W[0 * 3 + 0]) / dtdx;
        flow_data->residual[0 * 3 + 1] = -(rho[0] * u[0] - flow_data->W[0 * 3 + 1]) / dtdx;
        flow_data->residual[0 * 3 + 2] = -(e[0] - flow_data->W[0 * 3 + 2]) / dtdx;

        flow_data->W[0 * 3 + 0] = rho[0];
        flow_data->W[0 * 3 + 1] = rho[0] * u[0];
        flow_data->W[0 * 3 + 2] = e[0];
    } else {
        flow_data->residual[0 * 3 + 0] = 0;
        flow_data->residual[0 * 3 + 1] = 0;
        flow_data->residual[0 * 3 + 2] = 0;
    }
}

void outletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data)
{
	const int n_elem = flo_opts.n_elem;
	const double gam = flo_opts.gam;
    double avgc, avgu, dtdx, MachOut;
    double eigenvalues[3], Ri[3];
    double dpdx, dudx, du, drho, dp, T;
    double rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[(i + (n_elem - 2)) * 3 + 0];
        u[i] = flow_data->W[(i + (n_elem - 2)) * 3 + 1] / rho[i];
        e[i] = flow_data->W[(i + (n_elem - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    avgu = (u[1] + u[0]) / 2;
    avgc = (c[1] + c[0]) / 2;
    dtdx = dt0 / dx0;
    eigenvalues[0] = avgu * dtdx;
    eigenvalues[1] = (avgu + avgc) * dtdx;
    eigenvalues[2] = (avgu - avgc) * dtdx;

    dpdx = p[1] - p[0];
    dudx = u[1] - u[0];

    Ri[0] = -eigenvalues[0] * ( rho[1] - rho[0] - dpdx / pow(c[1], 2) );
    Ri[1] = -eigenvalues[1] * ( dpdx + rho[1] * c[1] * dudx );
    Ri[2] = -eigenvalues[2] * ( dpdx - rho[1] * c[1] * dudx );

    MachOut = avgu / avgc;
    if (MachOut > 1) {
        dp = 0.5 * (Ri[1] + Ri[2]);
    } else {
        dp = 0;
    }

    drho = Ri[0] + dp / (pow(c[1], 2));
    du = (Ri[1] - dp) / (rho[1] * c[1]);

    u[1] = u[1] + du;
    rho[1] = rho[1] + drho;
    p[1] = p[1] + dp;
    T = p[1] / (rho[1] * flo_opts.R);
    e[1] = rho[1] * (flo_opts.Cv * T + 0.5 * pow(u[1], 2));

    flow_data->residual[(n_elem - 1) * 3 + 0] = (flow_data->W[(n_elem - 1) * 3 + 0] - rho[1]) / dtdx;
    flow_data->residual[(n_elem - 1) * 3 + 1] = (flow_data->W[(n_elem - 1) * 3 + 1] - rho[1] * u[1]) / dtdx;
    flow_data->residual[(n_elem - 1) * 3 + 2] = (flow_data->W[(n_elem - 1) * 3 + 2] - e[1]) / dtdx;

    flow_data->W[(n_elem - 1) * 3 + 0] = rho[1];
    flow_data->W[(n_elem - 1) * 3 + 1] = rho[1] * u[1];
    flow_data->W[(n_elem - 1) * 3 + 2] = e[1];
}
