#include "structures.h"
#include "convert.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

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

