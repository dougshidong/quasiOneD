#include"boundary_conditions.hpp"
#include "structures.hpp"
#include "convert.hpp"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include<adolc/adolc.h>

template<typename dreal>
void inletBC(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data)
{
	const dreal a2  = flo_opts.a2;
	const dreal gam = flo_opts.gam;
    dreal rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[i*3+0];
        u[i] = flow_data->W[i*3+1] / rho[i];
        e[i] = flow_data->W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
    if (u[0] < c[0]) {
        const dreal dpdu = flo_opts.inlet_total_p * (gam / (gam - 1.0))
             * pow(1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2,
                   1.0 / (gam - 1.0))
             * ( - 2.0 * ((gam - 1.0) / (gam + 1.0)) * u[0] / a2);
        const dreal dtdx = dt0 / dx0;
        const dreal eigenvalue = ((u[1] + u[0] - c[1] - c[0]) / 2.0) * dtdx;

        const dreal dpdx = p[1] - p[0];
        const dreal dudx = u[1] - u[0];
        const dreal du = -eigenvalue * (dpdx - rho[0] * c[0] * dudx)
                / (dpdu - rho[0] * c[0]);

//      flow_data->residual[0 * 3 + 1] = -((u[0] + du) - u[0]) / dtdx;
        u[0] = u[0] + du;

        const dreal T0 = flo_opts.inlet_total_T * (1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2);
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
template void inletBC(const Flow_options &flo_opts, const double dt0, const double dx0, class Flow_data<double>* const flow_data);
template void inletBC(const Flow_options &flo_opts, const adouble dt0, const adouble dx0, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void outletBC(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data)
{
	//const int n_elem = flo_opts.n_elem;
	//const int n_elem_ghost = flo_opts.n_elem+2;
	const int n_elem_ghost = flow_data->W.size()/3;
	const dreal gam = flo_opts.gam;
    dreal eigenvalues[3], Ri[3];
    dreal rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 0];
        u[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 1] / rho[i];
        e[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    const dreal avgu = (u[1] + u[0]) / 2;
    const dreal avgc = (c[1] + c[0]) / 2;
    const dreal dtdx = dt0 / dx0;
    eigenvalues[0] = avgu * dtdx;
    eigenvalues[1] = (avgu + avgc) * dtdx;
    eigenvalues[2] = (avgu - avgc) * dtdx;

    const dreal dpdx = p[1] - p[0];
    const dreal dudx = u[1] - u[0];

    Ri[0] = -eigenvalues[0] * ( rho[1] - rho[0] - dpdx / pow(c[1], 2) );
    Ri[1] = -eigenvalues[1] * ( dpdx + rho[1] * c[1] * dudx );
    Ri[2] = -eigenvalues[2] * ( dpdx - rho[1] * c[1] * dudx );

    const dreal MachOut = avgu / avgc;
    dreal dp;
    //dreal zero = 0.0;
    //dreal sonic = MachOut-1.0;
    //dreal avgRi = 0.5*(Ri[1] + Ri[2]);
    //condassign(dp, sonic, avgRi, zero);
    if (MachOut > 1) {
        dp = 0.5 * (Ri[1] + Ri[2]);
    } else {
        dp = 0;
    }

    const dreal drho = Ri[0] + dp / (pow(c[1], 2));
    const dreal du = (Ri[1] - dp) / (rho[1] * c[1]);

    u[1] = u[1] + du;
    rho[1] = rho[1] + drho;
    p[1] = p[1] + dp;
    const dreal T = p[1] / (rho[1] * flo_opts.R);
    e[1] = rho[1] * (flo_opts.Cv * T + 0.5 * pow(u[1], 2));

    flow_data->residual[(n_elem_ghost - 1) * 3 + 0] = (flow_data->W[(n_elem_ghost - 1) * 3 + 0] - rho[1]) / dtdx;
    flow_data->residual[(n_elem_ghost - 1) * 3 + 1] = (flow_data->W[(n_elem_ghost - 1) * 3 + 1] - rho[1] * u[1]) / dtdx;
    flow_data->residual[(n_elem_ghost - 1) * 3 + 2] = (flow_data->W[(n_elem_ghost - 1) * 3 + 2] - e[1]) / dtdx;

    flow_data->W[(n_elem_ghost - 1) * 3 + 0] = rho[1];
    flow_data->W[(n_elem_ghost - 1) * 3 + 1] = rho[1] * u[1];
    flow_data->W[(n_elem_ghost - 1) * 3 + 2] = e[1];
}
template void outletBC( const Flow_options &flo_opts, const double dt0, const double dx0, class Flow_data<double>* const flow_data);
template void outletBC( const Flow_options &flo_opts, const adouble dt0, const adouble dx0, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void outletBC(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data)
{
	//const int n_elem = flo_opts.n_elem;
	//const int n_elem_ghost = flo_opts.n_elem+2;
	const int n_elem_ghost = flow_data->W.size()/3;
	const dreal gam = flo_opts.gam;
    dreal eigenvalues[3], Ri[3];
    dreal rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 0];
        u[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 1] / rho[i];
        e[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    const dreal avgu = (u[1] + u[0]) / 2;
    const dreal avgc = (c[1] + c[0]) / 2;
    const dreal dtdx = dt0 / dx0;
    eigenvalues[0] = avgu * dtdx;
    eigenvalues[1] = (avgu + avgc) * dtdx;
    eigenvalues[2] = (avgu - avgc) * dtdx;

    const dreal dpdx = p[1] - p[0];
    const dreal dudx = u[1] - u[0];

    Ri[0] = -eigenvalues[0] * ( rho[1] - rho[0] - dpdx / pow(c[1], 2) );
    Ri[1] = -eigenvalues[1] * ( dpdx + rho[1] * c[1] * dudx );
    Ri[2] = -eigenvalues[2] * ( dpdx - rho[1] * c[1] * dudx );

    const dreal MachOut = avgu / avgc;
    dreal dp;
    //dreal zero = 0.0;
    //dreal sonic = MachOut-1.0;
    //dreal avgRi = 0.5*(Ri[1] + Ri[2]);
    //condassign(dp, sonic, avgRi, zero);
    if (MachOut > 1) {
        dp = 0.5 * (Ri[1] + Ri[2]);
    } else {
        dp = 0;
    }

    const dreal drho = Ri[0] + dp / (pow(c[1], 2));
    const dreal du = (Ri[1] - dp) / (rho[1] * c[1]);

    u[1] = u[1] + du;
    rho[1] = rho[1] + drho;
    p[1] = p[1] + dp;
    const dreal T = p[1] / (rho[1] * flo_opts.R);
    e[1] = rho[1] * (flo_opts.Cv * T + 0.5 * pow(u[1], 2));

    flow_data->residual[(n_elem_ghost - 1) * 3 + 0] = (flow_data->W[(n_elem_ghost - 1) * 3 + 0] - rho[1]) / dtdx;
    flow_data->residual[(n_elem_ghost - 1) * 3 + 1] = (flow_data->W[(n_elem_ghost - 1) * 3 + 1] - rho[1] * u[1]) / dtdx;
    flow_data->residual[(n_elem_ghost - 1) * 3 + 2] = (flow_data->W[(n_elem_ghost - 1) * 3 + 2] - e[1]) / dtdx;

    flow_data->W[(n_elem_ghost - 1) * 3 + 0] = rho[1];
    flow_data->W[(n_elem_ghost - 1) * 3 + 1] = rho[1] * u[1];
    flow_data->W[(n_elem_ghost - 1) * 3 + 2] = e[1];
}
template void outletBC( const Flow_options &flo_opts, const double dt0, const double dx0, class Flow_data<double>* const flow_data);
template void outletBC( const Flow_options &flo_opts, const adouble dt0, const adouble dx0, class Flow_data<adouble>* const flow_data);
