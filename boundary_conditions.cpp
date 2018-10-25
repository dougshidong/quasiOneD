#include"boundary_conditions.hpp"
#include "structures.hpp"
#include "convert.hpp"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include<adolc/adolc.h>
#include <complex>

template<typename dreal>
void inletBC(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data)
{
	const double gam = flo_opts.gam;
    dreal rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[i*3+0];
        u[i] = flow_data->W[i*3+1] / rho[i];
        e[i] = flow_data->W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
    if (u[1] < c[1]) {
		const dreal normal = 1.0;
		const dreal U_i    = u[1] * normal;
		//const dreal total_enthalpy = p[1]/rho[1] * (gam/(gam-1.0)) + 0.5*u[1]*u[1];
		const dreal total_enthalpy = (e[1] + p[1])/rho[1];
		//const dreal riemman_plus   = -U_i - 2.0*c[1]/(gam-1.0);
		const dreal riemman_plus   = U_i + 2.0*c[1]/(gam-1.0);

		const dreal a = 1.0+2.0/(gam-1.0);
		const dreal b = -2.0*riemman_plus;
		const dreal c = 0.5*(gam-1.0) * (riemman_plus*riemman_plus - 2.0*total_enthalpy);
		
		const dreal term1 = -0.5*b/a;
		const dreal term2= 0.5*sqrt(b*b-4.0*a*c)/a;
		const dreal c_b1 = term1 + term2;
		const dreal c_b2 = term1 - term2;
		const dreal c_b  = fmax(c_b1, c_b2);

		const dreal U  = riemman_plus - 2.0*c_b/(gam-1.0);
		const dreal M_b = U/c_b;

		const dreal T_b = flo_opts.inlet_total_T / (1.0+0.5*(gam-1.0)*M_b*M_b);
		const dreal p_b = flo_opts.inlet_total_p * pow(T_b / flo_opts.inlet_total_T, gam/(gam-1.0));
		//const dreal p_b = flo_opts.inlet_total_p * pow(1.0 + 0.5*(gam-1.0) * M_b*M_b, gam/(gam-1.0));
		//const dreal T_b = flo_opts.inlet_total_T * pow(p_b / flo_opts.inlet_total_p, (gam-1.0)/gam);

		const dreal rho_b = p_b/(flo_opts.R * T_b);
		const dreal u_b   = U;//*normal;

        rho[0] = rho_b;
        u[0] = u_b;
        p[0] = p_b;

		e[0] = p[0]/(gam - 1.0) + rho[0] * u[0] * u[0] / 2.0;

        flow_data->W[0 * 3 + 0] = rho[0];
        flow_data->W[0 * 3 + 1] = rho[0] * u[0];
        flow_data->W[0 * 3 + 2] = e[0];
    } else {
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
    }
}
template void inletBC(const Flow_options &flo_opts, class Flow_data<double>* const flow_data);
template void inletBC(const Flow_options &flo_opts, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void outletBC(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data)
{
	//const int n_elem = flo_opts.n_elem;
	//const int n_elem_ghost = flo_opts.n_elem+2;
	const int n_elem_ghost = flow_data->W.size()/3;
	const double gam = flo_opts.gam;
    dreal rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 0];
        u[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 1] / rho[i];
        e[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    const dreal avgu = (u[1] + u[0]) / 2.0;
    const dreal avgc = (c[1] + c[0]) / 2.0;
    const dreal MachOut = avgu / avgc;
    if (MachOut > 1.0) {
		flow_data->W[(n_elem_ghost - 1) * 3 + 0] = flow_data->W[(n_elem_ghost - 2) * 3 + 0];
		flow_data->W[(n_elem_ghost - 1) * 3 + 1] = flow_data->W[(n_elem_ghost - 2) * 3 + 1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 2] = flow_data->W[(n_elem_ghost - 2) * 3 + 2];
		return;
    } else {
		const dreal T_domain = p[0] / (rho[0] * flo_opts.R); // Extrapolate
        rho[1] = gam * flo_opts.outlet_p / T_domain;
        u[1] = u[0]; // Extrapolate
        p[1] = flo_opts.outlet_p; // Specify
		e[1] = p[1]/(gam - 1) + rho[1] * u[1] * u[1] / 2.0;

		flow_data->W[(n_elem_ghost - 1) * 3 + 0] = rho[1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 1] = rho[1] * u[1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 2] = e[1];
    }

	return;
}
template void outletBC( const Flow_options &flo_opts, class Flow_data<double>* const flow_data);
template void outletBC( const Flow_options &flo_opts, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void outletBC_old(
	const Flow_options &flo_opts,
	const dreal dt0,
	const dreal dx0,
	class Flow_data<dreal>* const flow_data)
{
	//const int n_elem = flo_opts.n_elem;
	//const int n_elem_ghost = flo_opts.n_elem+2;
	const int n_elem_ghost = flow_data->W.size()/3;
	const double gam = flo_opts.gam;
	dreal eigenvalues[3], Ri[3];
	dreal rho[2], u[2], e[2], p[2], c[2];

	for (int i = 0; i < 2; i++) {
		rho[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 0];
		u[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 1] / rho[i];
		e[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 2];
		p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
		c[i] = sqrt( gam * p[i] / rho[i] );
	}

	// Exit boundary condition
	const dreal avgu = (u[1] + u[0]) / 2.0;
	const dreal avgc = (c[1] + c[0]) / 2.0;
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

		//flow_data->W[(n_elem_ghost - 1) * 3 + 0] = flow_data->W[(n_elem_ghost - 2) * 3 + 0];
		//flow_data->W[(n_elem_ghost - 1) * 3 + 1] = flow_data->W[(n_elem_ghost - 2) * 3 + 1];
		//flow_data->W[(n_elem_ghost - 1) * 3 + 2] = flow_data->W[(n_elem_ghost - 2) * 3 + 2];
		//return;
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
template void outletBC_old( const Flow_options &flo_opts, const double dt0, const double dx0, class Flow_data<double>* const flow_data);
template void outletBC_old( const Flow_options &flo_opts, const adouble dt0, const adouble dx0, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void inletBC_old(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data)
{
	const dreal a2  = flo_opts.a2;
	const double gam = flo_opts.gam;
    dreal rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[i*3+0];
        u[i] = flow_data->W[i*3+1] / rho[i];
        e[i] = flow_data->W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
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
template void inletBC_old(const Flow_options &flo_opts, const double dt0, const double dx0, class Flow_data<double>* const flow_data);
template void inletBC_old(const Flow_options &flo_opts, const adouble dt0, const adouble dx0, class Flow_data<adouble>* const flow_data);

template<>
void inletBC<std::complex<double>>(
    const Flow_options &flo_opts,
    class Flow_data<std::complex<double>>* const flow_data)
{
	const double gam = flo_opts.gam;
    std::complex<double> rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[i*3+0];
        u[i] = flow_data->W[i*3+1] / rho[i];
        e[i] = flow_data->W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
    if (real(u[1]) < real(c[1])) {
		const std::complex<double> normal = 1.0;
		const std::complex<double> U_i    = u[1] * normal;
		//const std::complex<double> total_enthalpy = p[1]/rho[1] * (gam/(gam-1.0)) + 0.5*u[1]*u[1];
		const std::complex<double> total_enthalpy = (e[1] + p[1])/rho[1];
		//const std::complex<double> riemman_plus   = -U_i - 2.0*c[1]/(gam-1.0);
		const std::complex<double> riemman_plus   = U_i + 2.0*c[1]/(gam-1.0);

		const std::complex<double> a = 1.0+2.0/(gam-1.0);
		const std::complex<double> b = -2.0*riemman_plus;
		const std::complex<double> c = 0.5*(gam-1.0) * (riemman_plus*riemman_plus - 2.0*total_enthalpy);
		
		const std::complex<double> term1 = -0.5*b/a;
		const std::complex<double> term2= 0.5*sqrt(b*b-4.0*a*c)/a;
		const std::complex<double> c_b1 = term1 + term2;
		const std::complex<double> c_b2 = term1 - term2;

		std::complex<double> c_b;
		if (real(c_b1) > real(c_b2)) {
			c_b  = c_b1;
		} else {
			c_b  = c_b2;
		}

		const std::complex<double> U  = riemman_plus - 2.0*c_b/(gam-1.0);
		const std::complex<double> M_b = U/c_b;

		const std::complex<double> T_b = flo_opts.inlet_total_T / (1.0+0.5*(gam-1.0)*M_b*M_b);
		const std::complex<double> p_b = flo_opts.inlet_total_p * pow(T_b / flo_opts.inlet_total_T, gam/(gam-1.0));
		//const std::complex<double> p_b = flo_opts.inlet_total_p * pow(1.0 + 0.5*(gam-1.0) * M_b*M_b, gam/(gam-1.0));
		//const std::complex<double> T_b = flo_opts.inlet_total_T * pow(p_b / flo_opts.inlet_total_p, (gam-1.0)/gam);

		const std::complex<double> rho_b = p_b/(flo_opts.R * T_b);
		const std::complex<double> u_b   = U;//*normal;

        rho[0] = rho_b;
        u[0] = u_b;
        p[0] = p_b;

		e[0] = p[0]/(gam - 1.0) + rho[0] * u[0] * u[0] / 2.0;

        flow_data->W[0 * 3 + 0] = rho[0];
        flow_data->W[0 * 3 + 1] = rho[0] * u[0];
        flow_data->W[0 * 3 + 2] = e[0];
    } else {
		std::complex<double> inlet_T = isenT(gam, flo_opts.inlet_total_T, flo_opts.inlet_mach);
		const std::complex<double> p_inlet = isenP(gam, flo_opts.inlet_total_p, flo_opts.inlet_mach);
		std::complex<double> p = p_inlet;
		std::complex<double> rho = p / (flo_opts.R * inlet_T);
		std::complex<double> c = sqrt(gam * p / rho);
		std::complex<double> u = flo_opts.inlet_mach * c;
		std::complex<double> e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

		flow_data->W[0*3+0] = rho;
		flow_data->W[0*3+1] = rho * u;
		flow_data->W[0*3+2] = e;
    }
}

template<>
void outletBC<std::complex<double>>(
    const Flow_options &flo_opts,
    class Flow_data<std::complex<double>>* const flow_data)
{
	//const int n_elem = flo_opts.n_elem;
	//const int n_elem_ghost = flo_opts.n_elem+2;
	const int n_elem_ghost = flow_data->W.size()/3;
	const double gam = flo_opts.gam;
    std::complex<double> rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 0];
        u[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 1] / rho[i];
        e[i] = flow_data->W[(i + (n_elem_ghost - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    const std::complex<double> avgu = (u[1] + u[0]) / 2.0;
    const std::complex<double> avgc = (c[1] + c[0]) / 2.0;
    const std::complex<double> MachOut = avgu / avgc;
    if (real(MachOut) > 1.0) {
		flow_data->W[(n_elem_ghost - 1) * 3 + 0] = flow_data->W[(n_elem_ghost - 2) * 3 + 0];
		flow_data->W[(n_elem_ghost - 1) * 3 + 1] = flow_data->W[(n_elem_ghost - 2) * 3 + 1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 2] = flow_data->W[(n_elem_ghost - 2) * 3 + 2];
		return;
    } else {
		const std::complex<double> T_domain = p[0] / (rho[0] * flo_opts.R); // Extrapolate
        rho[1] = gam * flo_opts.outlet_p / T_domain;
        u[1] = u[0]; // Extrapolate
        p[1] = flo_opts.outlet_p; // Specify
		e[1] = p[1]/(gam - 1) + rho[1] * u[1] * u[1] / 2.0;

		flow_data->W[(n_elem_ghost - 1) * 3 + 0] = rho[1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 1] = rho[1] * u[1];
		flow_data->W[(n_elem_ghost - 1) * 3 + 2] = e[1];
    }

	return;
}
