#include"boundary_conditions.hpp"
#include "structures.hpp"
#include "convert.hpp"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include<Eigen/Core>

#include<adolc/adolc.h>
#include <complex>


template<typename dreal>
Eigen::Matrix<dreal,3,3> inletBC_gradient(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data)
{
	Eigen::Matrix<dreal,3,3> dWbdWd;
	const double gam = flo_opts.gam;
	const int i_d = 1;
	const dreal w1 = flow_data->W[i_d*3+0];
	const dreal w2 = flow_data->W[i_d*3+1];

	const dreal rho_d = flow_data->W[i_d*3+0];
	const dreal u_d   = flow_data->W[i_d*3+1] / rho_d;
	const dreal e_d   = flow_data->W[i_d*3+2];
	const dreal p_d   = (gam - 1) * ( e_d - rho_d * u_d * u_d / 2.0 );
	const dreal c_d   = sqrt( gam * p_d / rho_d );

	//const dreal d_rho_d_dw1 = 1.0;
	//const dreal d_rho_d_dw2 = 0.0;
	//const dreal d_rho_d_dw3 = 0.0;
	const dreal d_u_d_dw1 = -flow_data->W[i_d*3+1] / (rho_d*rho_d);
	const dreal d_u_d_dw2 = 1.0 / rho_d;
	const dreal d_u_d_dw3 = 0.0;
	//const dreal d_e_d_dw1 = 0.0;
	//const dreal d_e_d_dw2 = 0.0;
	const dreal d_e_d_dw3 = 1.0;
	const dreal d_p_d_dw1 = (gam - 1.0) * 0.5*w2*w2/(w1*w1);
	const dreal d_p_d_dw2 = -(gam - 1.0) * w2/w1;
	const dreal d_p_d_dw3 = (gam - 1.0);

	const dreal d_c_d_dp  = 0.5*sqrt(gam / (w1*p_d));
	const dreal d_c_d_dw1 = -0.5*sqrt(gam * p_d / (rho_d*rho_d*rho_d)) + d_c_d_dp * d_p_d_dw1;
	const dreal d_c_d_dw2 = d_c_d_dp * d_p_d_dw2;
	const dreal d_c_d_dw3 = d_c_d_dp * d_p_d_dw3;

    if (u_d < c_d) {
		const dreal normal = 1.0;
		const dreal U_i    = u_d * normal;
		
		const dreal d_U_i_dw1 = d_u_d_dw1 * normal;
		const dreal d_U_i_dw2 = d_u_d_dw2 * normal;
		const dreal d_U_i_dw3 = d_u_d_dw3 * normal;

		const dreal total_enthalpy = (e_d + p_d)/rho_d;
		const dreal d_total_enthalpy_dw1 = -(e_d + p_d)/(rho_d*rho_d) + d_p_d_dw1/rho_d;
		const dreal d_total_enthalpy_dw2 = d_p_d_dw2/rho_d;
		const dreal d_total_enthalpy_dw3 = (d_e_d_dw3+d_p_d_dw3)/rho_d;

		const dreal riemman_plus   = U_i + 2.0*c_d/(gam-1.0);

		const dreal d_riemman_plus_dw1 = d_U_i_dw1 + 2.0/(gam-1.0) * d_c_d_dw1;
		const dreal d_riemman_plus_dw2 = d_U_i_dw2 + 2.0/(gam-1.0) * d_c_d_dw2;
		const dreal d_riemman_plus_dw3 = d_U_i_dw3 + 2.0/(gam-1.0) * d_c_d_dw3;

		const dreal a = 1.0+2.0/(gam-1.0);
		const dreal b = -2.0*riemman_plus;
		const dreal c = 0.5*(gam-1.0) * (riemman_plus*riemman_plus - 2.0*total_enthalpy);

		const dreal d_b_dw1 =  -2.0*d_riemman_plus_dw1;
		const dreal d_b_dw2 =  -2.0*d_riemman_plus_dw2;
		const dreal d_b_dw3 =  -2.0*d_riemman_plus_dw3;

		const dreal d_c_dw1 =  0.5*(gam-1.0) * (2.0*riemman_plus + d_riemman_plus_dw1 - 2.0*d_total_enthalpy_dw1);
		const dreal d_c_dw2 =  0.5*(gam-1.0) * (2.0*riemman_plus + d_riemman_plus_dw2 - 2.0*d_total_enthalpy_dw2);
		const dreal d_c_dw3 =  0.5*(gam-1.0) * (2.0*riemman_plus + d_riemman_plus_dw3 - 2.0*d_total_enthalpy_dw3);
		
		const dreal term1 = -0.5*b/a;
		const dreal term2= 0.5*sqrt(b*b-4.0*a*c)/a;

		const dreal d_term1_dw1 = -0.5*d_b_dw1/a;
		const dreal d_term1_dw2 = -0.5*d_b_dw2/a;
		const dreal d_term1_dw3 = -0.5*d_b_dw3/a;

		const dreal d_term2_db = 0.5*b/sqrt(b*b-4.0*a*c)/a;
		const dreal d_term2_dc = -1.0/sqrt(b*b-4.0*a*c);

		const dreal d_term2_dw1 = d_term2_db*d_b_dw1 + d_term2_dc*d_c_dw1;
		const dreal d_term2_dw2 = d_term2_db*d_b_dw2 + d_term2_dc*d_c_dw2;
		const dreal d_term2_dw3 = d_term2_db*d_b_dw3 + d_term2_dc*d_c_dw3;

		const dreal c_b1 = term1 + term2;
		const dreal c_b2 = term1 - term2;

		const dreal c_b  = fmax(c_b1, c_b2);
		
		dreal d_c_b_dw1, d_c_b_dw2, d_c_b_dw3;
		// condassign here******************************
		const dreal dcb = c_b1 - c_b2;
		condassign(d_c_b_dw1, dcb, d_term1_dw1 + d_term2_dw1, d_term1_dw1 - d_term2_dw1);
		condassign(d_c_b_dw2, dcb, d_term1_dw2 + d_term2_dw2, d_term1_dw2 - d_term2_dw2);
		condassign(d_c_b_dw3, dcb, d_term1_dw3 + d_term2_dw3, d_term1_dw3 - d_term2_dw3);

		const dreal U  = riemman_plus - 2.0*c_b/(gam-1.0);

		const dreal d_U_dw1 = d_riemman_plus_dw1 - 2.0/(gam-1.0) * d_c_b_dw1;
		const dreal d_U_dw2 = d_riemman_plus_dw2 - 2.0/(gam-1.0) * d_c_b_dw2;
		const dreal d_U_dw3 = d_riemman_plus_dw3 - 2.0/(gam-1.0) * d_c_b_dw3;

		const dreal M_b = U/c_b;

		const dreal d_M_b_dw1 = d_U_dw1 / c_b - U/(2.0*c_b*c_b) * d_c_b_dw1;
		const dreal d_M_b_dw2 = d_U_dw2 / c_b - U/(2.0*c_b*c_b) * d_c_b_dw2;
		const dreal d_M_b_dw3 = d_U_dw3 / c_b - U/(2.0*c_b*c_b) * d_c_b_dw3;
		

		const dreal T_b = flo_opts.inlet_total_T / (1.0+0.5*(gam-1.0)*M_b*M_b);

		const dreal d_T_b_dM_b = flo_opts.inlet_total_T * (gam-1.0) * M_b / pow((1.0+0.5*(gam-1.0)*M_b*M_b),2);
		const dreal d_T_b_dw1 = d_T_b_dM_b * d_M_b_dw1;
		const dreal d_T_b_dw2 = d_T_b_dM_b * d_M_b_dw2;
		const dreal d_T_b_dw3 = d_T_b_dM_b * d_M_b_dw3;

		const dreal p_b = flo_opts.inlet_total_p * pow(T_b / flo_opts.inlet_total_T, gam/(gam-1.0));

		const dreal d_p_b_dT_b = flo_opts.inlet_total_p * gam/(gam-1.0) / T_b * pow(T_b / flo_opts.inlet_total_T, gam/(gam-1.0));
		const dreal d_p_b_dw1 = d_p_b_dT_b * d_T_b_dw1;
		const dreal d_p_b_dw2 = d_p_b_dT_b * d_T_b_dw2;
		const dreal d_p_b_dw3 = d_p_b_dT_b * d_T_b_dw3;

		const dreal rho_b = p_b/(flo_opts.R * T_b);
		
		const dreal d_rho_b_dp_b = 1.0/(flo_opts.R * T_b);
		const dreal d_rho_b_dT_b = -p_b/(flo_opts.R * T_b*T_b);

		const dreal d_rho_b_dw1 = d_rho_b_dp_b*d_p_b_dw1 + d_rho_b_dT_b*d_T_b_dw1;
		const dreal d_rho_b_dw2 = d_rho_b_dp_b*d_p_b_dw2 + d_rho_b_dT_b*d_T_b_dw2;
		const dreal d_rho_b_dw3 = d_rho_b_dp_b*d_p_b_dw3 + d_rho_b_dT_b*d_T_b_dw3;

		const dreal u_b   = U;//*normal;

		const dreal d_u_b_dw1 = d_U_dw1;
		const dreal d_u_b_dw2 = d_U_dw2;
		const dreal d_u_b_dw3 = d_U_dw3;

		//const dreal e_b = p_b/(gam - 1.0) + rho_b * u_b * u_b / 2.0;

		const dreal d_e_b_d_p_b = 1.0/(gam - 1.0);
		const dreal d_e_b_d_rho_b = u_b * u_b / 2.0;
		const dreal d_e_b_d_u_b = rho_b * u_b;

		const dreal d_e_b_dw1 = d_e_b_d_p_b   * d_p_b_dw1 
					     	 + d_e_b_d_rho_b * d_rho_b_dw1
					     	 + d_e_b_d_u_b   * d_u_b_dw1;
		const dreal d_e_b_dw2 = d_e_b_d_p_b   * d_p_b_dw2 
					     	 + d_e_b_d_rho_b * d_rho_b_dw2
					     	 + d_e_b_d_u_b   * d_u_b_dw2;
		const dreal d_e_b_dw3 = d_e_b_d_p_b   * d_p_b_dw3 
							 + d_e_b_d_rho_b * d_rho_b_dw3
							 + d_e_b_d_u_b   * d_u_b_dw3;


        //flow_data->W[0 * 3 + 0] = rho_b;
        //flow_data->W[0 * 3 + 1] = rho_b * u_b;
        //flow_data->W[0 * 3 + 2] = e_b;

		dWbdWd(0,0) = d_rho_b_dw1;
		dWbdWd(0,1) = d_rho_b_dw2;
		dWbdWd(0,2) = d_rho_b_dw3;
                 
		dWbdWd(1,0) = d_rho_b_dw1 * u_b + rho_b * d_u_b_dw1;
		dWbdWd(1,1) = d_rho_b_dw2 * u_b + rho_b * d_u_b_dw2;
		dWbdWd(1,2) = d_rho_b_dw3 * u_b + rho_b * d_u_b_dw3;
                 
		dWbdWd(2,0) = d_e_b_dw1;
		dWbdWd(2,1) = d_e_b_dw2;
		dWbdWd(2,2) = d_e_b_dw3;

		//const dreal dWb1dWd1 = d_rho_b_dw1;
		//const dreal dWb1dWd2 = d_rho_b_dw2;
		//const dreal dWb1dWd3 = d_rho_b_dw3;

		//const dreal dWb2dWd1 = d_rho_b_dw1 * u_b + rho_b * d_u_b_dw1;
		//const dreal dWb2dWd2 = d_rho_b_dw2 * u_b + rho_b * d_u_b_dw2;
		//const dreal dWb2dWd3 = d_rho_b_dw3 * u_b + rho_b * d_u_b_dw3;

		//const dreal dWb3dWd1 = d_e_b_dw1;
		//const dreal dWb3dWd2 = d_e_b_dw2;
		//const dreal dWb3dWd3 = d_e_b_dw3;


    } else {
		const dreal inlet_T = isenT(gam, flo_opts.inlet_total_T, flo_opts.inlet_mach);
		const dreal p_inlet = isenP(gam, flo_opts.inlet_total_p, flo_opts.inlet_mach);
		const dreal p = p_inlet;
		const dreal rho = p / (flo_opts.R * inlet_T);
		const dreal c = sqrt(gam * p / rho);
		const dreal u = flo_opts.inlet_mach * c;
		const dreal e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

		flow_data->W[0*3+0] = rho;
		flow_data->W[0*3+1] = rho * u;
		flow_data->W[0*3+2] = e;

		//const dreal dWb1dWd1 = 0.0;
		//const dreal dWb1dWd2 = 0.0;
		//const dreal dWb1dWd3 = 0.0;

		//const dreal dWb2dWd1 = 0.0;
		//const dreal dWb2dWd2 = 0.0;
		//const dreal dWb2dWd3 = 0.0;

		//const dreal dWb3dWd1 = 0.0;
		//const dreal dWb3dWd2 = 0.0;
		//const dreal dWb3dWd3 = 0.0;

		dWbdWd(0,0) = 0.0; 
		dWbdWd(0,1) = 0.0; 
		dWbdWd(0,2) = 0.0; 
                 
		dWbdWd(1,0) = 0.0; 
		dWbdWd(1,1) = 0.0; 
		dWbdWd(1,2) = 0.0; 
                 
		dWbdWd(2,0) = 0.0; 
		dWbdWd(2,1) = 0.0; 
		dWbdWd(2,2) = 0.0; 
    }
	return dWbdWd;
}
template Eigen::Matrix<double,3,3> inletBC_gradient(const Flow_options &flo_opts, class Flow_data<double>* const flow_data);
template Eigen::Matrix<adouble,3,3> inletBC_gradient(const Flow_options &flo_opts, class Flow_data<adouble>* const flow_data);

template<typename dreal>
Eigen::Matrix<dreal,3,3> outletBC_gradient(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data)
{
	Eigen::Matrix<dreal,3,3> dWbdWd;
	const int n_elem = flo_opts.n_elem;
	const double gam = flo_opts.gam;
	//const int i_b = n_elem+1;
	const int i_d = n_elem;

    const dreal rho_d = flow_data->W[i_d*3+0];
    const dreal u_d = flow_data->W[i_d*3+1] / rho_d;
    const dreal e_d = flow_data->W[i_d*3+2];
    const dreal p_d = (gam - 1.0) * ( e_d - rho_d * u_d * u_d / 2.0 );
    const dreal c_d = sqrt( gam * p_d / rho_d );

	const dreal w1 = flow_data->W[i_d*3+0];
	const dreal w2 = flow_data->W[i_d*3+1];
	const dreal d_rho_d_dw1 = 1.0;
	const dreal d_rho_d_dw2 = 0.0;
	const dreal d_rho_d_dw3 = 0.0;
	const dreal d_u_d_dw1 = -flow_data->W[i_d*3+1] / (rho_d*rho_d);
	const dreal d_u_d_dw2 = 1.0 / rho_d;
	const dreal d_u_d_dw3 = 0.0;
	//const dreal d_e_d_dw1 = 0.0;
	//const dreal d_e_d_dw2 = 0.0;
	//const dreal d_e_d_dw3 = 1.0;
	const dreal d_p_d_dw1 = (gam - 1.0) * 0.5*w2*w2/(w1*w1);
	const dreal d_p_d_dw2 = -(gam - 1.0) * w2/w1;
	const dreal d_p_d_dw3 = (gam - 1.0);

	//const dreal d_c_d_dp  = 0.5*sqrt(gam / (w1*p_d));
	//const dreal d_c_d_dw1 = -0.5*sqrt(gam * p_d / (rho_d*rho_d*rho_d)) + d_c_d_dp * d_p_d_dw1;
	//const dreal d_c_d_dw2 = d_c_d_dp * d_p_d_dw2;
	//const dreal d_c_d_dw3 = d_c_d_dp * d_p_d_dw3;


	const dreal T_domain = p_d / (rho_d * flo_opts.R); // Extrapolate

	const dreal d_T_domain_d_p_d = 1.0 / (rho_d * flo_opts.R);
	const dreal d_T_domain_d_rho_d = -0.5*p_d / (rho_d * rho_d * flo_opts.R);
	const dreal d_T_domain_dw1 = d_T_domain_d_p_d * d_p_d_dw1 + d_T_domain_d_rho_d * d_rho_d_dw1;
	const dreal d_T_domain_dw2 = d_T_domain_d_p_d * d_p_d_dw2 + d_T_domain_d_rho_d * d_rho_d_dw2;
	const dreal d_T_domain_dw3 = d_T_domain_d_p_d * d_p_d_dw3 + d_T_domain_d_rho_d * d_rho_d_dw3;

	const dreal rho_b = gam * flo_opts.outlet_p / T_domain;

	const dreal d_rho_b_dw1 = -gam * flo_opts.outlet_p / (T_domain*T_domain) * d_T_domain_dw1;
	const dreal d_rho_b_dw2 = -gam * flo_opts.outlet_p / (T_domain*T_domain) * d_T_domain_dw2;
	const dreal d_rho_b_dw3 = -gam * flo_opts.outlet_p / (T_domain*T_domain) * d_T_domain_dw3;

	const dreal u_b = u_d; // Extrapolate

	const dreal d_u_b_dw1 = d_u_d_dw1;
	const dreal d_u_b_dw2 = d_u_d_dw2;
	const dreal d_u_b_dw3 = d_u_d_dw3;

	//const dreal p_b = flo_opts.outlet_p; // Specify

	//const dreal e_b = p_b/(gam - 1) + rho_b * u_b * u_b / 2.0;
	
	const dreal d_e_b_d_rho_b = u_b * u_b / 2.0;
	const dreal d_e_b_d_u_b = rho_b * u_b;

	const dreal d_e_b_dw1 = d_e_b_d_rho_b * d_rho_b_dw1 + d_e_b_d_u_b   * d_u_b_dw1;
	const dreal d_e_b_dw2 = d_e_b_d_rho_b * d_rho_b_dw2 + d_e_b_d_u_b   * d_u_b_dw2;
	const dreal d_e_b_dw3 = d_e_b_d_rho_b * d_rho_b_dw3 + d_e_b_d_u_b   * d_u_b_dw3;


	const dreal mach_d = u_d + c_d;
    //if (mach_d > 1.0) {
	//	flow_data->W[i_b*3+0] = flow_data->W[i_d*3+0];
	//	flow_data->W[i_b*3+1] = flow_data->W[i_d*3+1];
	//	flow_data->W[i_b*3+2] = flow_data->W[i_d*3+2];
	//	return;
    //} else {

	//	flow_data->W[i_b*3+0] = rho_b;
	//	flow_data->W[i_b*3+1] = rho_b * u_b;
	//	flow_data->W[i_b*3+2] = e_b;
    //}
	dreal dWbdWd11,  dWbdWd12, dWbdWd13, dWbdWd21, dWbdWd22, dWbdWd23, dWbdWd31, dWbdWd32, dWbdWd33;
	const dreal one = 1.0, zero = 0.0;
	condassign(dWbdWd11, mach_d-1.0, one, d_rho_b_dw1);
	condassign(dWbdWd12, mach_d-1.0, zero, d_rho_b_dw2);
	condassign(dWbdWd13, mach_d-1.0, zero, d_rho_b_dw3);

	condassign(dWbdWd21, mach_d-1.0, zero, d_rho_b_dw1*u_b + rho_b * d_u_b_dw1);
	condassign(dWbdWd22, mach_d-1.0, one, d_rho_b_dw2*u_b + rho_b * d_u_b_dw2);
	condassign(dWbdWd23, mach_d-1.0, zero, d_rho_b_dw3*u_b + rho_b * d_u_b_dw3);

	condassign(dWbdWd31, mach_d-1.0, zero, d_e_b_dw1);
	condassign(dWbdWd32, mach_d-1.0, zero, d_e_b_dw2);
	condassign(dWbdWd33, mach_d-1.0, one, d_e_b_dw3);

	dWbdWd(0,0) = dWbdWd11;
	dWbdWd(0,1) = dWbdWd12;
	dWbdWd(0,2) = dWbdWd13;
	dWbdWd(1,0) = dWbdWd21;
	dWbdWd(1,1) = dWbdWd22;
	dWbdWd(1,2) = dWbdWd23;
	dWbdWd(2,0) = dWbdWd31;
	dWbdWd(2,1) = dWbdWd32;
	dWbdWd(2,2) = dWbdWd33;

	return dWbdWd;
}
template Eigen::Matrix<double,3,3> outletBC_gradient( const Flow_options &flo_opts, class Flow_data<double>* const flow_data);
template Eigen::Matrix<adouble,3,3> outletBC_gradient( const Flow_options &flo_opts, class Flow_data<adouble>* const flow_data);

