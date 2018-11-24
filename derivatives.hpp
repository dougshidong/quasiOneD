#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include"convert.hpp"
#include<adolc/adolc.h>
#include"adolc_eigen.hpp"
#include<Eigen/Core>

// p
template<typename dreal>
inline dreal get_dpdw1(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return (gamma - 1.0) * 0.5*rho_u*rho_u/(rho*rho);
};
template double get_dpdw1(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dpdw1(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dpdw2(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return -(gamma - 1.0) * rho_u/rho;
};
template double get_dpdw2(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dpdw2(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dpdw3(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return (gamma - 1.0);
};
template double get_dpdw3(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dpdw3(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

// u = rho_u/rho
template<typename dreal>
inline dreal get_dudw1(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return -rho_u/(rho*rho);
};
template double get_dudw1(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dudw1(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dudw2(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return 1.0/rho;
};
template double get_dudw2(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dudw2(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dudw3(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	dreal a = 0.0;
	return a;
};
template double get_dudw3(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dudw3(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

// c sqrt(gamma / rho * p);
template<typename dreal>
inline dreal get_dcdw1(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	const dreal c = get_c(gamma, rho, rho_u, e);
	const dreal p = get_p(gamma, rho, rho_u, e);
	const dreal dcdr = -0.5*c/rho;
	const dreal dcdp = 0.5*c/p;
	const dreal dpdr = get_dpdw1(gamma, rho, rho_u, e);
	const dreal dcdw1 = dcdr + dcdp*dpdr;
	return dcdw1;
};
template double get_dcdw1(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dcdw1(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dcdw2(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	const dreal c = get_c(gamma, rho, rho_u, e);
	const dreal p = get_p(gamma, rho, rho_u, e);
	const dreal dcdp = 0.5*c/p;
	const dreal dpdru = get_dpdw2(gamma, rho, rho_u, e);
	const dreal dcdw2 = dcdp*dpdru;
	return dcdw2;
};
template double get_dcdw2(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dcdw2(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline dreal get_dcdw3(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	const dreal c = get_c(gamma, rho, rho_u, e);
	const dreal p = get_p(gamma, rho, rho_u, e);
	const dreal dcdp = 0.5*c/p;
	const dreal dpde = get_dpdw3(gamma, rho, rho_u, e);
	const dreal dcdw3 = dcdp*dpde;
	return dcdw3;
};
template double get_dcdw3(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_dcdw3(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

template<typename dreal>
inline Eigen::Matrix<dreal, 3, 3> get_d_lambdahalfw_dw(
	const double gamma, int plus1,
	const dreal rho1, const dreal rho_u1, const dreal e1,
	const dreal rho2, const dreal rho_u2, const dreal e2)
{
	Eigen::Matrix<dreal, 3, 3> dFdW;
	// Derivative of lambda_{i+1/2} (W_{i+1} - W_{i}))
	const dreal u1 = rho_u1/rho1;
	const dreal u2 = rho_u2/rho2;
	const dreal c1 = get_c(gamma, rho1, rho_u1, e1);
	const dreal c2 = get_c(gamma, rho2, rho_u2, e2);
	const dreal lambda = 0.5*(u1+c1+u2+c2);

	const dreal delta_w1 = rho2-rho1;
	const dreal delta_w2 = rho_u2-rho_u1;
	const dreal delta_w3 = e2-e1;
	if (plus1 == 0) { // d/dW_{i}
		const dreal dldw1 = 0.5*(get_dudw1(gamma, rho1, rho_u1, e1) + get_dcdw1(gamma, rho1, rho_u1, e1));
		const dreal dldw2 = 0.5*(get_dudw2(gamma, rho1, rho_u1, e1) + get_dcdw2(gamma, rho1, rho_u1, e1));
		const dreal dldw3 = 0.5*(get_dudw3(gamma, rho1, rho_u1, e1) + get_dcdw3(gamma, rho1, rho_u1, e1));

		dFdW(0,0) = dldw1*delta_w1 - lambda;
		dFdW(0,1) = dldw2*delta_w1;
		dFdW(0,2) = dldw3*delta_w1;
		dFdW(1,0) = dldw1*delta_w2;
		dFdW(1,1) = dldw2*delta_w2 - lambda;
		dFdW(1,2) = dldw3*delta_w2;
		dFdW(2,0) = dldw1*delta_w3;
		dFdW(2,1) = dldw2*delta_w3;
		dFdW(2,2) = dldw3*delta_w3 - lambda;

	} else { // d/dW_{i+1}
		const dreal dldw1 = 0.5*(get_dudw1(gamma, rho2, rho_u2, e2) + get_dcdw1(gamma, rho2, rho_u2, e2));
		const dreal dldw2 = 0.5*(get_dudw2(gamma, rho2, rho_u2, e2) + get_dcdw2(gamma, rho2, rho_u2, e2));
		const dreal dldw3 = 0.5*(get_dudw3(gamma, rho2, rho_u2, e2) + get_dcdw3(gamma, rho2, rho_u2, e2));

		dFdW(0,0) = dldw1*delta_w1 + lambda;
		dFdW(0,1) = dldw2*delta_w1;
		dFdW(0,2) = dldw3*delta_w1;
		dFdW(1,0) = dldw1*delta_w2;
		dFdW(1,1) = dldw2*delta_w2 + lambda;
		dFdW(1,2) = dldw3*delta_w2;
		dFdW(2,0) = dldw1*delta_w3;
		dFdW(2,1) = dldw2*delta_w3;
		dFdW(2,2) = dldw3*delta_w3 + lambda;
	}
	return dFdW;
};
template Eigen::Matrix<double, 3, 3> get_d_lambdahalfw_dw(
	const double gamma, int plus1,
	const double rho1, const double rho_u1, const double e1,
	const double rho2, const double rho_u2, const double e2);
template Eigen::Matrix<adouble, 3, 3> get_d_lambdahalfw_dw(
	const double gamma, int plus1,
	const adouble rho1, const adouble rho_u1, const adouble e1,
	const adouble rho2, const adouble rho_u2, const adouble e2);

#endif
