#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include"convert.hpp"
#include<adolc/adolc.h>

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

// u
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

// c
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
#endif
