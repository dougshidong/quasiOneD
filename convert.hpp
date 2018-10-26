#ifndef CONVERT_H
#define CONVERT_H

#include<vector>
#include<Eigen/Core>
#include<adolc/adolc.h>

template<typename dreal>
inline dreal norm2 (std::vector<dreal> vec) {
	dreal norm_val = 0;
	for (int i = 0; i < vec.size(); i++) { norm_val = norm_val + vec[i]*vec[i]; }
	return sqrt(norm_val);
}
template double norm2 (std::vector<double> vec);
template adouble norm2 (std::vector<adouble> vec);

// Get pressure
template<typename dreal>
inline dreal get_p(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	return (gamma - 1.0) * ( e - (pow(rho_u, 2.0) / rho+1e-14) / 2.0 );
};
template double get_p(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_p(const double gamma, const adouble rho, const adouble rho_u, const adouble e);
template std::complex<double> get_p(const double gamma, const std::complex<double> rho, const std::complex<double> rho_u, const std::complex<double> e);

template<typename dreal>
inline dreal get_c(const double gamma, const dreal rho, const dreal rho_u, const dreal e) {
	const dreal p = get_p(gamma, rho, rho_u, e);
	dreal a =  gamma / rho * p;
	//if(a<0.0) {
	//	//printf("Negative radicand of sqrt. Either p: %4.2f or rho: %4.2f", p, rho);
	//	abort();
	//}
	return sqrt(a);
};
template double get_c(const double gamma, const double rho, const double rho_u, const double e);
template adouble get_c(const double gamma, const adouble rho, const adouble rho_u, const adouble e);

//template<>
//std::complex<double> get_c<std::complex<double>(const double gamma, const std::complex<double> rho, const std::complex<double> rho_u, const std::complex<double> e) {
//	const std::complex<double> p = get_p(gamma, rho, rho_u, e);
//	std::complex<double> a =  gamma / rho * p;
//	return std::sqrt(a);
//};

template<typename dreal>
inline dreal isenP(const double gam, const double pt, const dreal M) {
	return pt * pow((1 + (gam - 1) / 2 * pow(M, 2)), ( - gam / (gam - 1)));
};
template double isenP(const double gam, const double pt, const double M);
template adouble isenP(const double gam, const double pt, const adouble M);

template<typename dreal>
inline dreal isenT(const double gam, const double Tt, const dreal M) {
	return Tt * pow((1 + (gam - 1) / 2 * pow(M, 2)), - 1);
};
template double isenT(const double gam, const double Tt, const double M);
template adouble isenT(const double gam, const double Tt, const adouble M);

void get_all_p(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &p);

void eval_dWpdW(
	const double gam,
	const double rho,
	const double rho_u,
    std::vector<double>* const dWpdW_vec);

Eigen::MatrixXd eval_dWpdW(
	const double gam,
	const double rho,
	const double rho_u);

void dWdWp(
	const double gam,
    std::vector<double> &M,
    std::vector<double> const &W,
    int k);

void WtoF_all(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &F);

void WtoQ(
	const double gam,
    const std::vector<double> &area,
    const std::vector<double> &W,
    std::vector<double> &Q);

std::vector <Eigen::Matrix3d> ddWpdWdWp(
    const double gam,
    const std::vector<double> &W,
    const int i);

template<typename dreal> Eigen::Matrix<dreal, 3, 3> analytical_flux_jacobian( const double gam, const Eigen::Matrix<dreal, 3, 1> &W);
template<typename dreal> Eigen::Matrix<dreal, 3, 1> VectorToEigen3(const dreal w1, const dreal w2, const dreal w3);
template<typename dreal> Eigen::Matrix<dreal, 3, 1> WtoF( const double gam, const Eigen::Matrix<dreal, 3, 1> &W);

#endif
