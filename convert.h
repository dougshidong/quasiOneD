#ifndef CONVERT_H
#define CONVERT_H

#include<vector>
#include<Eigen/Core>

//inline double get_p(const double rho, const double rho_u, const double e);
//inline double get_c(const double rho, const double rho_u, const double e);
// Get pressure
inline double get_p(const double gamma, const double rho, const double rho_u, const double e) {
	return (gamma - 1.0) * ( e - (pow(rho_u, 2.0) / rho+1e-14) / 2.0 );
};
inline double get_c(const double gamma, const double rho, const double rho_u, const double e) {
	double a =  gamma / rho * get_p(gamma, rho, rho_u, e);
	if(a<0.0) abort();
	return sqrt(a);
	//return sqrt( gamma / rho * get_p(gamma, rho, rho_u, e) );
};
void get_all_p(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &p);

void WtoP(
	const double gam,
	const double R,
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e);

void WtoP2(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &p);

void WtoP(
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &e,
    std::vector<double> &u,
    std::vector<double> &p,
    std::vector<double> &c,
    std::vector<double> &T);

void PtoW(
	const double gam,
    std::vector<double> &W,
    std::vector<double> const &Wp);

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

void WtoF(
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
#endif
