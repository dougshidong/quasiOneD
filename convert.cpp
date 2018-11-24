#include"convert.hpp"
#include<vector>
#include<Eigen/Core>
#include<math.h>
#include"convert.hpp"

#include<adolc/adolc.h>
#include"adolc_eigen.hpp"

#include<complex>


// Get pressure vector
void get_all_p(
	const double gam,
	std::vector<double> const &W,
	std::vector<double> &p)
{
	int n_elem = p.size();
    assert(n_elem == W.size()/3);
    for (int i = 0; i < n_elem; i++) {
		p[i] = get_p(gam, W[i*3+0], W[i*3+1], W[i*3+2]);
    }
}

template<typename dreal>
Eigen::Matrix<dreal, 3, 3> analytical_flux_jacobian(
	const double gam,
	const Eigen::Matrix<dreal, 3, 1> &W)
{
	Eigen::Matrix<dreal, 3, 3> dFdW;
	const dreal w1 = W(0);
	const dreal w2 = W(1);
	const dreal w3 = W(2);
	dFdW(0,0) = 0.0;
	dFdW(0,1) = 1.0;
	dFdW(0,2) = 0.0;
	dFdW(1,0) = 0.5 * w2*w2 * (gam-3.0) / (w1*w1);
	dFdW(1,1) = -w2*(gam-3.0)/w1;
	dFdW(1,2) = (gam-1.0);
	dFdW(2,0) = (w2*w2*w2*(gam-1.0) - w1*w2*w3*gam) / (w1*w1*w1);
	dFdW(2,1) = w3*gam/w1 - 3.0*w2*w2*(gam-1.0) / (2.0*w1*w1);
	dFdW(2,2) = w2*gam/w1;

	return dFdW;
}
template Eigen::Matrix<double, 3, 3> analytical_flux_jacobian( const double gam, const Eigen::Matrix<double, 3, 1> &W);
template Eigen::Matrix<adouble, 3, 3> analytical_flux_jacobian( const double gam, const Eigen::Matrix<adouble, 3, 1> &W);


template<typename dreal>
Eigen::Matrix<dreal, 3, 1> VectorToEigen3(const dreal w1, const dreal w2, const dreal w3) {
	Eigen::Matrix<dreal, 3, 1> W;
	W(0) = w1;
	W(1) = w2;
	W(2) = w3;
	return W;
}
template Eigen::Matrix<double, 3, 1> VectorToEigen3(const double w1, const double w2, const double w3);
template Eigen::Matrix<adouble, 3, 1> VectorToEigen3(const adouble w1, const adouble w2, const adouble w3);

template<typename dreal>
Eigen::Matrix<dreal, 3, 1> WtoF(
	const double gam,
	const Eigen::Matrix<dreal, 3, 1> &W)
{
	const dreal w1 = W(0);
	const dreal w2 = W(1);
	const dreal w3 = W(2);

	Eigen::Matrix<dreal, 3, 1> F;
	F(0) = w2;
	F(1) = w2 * w2 / w1 + (gam - 1.0) * ( w3 - w2 * w2 / (2.0 * w1) );
	F(2) = ( w3 + (gam - 1.0) * (w3 - w2 * w2 / (2.0 * w1)) ) * w2 / w1;
	return F;
}
template Eigen::Matrix<double, 3, 1> WtoF( const double gam, const Eigen::Matrix<double, 3, 1> &W);
template Eigen::Matrix<adouble, 3, 1> WtoF( const double gam, const Eigen::Matrix<adouble, 3, 1> &W);
template Eigen::Matrix<std::complex<double>, 3, 1> WtoF( const double gam, const Eigen::Matrix<std::complex<double>, 3, 1> &W);

void WtoF_all(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &F)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
        double w1 = W[i*3+0];
        double w2 = W[i*3+1];
        double w3 = W[i*3+2];
        F[i*3+0] = w2;
        F[i*3+1] = w2 * w2 / w1 + (gam - 1.0) * ( w3 - w2 * w2 / (2.0 * w1) );
        F[i*3+2] = ( w3 + (gam - 1.0) * (w3 - w2 * w2 / (2.0 * w1)) ) * w2 / w1;
    }
}

// get Q
void WtoQ(
	const double gam,
    const std::vector<double> &area,
    const std::vector<double> &W,
    std::vector<double> &Q)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
		double p = get_p(gam,W[i*3+0], W[i*3+1], W[i*3+2]);
        Q[i*3+0] = 0;
        Q[i*3+1] = p * (area[i + 1] - area[i]);
        Q[i*3+2] = 0;
    }
}

// Inverse[dW/dWp] or dWp/dW
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void eval_dWpdW(
	const double gam,
	const double rho,
	const double rho_u,
    std::vector<double>* const dWpdW_vec)
{
    const double u = rho_u / rho;
    (*dWpdW_vec)[0] = 1.0;
    (*dWpdW_vec)[1] = 0.0;
    (*dWpdW_vec)[2] = 0.0;
    (*dWpdW_vec)[3] = -u / rho;
    (*dWpdW_vec)[4] = 1.0 / rho;
    (*dWpdW_vec)[5] = 0.0;
    (*dWpdW_vec)[6] = u * u * (gam - 1.0) / 2.0;
    (*dWpdW_vec)[7] = u * (1.0 - gam);
    (*dWpdW_vec)[8] = gam - 1.0;
}

Eigen::MatrixXd eval_dWpdW(
	const double gam,
	const double rho,
	const double rho_u)
{
    Eigen::MatrixXd M(3, 3);
    const double u = rho_u / rho;
    M(0,0) = 1.0;
    M(0,1) = 0.0;
    M(0,2) = 0.0;
    M(1,0) = -u / rho;
    M(1,1) = 1.0 / rho;
    M(1,2) = 0.0;
    M(2,0) = u * u * (gam - 1.0) / 2.0;
    M(2,1) = u * (1.0 - gam);
    M(2,2) = gam - 1.0;
    return M;
}



// dW/dWp
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void dWdWp(
	const double gam,
    std::vector<double> &M,
    std::vector<double> const &W,
    int i)
{
    double rho, u;
    rho = W[i*3+0];
    u = W[i*3+1] / rho;
    M[0] = 1.0;
    M[1] = 0.0;
    M[2] = 0.0;
    M[3] = u;
    M[4] = rho;
    M[5] = 0.0;
    M[6] = u * u / 2.0;
    M[7] = u * rho;
    M[8] = 1.0 / (gam - 1.0);
}

std::vector <Eigen::Matrix3d> ddWpdWdWp(
    const double gam,
    const std::vector<double> &W,
    const int i)
{
    double rho, u;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;

    std::vector <Eigen::Matrix3d> M(3);
    M[0].setZero();
    M[1].setZero();
    M[2].setZero();

    M[1](0,0) = u / (rho * rho);
    M[1](0,1) = -1.0 / rho;
    M[1](1,0) = -1.0 / (rho * rho);

    M[2](0,1) = u * (gam - 1.0);
    M[2](1,1) = 1.0 - gam;

    return M;
}
