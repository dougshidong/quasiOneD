// Convert W to F or primitive variables

#include<vector>
#include<Eigen/Core>
#include<math.h>
#include"convert.h"

// Get pressure vector
void get_all_p(
	const double gam,
	std::vector<double> const &W,
	std::vector<double> &p)
{
	int n_elem = p.size();
    for (int i = 0; i < n_elem; i++) {
		p[i] = get_p(gam, W[i*3+0], W[i*3+1], W[i*3+2]);
    }
}


void WtoF(
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
        F[i*3+1] = w2 * w2 / w1 + (gam - 1) * ( w3 - w2 * w2 / (2 * w1) );
        F[i*3+2] = ( w3 + (gam - 1) * (w3 - w2 * w2 / (2 * w1)) ) * w2 / w1;
    }
}

// get Q
void WtoQ(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &Q,
    std::vector<double> const &area)
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

// Get primitive variables
void WtoP(
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
        rho[i] = W[i*3+0];
        u[i] = W[i*3+1] / rho[i];
        e[i] = W[i*3+2];
    }
}

// Get other primitive variables
void WtoP2(
	const double gam,
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &p)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
        rho[i] = W[i*3+0];
        u[i] = W[i*3+1] / rho[i];
        p[i] = (gam - 1.0) * ( W[i*3+2] - (pow(W[i*3+1], 2.0) / W[i*3+0]) / 2.0 );
    }
}

// Given rho, u p, get W
void PtoW(
	const double gam,
    std::vector<double> &W,
    std::vector<double> const &Wp)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
        W[i*3+0] = Wp[i*3+0];
        W[i*3+1] = Wp[i*3+0] * Wp[i*3+1];
        W[i*3+2] = Wp[i*3+2] / (gam - 1.0) + 
                       (pow(Wp[i*3+1], 2.0) * Wp[i*3+0]) / 2.0;
    }
}

// Get more primitive variables
void WtoP(
	const double gam,
	const double R,
    std::vector<double> const &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e,
    std::vector<double> &p,
    std::vector<double> &c,
    std::vector<double> &T)
{
	int n_elem = W.size()/3;
    for (int i = 0; i < n_elem; i++)
    {
        rho[i] = W[i*3+0];
        u[i] = W[i*3+1] / rho[i];
        e[i] = W[i*3+2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
        T[i] = p[i] / (rho[i] * R);
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
    std::vector<double> const &W,
    int i)
{
    Eigen::MatrixXd M(3, 3);
    double rho, u;
    rho = W[i*3+0];
    u = W[i*3+1] / rho;
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
