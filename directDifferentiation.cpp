#include<math.h>
#include<vector>
#include<iostream>
#include<iomanip>
#include<Eigen/Eigen>
#include "globals.h"
#include "adjoint.h"
#include "convert.h"
#include "flux.h"
#include "residuald1.h"
#include "parametrization.h"
#include "objectiveDerivatives.h"


using namespace Eigen;
MatrixXd evaldWdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> W,
    std::vector<double> designVar)
{
    // DR   dR   dR DW           DW   -( dR ) ^ (-1) ( dR )
    // -- = -- + -- -- = 0  -->  -- =  ( -- )        ( -- )
    // DS   dS   dW DS           DS    ( dW )        ( dS )

    //Get Primitive Variables
    std::vector<double> rho(n_elem), u(n_elem), e(n_elem);
    std::vector<double> T(n_elem), p(n_elem), c(n_elem), Mach(n_elem);
    WtoP(W, rho, u, e, p, c, T);

    // Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);

    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);

    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);

    //Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;

    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);

    // Solve DWDS
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    int eig_solv = 0;
    dWdDes = solveSparseAXB(-dRdW, dRdDes, eig_solv);

    return dWdDes;
}

VectorXd directDifferentiation(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> W,
    std::vector<double> designVar)
{
    // Direct Differentiation
    // I = Ic(W, area)
    // R = R(W, area) = 0 @ SS
    // W = W(W0, area)
    // DI   dIc   dIc DW
    // -- = --- + --- --
    // DS   dS    dW  DS
    //
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS.setZero();
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;

    // Evaluate dWdDes
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    dWdDes = evaldWdDes(x, dx, area, W, designVar);
    // Evaluate dIcdDes
    VectorXd dIdDes(nDesVar);
    dIdDes = (dIcdDes.transpose() + dIcdW.transpose() * dWdDes);


    VectorXd grad(designVar.size());
    grad = dIdDes;
    std::cout<<"Gradient from Direct Differentiation:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return grad;
}
