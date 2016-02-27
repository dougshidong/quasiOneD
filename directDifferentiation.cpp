#include<math.h>
#include<vector>
#include<iostream>
#include<iomanip>
#include<Eigen/Eigen>
#include "globals.h"
#include "adjoint.h"
#include "convert.h"
#include "flux.h"

using namespace Eigen;
VectorXd directDifferentiation(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar)
{
    // Direct Differentiation
    // I = Ic(W, S)
    // R = R(W, S) = 0 @ SS
    // W = W(W0, S)
    // DI   dIc   dIc DW
    // -- = --- + --- --
    // DS   dS    dW  DS
    //
    // DR   dR   dR DW           DW   -( dR ) ^ (-1) ( dR )
    // -- = -- + -- -- = 0  -->  -- =  ( -- )        ( -- )
    // DS   dS   dW DS           DS    ( dW )        ( dS )

    //Get Primitive Variables
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T); 

    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();

    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);

    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS_FD(Flux, S, W);

    // Evaluate dQdW
    std::vector <double> dQdW(3 * nx, 0);
    evaldQdW(dQdW, W, S);

    // Get Jacobians
    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);
    if(FluxScheme == 1) ScalarJac(W, Ap_list, An_list);

    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(Ap_list, An_list, W, dQdW, dx, dt, S, u[0]/c[0]);

    // Solve DWDS
    MatrixXd DWDS(3 * nx, nx + 1);
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    int eig_solv = 0;
    DWDS = solveSparseAXB(-dRdW, dRdS, eig_solv);

    // Evaluate DIDS
    VectorXd DIDS(nx + 1);
    DIDS = (dIcdS.transpose() + dIcdW.transpose() * DWDS);

    // Evaluate dSdDesign
    MatrixXd dSdDesign(nx + 1, designVar.size());
    dSdDesign = evaldSdDesign(x, dx, designVar);
    
    // Evaluate dIdDesign
    VectorXd dIdDesign(designVar.size());
    dIdDesign = DIDS.transpose() * dSdDesign;

    VectorXd grad(designVar.size());
    grad = dIdDesign;
    std::cout<<"Gradient from Direct Differentiation:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return grad;
}
