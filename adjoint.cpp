// Calculates the discrete costate fluxes
#include<iostream>
#include<math.h>
#include<vector>
#include<Eigen/Core>
#include<Eigen/LU>
#include<Eigen/Sparse>
#include<Eigen/SVD>
#include<stdio.h>
#include<iomanip>
#include"adjoint.h"
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"residuald1.h"
#include"parametrization.h"
#include"objectiveDerivatives.h"
#include"output.h"
#include"petscGMRES.h"

using namespace Eigen;

VectorXd adjoint(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> area,
    std::vector <double> W,
    std::vector <double> designVar,
    VectorXd &psi)
{
    //Get Primitive Variables
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    // Evalutate dt and d(dt)dW
    std::vector <double> dt(nx, 1);

    // Build A matrix
    SparseMatrix <double> dRdWt, dRdW;
    SparseMatrix <double> matAFD, matAFD2;
    dRdW = evaldRdW(W, dx, dt, area);
//  matAFD2 = evaldRdW_FD(W, area, u[0]/c[0]);
    dRdWt = dRdW.transpose();
//  dRdWt = matAFD2.transpose();
    std::cout.precision(17);

    // Build B matrix
    // Evaluate dIcdW
    VectorXd bvec(3 * nx);
    bvec = evaldIcdW(W, dx);
//  std::cout<<"Vector B:"<<std::endl;
//  std::cout<<bvec<<std::endl;

    psi.setZero();
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    //int eig_solv = 0;
    //int directSolve = 1;
    //if(directSolve == 1)
    //{
    //    psi = solveSparseAXB(-dRdWt, bvec, eig_solv);
    //}
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
    slusolver.compute(-dRdWt);
    if(slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;
    psi = slusolver.solve(bvec);
//  psi = solveGMRES(-dRdWt, bvec);

    // If supersonic copy psi2 onto psi1 garbage
    if(u[0] > c[0])
    {
        psi(0) = psi(3);
        psi(1) = psi(4);
        psi(2) = psi(5);
    }
//  std::cout<<"Adjoint Result:"<<std::endl;
//  std::cout<<psi<<std::endl;

    // Save Adjoint
    outVec("adjoint.dat", "w", psi);

    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS = evaldIcdS();

    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate psi * dRdS
    VectorXd psidRdS(nx + 1);
    psidRdS = evalpsidRdS(psi, Flux, p);

    // Finite Difference dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    MatrixXd dRdSFD(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, area, W);
    dRdSFD = evaldRdS_FD(Flux, area, W);

    VectorXd psidRdSFD(nx + 1);
    psidRdSFD.setZero();
    psidRdSFD = psi.transpose() * dRdS;

    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);

    // Evaluate dIdDes
    VectorXd grad(designVar.size());
    grad = psidRdS.transpose() * dSdDes;

    std::cout<<"Gradient from Adjoint:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return grad;
}

VectorXd buildbMatrix(std::vector <double> dIcdW)
{
    VectorXd matb(3 * nx);

    for(int i = 0; i < nx; i++)
    for(int k = 0; k < 3; k++)
        matb(i * 3 + k) = -dIcdW[i * 3 + k];

    return matb;
}

VectorXd evalpsidRdS(
    VectorXd psi,
    std::vector <double> Flux,
    std::vector <double> p)
{
    VectorXd psidRdS(nx + 1);
    psidRdS.setZero();
    for(int i = 2; i < nx - 1; i++)
    for(int k = 0; k < 3; k++)
    {
        psidRdS(i) += psi((i - 1) * 3 + k) * Flux[i * 3 + k];
        psidRdS(i) -= psi(i * 3 + k) * Flux[i * 3 + k];
        if(k == 1)
        {
            psidRdS(i) -= psi((i - 1) * 3 + k) * p[i - 1];
            psidRdS(i) += psi(i * 3 + k) * p[i];
        }
    }

    // Evaluate psi * dRdS neat the Boundaries
    for(int k = 0; k < 3; k++)
    {
        // Cell 0 Inlet is not a function of the shape

        // Cell 1
        psidRdS(1) -= psi(1 * 3 + k) * Flux[1 * 3 + k];

        // Cell nx - 1
        psidRdS(nx - 1) += psi((nx - 2) * 3 + k) * Flux[(nx - 1) * 3 + k];

        // Cell nx Outlet is not a function of the shape

        if(k == 1)
        {
            psidRdS(1) += psi(1 * 3 + k) * p[1];
            psidRdS(nx - 1) -= psi((nx - 2) * 3 + k) * p[nx - 1];
        }
    }
    return psidRdS;
}


MatrixXd solveSparseAXB(
    SparseMatrix <double> A,
    MatrixXd B, int eig_solv)
{
    MatrixXd X(A.rows(), B.cols());
    X.setZero();
    MatrixXd matAdense(3 * nx, 3 * nx);
    MatrixXd eye(3 * nx, 3 * nx);
    eye.setIdentity();
    matAdense = A * eye;

    double offset = 0;//0.00001;
    matAdense = matAdense + eye * offset;

//  JacobiSVD<MatrixXd> svd(matAdense);
//  double svdmax = svd.singularValues()(0);
//  double svdmin = svd.singularValues()(svd.singularValues().size()-1);
//  double cond = svdmax / svdmin;
//  std::cout<<"Condition Number SVD"<<std::endl;
//  std::cout<<cond<<std::endl;
//  std::cout<<"Max/Min Singular Values"<<std::endl;
//  std::cout<<svdmax<< " / "<<svdmin<<std::endl;

    // Sparse LU
    if(eig_solv == 0)
    {
        SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
        slusolver.analyzePattern(A);
        slusolver.factorize(A);

        if(slusolver.info() != 0)
            std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

        // Solve for X
        X = slusolver.solve(B);
    }
    // Dense LU full pivoting
    if(eig_solv == 1)
    {
        // Full Pivoting LU Factorization
        X = matAdense.fullPivLu().solve(B);
    }
    // Iterative LU
    if(eig_solv == 2)
    {
        BiCGSTAB<SparseMatrix <double> > itsolver;
        itsolver.compute(A);
        if(itsolver.info() == 0)
            std::cout<<"Iterative Factorization success"<<std::endl;
        else
            std::cout<<"Factorization failed. Error: "<<itsolver.info()<<std::endl;
        std::cout << "#iterations:     " << itsolver.iterations() << std::endl;
        std::cout << "estimated error: " << itsolver.error()      << std::endl;
        X = itsolver.solve(B);
    }
//  std::cout<<"||Ax - B||"<<std::endl;
//  std::cout<<(matAdense * X - B).norm()<<std::endl;

    return X;
}

