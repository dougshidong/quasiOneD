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
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> W,
    std::vector<double> designVar,
    VectorXd &psi)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS = evaldIcdS();
    // Evaluate dIcdDes
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;

    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    //// Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);

    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    //VectorXd psi(3 * n_elem);
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if (slusolver1.info() != 0)
		std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    psi = slusolver1.solve(dIcdW);

    VectorXd dIdDes = dIcdDes.transpose() + psi.transpose()*dRdDes;

	return dIdDes;
}

MatrixXd solveSparseAXB(
    SparseMatrix<double> A,
    MatrixXd B, int eig_solv)
{
    MatrixXd X(A.rows(), B.cols());
    X.setZero();

//  JacobiSVD<MatrixXd> svd(matAdense);
//  double svdmax = svd.singularValues()(0);
//  double svdmin = svd.singularValues()(svd.singularValues().size()-1);
//  double cond = svdmax / svdmin;
//  std::cout<<"Condition Number SVD"<<std::endl;
//  std::cout<<cond<<std::endl;
//  std::cout<<"Max/Min Singular Values"<<std::endl;
//  std::cout<<svdmax<< " / "<<svdmin<<std::endl;

    // Sparse LU
    if (eig_solv == 0)
    {
        SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver;
        slusolver.analyzePattern(A);
        slusolver.factorize(A);

        if (slusolver.info() != 0)
            std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

        // Solve for X
        X = slusolver.solve(B);
    }
    // Dense LU full pivoting
    if (eig_solv == 1)
    {
		MatrixXd matAdense(3 * n_elem, 3 * n_elem);
		MatrixXd eye(3 * n_elem, 3 * n_elem);
		eye.setIdentity();
		matAdense = A * eye;

		double offset = 0;//0.00001;
		matAdense = matAdense + eye * offset;
        // Full Pivoting LU Factorization
        X = matAdense.fullPivLu().solve(B);
    }
    // Iterative LU
    if (eig_solv == 2)
    {
        BiCGSTAB<SparseMatrix<double> > itsolver;
        itsolver.compute(A);
        if (itsolver.info() == 0)
            std::cout<<"Iterative Factorization success"<<std::endl;
        else
            std::cout<<"Factorization failed. Error: "<<itsolver.info()<<std::endl;
        std::cout << "#iterations:     " << itsolver.iterations() << std::endl;
        std::cout << "estimated error: " << itsolver.error()      << std::endl;
        X = itsolver.solve(B);
    }

    return X;
}

