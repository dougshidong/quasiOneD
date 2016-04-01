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
#include"quasiOneD.h"
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"residuald1.h"
#include"parametrization.h"
#include"objectiveDerivatives.h"

using namespace Eigen;

VectorXd adjoint(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> &psi,
    std::vector <double> designVar)
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
    dRdW = evaldRdW(W, dx, dt, S);
//  matAFD2 = evaldRdW_FD(W, S, u[0]/c[0]);
    dRdWt = dRdW.transpose();
//  dRdWt = matAFD2.transpose();
    std::cout.precision(17);
//  std::cout<<dRdW<<std::endl;
//  std::cout<<matAFD2<<std::endl;
//  std::cout<<"(matAFD2 - dRdW).norm() / dRdW.norm():"<<std::endl;
//  std::cout<<(matAFD2 - dRdW).norm() / dRdW.norm()<<std::endl;

//  dRdWt.coeffRef(dRdWt.rows() - 3, dRdWt.cols() - 3) += 0.00001;
//  dRdWt.coeffRef(dRdWt.rows() - 2, dRdWt.cols() - 2) += 0.00001;
//  dRdWt.coeffRef(dRdWt.rows() - 1, dRdWt.cols() - 1) += 0.00001;

    // Build B matrix
    // Evaluate dIcdW
    VectorXd bvec(3 * nx);
    bvec = -evaldIcdW(W, dx);
//  std::cout<<"Vector B:"<<std::endl;
//  std::cout<<bvec<<std::endl;

    VectorXd psiV(3 * nx);
    psiV.setZero();
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    int eig_solv = 0;
    int directSolve = 1;
    if(directSolve == 1)
    {
        psiV = solveSparseAXB(dRdWt, bvec, eig_solv);
    }
    else
    {
        psiV = itSolve(dRdWt, bvec);
    }

    // If supersonic copy psi2 onto psi1 garbage
    if(u[0] > c[0])
    {
        psiV(0) = psiV(3);
        psiV(1) = psiV(4);
        psiV(2) = psiV(5);
    }
//  std::cout<<"Adjoint Result:"<<std::endl;
//  std::cout<<psiV<<std::endl;

    // Save Adjoint
    FILE *Results;
    Results = fopen("Adjoint.dat", "w");
    fprintf(Results, "%d\n", nx);
    for(int k = 0; k < 3; k++)
    for(int i = 0; i < nx; i++)
        fprintf(Results, "%.15f\n", psiV(i * 3 + k));

    fclose(Results);

    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS = evaldIcdS();

    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate psiV * dRdS
    VectorXd psidRdS(nx + 1);
    psidRdS = evalpsidRdS(psiV, Flux, p);

    // Finite Difference dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    MatrixXd dRdSFD(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    dRdSFD = evaldRdS_FD(Flux, S, W);
    std::cout<<"(dRdSFD - dRdS).norm() / dRdS.norm()"<<std::endl;
    std::cout<<(dRdSFD - dRdS).norm() / dRdS.norm()<<std::endl;

    VectorXd psidRdSFD(nx + 1);
    psidRdSFD.setZero();
    psidRdSFD = psiV.transpose() * dRdS;
    std::cout<<"(psidRdSFD - psidRdS).norm() / psidRdS.norm()"<<std::endl;
    std::cout<<(psidRdSFD - psidRdS).norm() / psidRdS.norm()<<std::endl;
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);

    // Evaluate dIdDes
    VectorXd grad(designVar.size());
    grad = psidRdSFD.transpose() * dSdDes;

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
    VectorXd psiV,
    std::vector <double> Flux,
    std::vector <double> p)
{
    VectorXd psidRdS(nx + 1);
    psidRdS.setZero();
    for(int i = 2; i < nx - 1; i++)
    for(int k = 0; k < 3; k++)
    {
        psidRdS(i) += psiV((i - 1) * 3 + k) * Flux[i * 3 + k];
        psidRdS(i) -= psiV(i * 3 + k) * Flux[i * 3 + k];
        if(k == 1)
        {
            psidRdS(i) -= psiV((i - 1) * 3 + k) * p[i - 1];
            psidRdS(i) += psiV(i * 3 + k) * p[i];
        }
    }

    // Evaluate psiV * dRdS neat the Boundaries
    for(int k = 0; k < 3; k++)
    {
        // Cell 0 Inlet is not a function of the shape

        // Cell 1
        psidRdS(1) -= psiV(1 * 3 + k) * Flux[1 * 3 + k];

        // Cell nx - 1
        psidRdS(nx - 1) += psiV((nx - 2) * 3 + k) * Flux[(nx - 1) * 3 + k];

        // Cell nx Outlet is not a function of the shape

        if(k == 1)
        {
            psidRdS(1) += psiV(1 * 3 + k) * p[1];
            psidRdS(nx - 1) -= psiV((nx - 2) * 3 + k) * p[nx - 1];
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

    JacobiSVD<MatrixXd> svd(matAdense);
    double svdmax = svd.singularValues()(0);
    double svdmin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = svdmax / svdmin;
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

VectorXd itSolve(SparseMatrix <double> A, VectorXd b)
{
    double resi1 = 1, resi2 = 1, resi3 = 1;
    // Directly Solve the Linear System Iteratively
    // Using Sub-Matrices
    //  --------------
    // |  A1   |  A2  |   | b1 |
    // |--------------| = |    |
    // |  A3   |  A4  |   | b2 |
    //  --------------
    MatrixXd A1(3 * (nx - 1), 3 * (nx - 1));
    MatrixXd A2(3 * (nx - 1), 3 * (nx - 1));
    MatrixXd A3(3 * (nx - 1), 3 * (nx - 1));
    MatrixXd A4(3 * (nx - 1), 3 * (nx - 1));
    A1 = MatrixXd(A.block(0, 0, 3 * (nx - 1), 3 * (nx - 1)));
    A2 = MatrixXd(A.block(3 * (nx - 2), 3 * (nx - 1), 3, 3));
    A3 = MatrixXd(A.block(3 * (nx - 1), 3 * (nx - 2), 3, 3));
    A4 = MatrixXd(A.block(3 * (nx - 1), 3 * (nx - 1), 3, 3));
    A4 = A4 + MatrixXd(3, 3).setIdentity() * 0.000001;
//  std::cout<<A<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;


    VectorXd b1(3 * (nx - 1)), b2(3);
    b1 = b.head(3 * (nx - 1));
    b2 = b.tail(3);

    VectorXd b1mod(3 * (nx - 1)), b2mod(3);
    b1mod = b1;
    b2mod = b2;

    VectorXd fullX(3 * nx);
    fullX.setZero();
    fullX.setOnes();
    fullX = MatrixXd(A).fullPivLu().solve(b);

    VectorXd x1(3 * (nx - 1));
    VectorXd x2(3);
    x1 = fullX.head(3 * (nx - 1));
    x2 = fullX.tail(3);

//  b1mod.tail(3) = b1.tail(3) - A2 * x2.tail(3);
//  b2mod.tail(3) = b2.tail(3) - A3 * x1.tail(3);

    double tol1 = 5e-13;
    double tol2 = tol1;
    double tol3 = 1;//tol1;
    int it = 0;
    while(resi1 > tol1 || resi2 > tol2 || resi3 > tol3)
    {
        it++;
        x1 = A1.fullPivLu().solve(b1mod);
        x2 = A4.fullPivLu().solve(b2mod);
        b1mod.tail(3) = b1.tail(3) - A2 * x2.tail(3);
        b2mod.tail(3) = b2.tail(3) - A3 * x1.tail(3);

        resi1 = (A1 * x1 - b1mod).norm();
        resi2 = (A4 * x2 - b2mod).norm();
        fullX.head(3 * (nx - 1)) = x1;
        fullX.tail(3) = x2;
        resi3 = (A * fullX - b).norm();

        std::cout<<"Iteration: "<<it
                 <<" resi1: "<<resi1
                 <<" resi2: "<<resi2
                 <<" resi3: "<<resi3
                 <<std::endl;
    }
    JacobiSVD<MatrixXd> svd(A1);
    double svdmax = svd.singularValues()(0);
    double svdmin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = svdmax / svdmin;
    std::cout<<"Condition Number A1"<<std::endl;
    std::cout<<cond<<std::endl;
    std::cout<<"Max/Min Singular Values"<<std::endl;
    std::cout<<svdmax<< " / "<<svdmin<<std::endl;
    JacobiSVD<MatrixXd> svd2(A4);
    svdmax = svd2.singularValues()(0);
    svdmin = svd2.singularValues()(svd2.singularValues().size()-1);
    cond = svdmax / svdmin;
    std::cout<<"Condition Number A4"<<std::endl;
    std::cout<<cond<<std::endl;
    std::cout<<"Max/Min Singular Values"<<std::endl;
    std::cout<<svdmax<< " / "<<svdmin<<std::endl;

    std::cout<<"||Ax - b||"<<std::endl;
    std::cout<<(A * fullX - b).norm()<<std::endl;

    return fullX;
}
