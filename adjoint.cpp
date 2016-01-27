// Calculates the discrete costate fluxes
#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include<vector>
#include<Eigen/Eigen>
//#include<Eigen/Dense>
//#include<Eigen/Sparse>
//#include <Eigen/IterativeLinearSolvers>
//#include<Eigen/SparseLU>
//#include<Eigen/LU>
#include <stdio.h>
#include <iomanip>
#include "globals.h"
#include "flovar.h"
#include "convert.h"

using namespace Eigen;

SparseMatrix<double> buildAMatrix(std::vector <double> Ap,
                                  std::vector <double> An,
                                  std::vector <double> dBidWi,
                                  std::vector <double> dBidWd,
                                  std::vector <double> dBodWd,
                                  std::vector <double> dBodWo,
                                  std::vector <double> dQdW,
                                  std::vector <double> dx,
                                  std::vector <double> dt,
                                  std::vector <double> S,
                                  double Min);

VectorXd buildbMatrix(std::vector <double> dIcdW);

void StegerJac(std::vector <double> W,
               std::vector <double> &Ap_list,
               std::vector <double> &An_list,
               std::vector <double> &Flux);

void BCJac(std::vector <double> W,
           std::vector <double> dt,
           std::vector <double> dx,
           std::vector <double> ddtdW,
           std::vector <double> &dBidWi,
           std::vector <double> &dBidWd,
           std::vector <double> &dBodWd,
           std::vector <double> &dBodWo,
           std::vector <double> &B1,
           std::vector <double> &BN);

void evaldIcdW(std::vector <double> &dIcdW,
               std::vector <double> W,
               std::vector <double> S);

void evaldQdW(std::vector <double> &dQdW,
                   std::vector <double> W,
                   std::vector <double> S);

void evalddtdW(std::vector <double> &ddtdW,
                 std::vector <double> rho,
                 std::vector <double> u,
                 std::vector <double> p,
                 std::vector <double> c,
                 std::vector <double> dx,
                 std::vector <double> S);

std::vector <double> adjoint(std::vector <double> x, 
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
    std::vector <double> dt(nx, 1), ddtdW(nx, 0);
//  for(int i = 0; i < nx; i++)
//      dt[i] = CFL * dx[i] / ( fabs(u[i] + c[i]) );
//  evalddtepdW(ddtdW, rho, u, p, c, S);

    // Evaluate dQdW
    std::vector <double> dQdW(3 * nx, 0);
    evaldQdW(dQdW, W, S);

    // Get Jacobians and Fluxes
    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);
    std::vector <double> Flux(3 * (nx + 1), 0);
    StegerJac(W, Ap_list, An_list, Flux);
    
    // Transposed Boundary Flux Jacobians
    std::vector <double> dBidWi(3 * 3, 0);
    std::vector <double> dBidWd(3 * 3, 0);
    std::vector <double> dBodWd(3 * 3, 0);
    std::vector <double> dBodWo(3 * 3, 0);
    std::vector <double> B1(3, 0);
    std::vector <double> BN(3, 0);
    BCJac(W, dt, dx, ddtdW, dBidWi, dBidWd, dBodWd, dBodWo, B1, BN);

    // Build A matrix
    SparseMatrix <double> matA;
    matA = buildAMatrix(Ap_list, An_list, dBidWi, dBidWd,
                        dBodWd, dBodWo, dQdW, dx, dt, S, u[0]/c[0]);
    // Transpose
    SparseMatrix <double> matAtranspose;
    matAtranspose = matA.transpose();
//    matAtranspose.makeCompressed();
    std::cout.precision(17);
    std::cout<<matAtranspose<<std::endl;

    // Evaluate dIcdW
    std::vector <double> dIcdW(3 * nx, 0);
    evaldIcdW(dIcdW, W, dx);

    // Build B matrix
    VectorXd bvec(3 * nx);
    bvec.setZero();

    bvec = buildbMatrix(dIcdW);
    std::cout<<"Vector B"<<std::endl;
    std::cout<<bvec<<std::endl;

    double eig_solv = 0;
    VectorXd xvec(3 * nx);
    xvec.setZero();
    MatrixXd matAdense(3 * nx, 3 * nx);
    MatrixXd eye(3 * nx, 3 * nx);
    eye.setIdentity();
    matAdense = matAtranspose * eye;
    if(eig_solv == 0)
    {
    // Setup Solver and Factorize A
        SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
        slusolver.analyzePattern(matAtranspose);
        slusolver.factorize(matAtranspose);
        
        if(slusolver.info() == 0)
            std::cout<<"Factorization success"<<std::endl;
        else
            std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;
    
        // Solve for X
        xvec = slusolver.solve(bvec);
        std::cout<<"Adjoint Result:"<<std::endl;
        std::cout<<xvec<<std::endl;
    }
    if(eig_solv == 1)
    {
        // Full Pivoting LU Factorization
        xvec = matAdense.fullPivLu().solve(bvec);
        std::cout<<"Adjoint Result:"<<std::endl;
        std::cout<<xvec<<std::endl;
        std::cout<<"Identity"<<std::endl;
        std::cout<<(matAdense.inverse() * matAdense).norm()<<std::endl;
    }
    if(eig_solv == 2)
    {
        BiCGSTAB<SparseMatrix <double> > itsolver;
        itsolver.compute(matAtranspose);
        if(itsolver.info() == 0)
            std::cout<<"Iterative Factorization success"<<std::endl;
        else
            std::cout<<"Factorization failed. Error: "<<itsolver.info()<<std::endl;
        std::cout << "#iterations:     " << itsolver.iterations() << std::endl;
        std::cout << "estimated error: " << itsolver.error()      << std::endl;
        xvec = itsolver.solve(bvec);
        std::cout<<"Adjoint Result:"<<std::endl;
        std::cout<<xvec<<std::endl;
    }
    std::cout<<"Condition Number"<<std::endl;
    std::cout<<matAdense.inverse().norm() * matAdense.norm()<<std::endl;
    std::cout<<"Condition Number SVD"<<std::endl;
    JacobiSVD<MatrixXd> svd(matAdense);
    double cond = svd.singularValues()(0) 
    / svd.singularValues()(svd.singularValues().size()-1);
    std::cout<<cond<<std::endl;
    std::cout<<"||Ax - b||"<<std::endl;
    std::cout<<(matAdense * xvec - bvec).norm()<<std::endl;


    if(u[0] > c[0])
    {
        xvec(0) = xvec(3);
        xvec(1) = xvec(4);
        xvec(2) = xvec(5);
    }

    // Print out Adjoint 
    FILE *Results;
    Results = fopen("Adjoint.dat", "w");
    fprintf(Results, "%d\n", nx);
    for(int k = 0; k < 3; k++)
    for(int i = 0; i < nx; i++)
        fprintf(Results, "%.15f\n", xvec(i * 3 + k));

    fclose(Results);

    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();

    // Evaluate psi * dRdS
    VectorXd psidRdS(nx + 1);
    psidRdS.setZero();
    for(int i = 2; i < nx - 1;  i++)
    for(int k = 0; k < 3; k++)
    {
        psidRdS(i) += xvec((i - 1) * 3 + k) * Flux[i * 3 + k];
        psidRdS(i) -= xvec(i * 3 + k) * Flux[i * 3 + k];
        if(k == 1)
        {
            psidRdS(i) -= xvec((i - 1) * 3 + k) * p[i - 1];
            psidRdS(i) += xvec(i * 3 + k) * p[i];
        }
    }


    for(int k = 0; k < 3; k++)
    {
        // 0 Inlet is not a function of the shape
//      psidRdS(0) -= xvec(0 * 3 + k) * B1[k] * dx[0] * S[0] / 2.0;
//      psidRdS(0) -= xvec(0 * 3 + k) * B1[k];

        // 1
        psidRdS(1) -= xvec(1 * 3 + k) * Flux[1 * 3 + k];
//        psidRdS(1) -= xvec(1 * 3 + k) * B1[k] * dx[0] * S[1] / 2.0;

        // nx - 1
//        psidRdS(nx - 1) -= xvec((nx - 1) * 3 + k) * BN[k] * dx[nx - 1] * S[nx - 1] / 2.0;
        psidRdS(nx - 1) += xvec((nx - 2) * 3 + k) * Flux[(nx - 1) * 3 + k];

        // nx Outlet is not a function of the shape
//      psidRdS(nx) -= xvec((nx - 1) * 3 + k) * BN[k] * dx[nx - 1] * S[nx] / 2.0;
//      psidRdS(nx) -= xvec((nx - 1) * 3 + k) * BN[k];

        if(k == 1)
        {
            psidRdS(nx - 1) -= xvec((nx - 2) * 3 + k) * p[nx - 1];
            psidRdS(1) += xvec(1 * 3 + k) * p[1];
        }
    }

    std::cout<<"psidRdS:"<<std::endl;
    for(int i = 0; i < nx + 1; i++)
        std::cout<<psidRdS(i)<<std::endl;
    std::cout<<"psidRdSd"<<std::endl;
    
    // Evaluate dSdDesign
    MatrixXd dSdDesign(nx + 1, 3);
    double d1 = designVar[0], d2 = designVar[1], d3 = designVar[2];
    std::cout<<d1<<std::endl;
    std::cout<<d2<<std::endl;
    std::cout<<d3<<std::endl;
    double xh;
    for(int i = 0; i < nx + 1; i++)
    {
        xh = fabs(x[i] - dx[i] / 2.0);

        if(i == 0 || i == nx)
        {
            dSdDesign(i, 0) = 0;
            dSdDesign(i, 1) = 0;
            dSdDesign(i, 2) = 0;
        }
        else
        {
            dSdDesign(i, 0) = - pow(sin(PI * pow(xh, d2)), d3);
            dSdDesign(i, 1) = - d1 * d3 * PI * pow(xh, d2)
                              * cos(PI * pow(xh, d2)) * log(xh)
                              * pow(sin(PI * pow(xh, d2)), d3 - 1);
            dSdDesign(i, 2) = - d1 * log(sin(PI * pow(xh, d2)))
                              * pow(sin(PI * pow(xh, d2)), d3);
        }
    }

    VectorXd grad(designVar.size());
    grad = psidRdS.transpose() * dSdDesign;

    std::vector <double> gradient(designVar.size());
    for(int iDes = 0; iDes < 3; iDes++)
    {
       gradient[iDes] = grad(iDes);
    }
    std::cout<<"Gradient from Adjoint:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return gradient;
}


// Calculates Jacobian
// Steger-Warming Flux Splitting
void StegerJac(std::vector <double> W,
               std::vector <double> &Ap_list,
               std::vector <double> &An_list,
               std::vector <double> &Flux)
{
    double eps = 0.1;
    double gam = 1.4;
    double M[3][3] = {{0}},
           Minv[3][3] = {{0}},
           N[3][3] = {{0}},
           Ninv[3][3] = {{0}},
           lambdaP[3][3],
           lambdaN[3][3];
    double lambdaa[3];
    
    
    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
    
    std::vector <double> rho(nx), u(nx), p(nx), c(nx);

    double beta = gam - 1;

    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[0 * nx + i];
        u[i] = W[1 * nx + i] / rho[i];
        p[i] = (gam-1) * (W[2 * nx + i] - rho[i] * pow(u[i], 2) / 2);
        c[i] = sqrt( gam * p[i] / rho[i] );
    }


    for(int i = 0; i < nx; i++)
    {
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            Ap[row][col] = 0;
            An[row][col] = 0;
            tempP[row][col] = 0;
            tempN[row][col] = 0;
            prefix[row][col] = 0;
            suffix[row][col] = 0;
            lambdaP[row][col] = 0;
            lambdaN[row][col] = 0;
        }
    
        M[0][0] = 1.0;
        M[1][0] = -u[i] / rho[i];
        M[2][0] = 0.5 * u[i] * u[i] * beta;
        M[1][1] = 1.0 / rho[i];
        M[2][1] = -u[i] * beta;
        M[2][2] = beta;
        Minv[0][0] = 1.0;
        Minv[1][0] = u[i];
        Minv[2][0] = 0.5 * u[i] * u[i];
        Minv[1][1] = rho[i];
        Minv[2][1] = u[i] * rho[i];
        Minv[2][2] = 1.0 / beta;
        N[0][0] = 1.0;
        N[1][1] = rho[i] * c[i];
        N[2][1] = -rho[i] * c[i];
        N[0][2] = -1.0 / (c[i] * c[i]);
        N[1][2] = 1.0;
        N[2][2] = 1.0;
        Ninv[0][0] = 1.0;
        Ninv[0][1] = 1.0 / (2.0 * c[i] * c[i]);
        Ninv[0][2] = 1.0 / (2.0 * c[i] * c[i]);
        Ninv[1][1] = 1.0 / (2.0 * rho[i] * c[i]);
        Ninv[1][2] = -1.0 / (2.0 * rho[i] * c[i]);
        Ninv[2][1] = 0.5;
        Ninv[2][2] = 0.5;
        lambdaa[0] = u[i];
        lambdaa[1] = u[i] + c[i];
        lambdaa[2] = u[i] - c[i];
        
        for(int k = 0; k < 3; k++)
            if(lambdaa[k] > 0)
                lambdaP[k][k] = (lambdaa[k] + sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;
            else
                lambdaN[k][k] = (lambdaa[k] - sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;

        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                prefix[row][col]+= Minv[row][k] * Ninv[k][col];
                suffix[row][col]+= N[row][k] * M[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                tempP[row][col] += prefix[row][k] * lambdaP[k][col];
                tempN[row][col] += prefix[row][k] * lambdaN[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                Ap[row][col]+= tempP[row][k] * suffix[k][col];
                An[row][col]+= tempN[row][k] * suffix[k][col];
            }
        // could remove above loop and just use aplist and anlist
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int vec_pos = (i * 3 * 3) + (row * 3) + col;
            Ap_list[vec_pos] = Ap[row][col];
            An_list[vec_pos] = An[row][col];
        }

    }

    for(int i = 1; i < nx; i++)
    {
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
            int An_pos = (i * 3 * 3) + (row * 3) + col;
            Flux[i * 3 + row] += Ap_list[Ap_pos] * W[col * nx + (i - 1)]
                                 + An_list[An_pos] * W[col * nx + i];
        }
    }
    
    // Transpose the Jacobians
//  for(int i = 0; i < nx; i++)
//  for(int row = 0; row < 3; row++)
//  for(int col = 0; col < 3; col++)
//  {
//      // TRANSPOSED JACOBIANS
//      int vec_pos = (i * 3 * 3) + (row * 3) + col;
//      Ap_list[vec_pos] = Ap[col][row];
//      An_list[vec_pos] = An[col][row];
//  }

}

void evaldIcdW(std::vector <double> &dIcdW,
               std::vector <double> W,
               std::vector <double> dx)
{
    std::vector <double> ptarget(nx, 0);
    double dpdw[3], rho, u, p;
    ioTargetPressure(-1, ptarget);
    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;
        p = (gam - 1) * ( W[2 * nx + i] - rho * u * u / 2.0 );

        dpdw[0] = (gam - 1) / 2.0 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dIcdW[i * 3 + 0] = (p / ptin - ptarget[i]) * dpdw[0] * dx[i] / ptin;
        dIcdW[i * 3 + 1] = (p / ptin - ptarget[i]) * dpdw[1] * dx[i] / ptin;
        dIcdW[i * 3 + 2] = (p / ptin - ptarget[i]) * dpdw[2] * dx[i] / ptin;
    }
}

void evaldQdW(std::vector <double> &dQdW,
                   std::vector <double> W,
                   std::vector <double> S)
{
    double dpdw[3], rho, u, dS;
    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;

        dpdw[0] = (gam - 1) / 2.0 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dS = S[i + 1] - S[i];

        dQdW[i * 3 + 0] = dpdw[0] * dS;
        dQdW[i * 3 + 1] = dpdw[1] * dS;
        dQdW[i * 3 + 2] = dpdw[2] * dS;
    }
}

void BCJac(std::vector <double> W,
           std::vector <double> dt,
           std::vector <double> dx,
           std::vector <double> ddtdW,
           std::vector <double> &dBidWi,
           std::vector <double> &dBidWd,
           std::vector <double> &dBodWd,
           std::vector <double> &dBodWo,
           std::vector <double> &B1,
           std::vector <double> &BN)
{
    std::vector <double> rho(nx), u(nx), e(nx), p(nx), c(nx), T(nx);
    std::vector <double> dbdwp(9, 0), dwpdw(9);

    for(int i = 0; i < 9; i++)
    {
        dBidWi[i] = 0;
        dBidWd[i] = 0;
        dBodWd[i] = 0;
        dBodWo[i] = 0;
    }

    WtoP(W, rho, u, e, p, c, T);

    // ************************
    // OUTLET JACOBIANS
    // ************************

    double i1, i2;
    double r1, r2, p1, p2, u1, u2, c1, c2, t1;
    i1 = nx - 1;
    i2 = nx - 2;
    r1 = rho[i1];
    r2 = rho[i2];
    p1 = p[i1];
    p2 = p[i2];
    u1 = u[i1];
    u2 = u[i2];
    c1 = c[i1];
    c2 = c[i2];
    t1 = T[i1];

    // Shorthand
    double gamr, fu, drho, dp, du, cr, uu;
    drho = r1 - r2;
    dp = p1 - p2;
    du = u1 - u2;
    cr = r1 * c1;
    uu = u1 * u1;

    // Speed of Sound
    double dc1dr1, dc2dr2, dc1dp1, dc2dp2;
    dc1dr1 = - p1 * gam / (2.0 * cr * r1);
    dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
    dc1dp1 = gam / (2.0 * cr);
    dc2dp2 = gam / (2.0 * c2 * r2);

    double eig1, eig2, eig3;
    double deig1du1, deig1du2;
    double deig2dr1, deig2du1, deig2dp1, deig2dr2, deig2du2, deig2dp2;
    double deig3dr1, deig3du1, deig3dp1, deig3dr2, deig3du2, deig3dp2;
    // Eigenvalue
    eig1 = (u1 + u2) / 2.0;
    eig2 = eig1 + (c1 + c2) / 2.0;
    eig3 = eig1 - (c1 + c2) / 2.0;

    deig1du1 = 1.0 / 2.0;
    deig1du2 = 1.0 / 2.0;

    deig2dr1 = dc1dr1 / 2.0;
    deig2du1 = deig1du1;
    deig2dp1 = dc1dp1 / 2.0;
    deig2dr2 = dc2dr2 / 2.0;
    deig2du2 = deig1du2;
    deig2dp2 = dc2dp2 / 2.0;

    deig3dr1 = - dc1dr1 / 2.0;
    deig3du1 = deig1du1;
    deig3dp1 = - dc1dp1 / 2.0;
    deig3dr2 = - dc2dr2 / 2.0;
    deig3du2 = deig1du2;
    deig3dp2 = - dc2dp2 / 2.0;

    // Riemann invariants
    double R1, R2, R3;
    double dR1dr1, dR1du1, dR1dp1, dR1dr2, dR1du2, dR1dp2;
    double dR2dr1, dR2du1, dR2dp1, dR2dr2, dR2du2, dR2dp2;
    double dR3dr1, dR3du1, dR3dp1, dR3dr2, dR3du2, dR3dp2;
    R1 = - eig1 * (drho - dp / (c1 * c1));
    R2 = - eig2 * (dp + cr * du);
    R3 = - eig3 * (dp - cr * du);

    dR1dr1 = - eig1 * (1.0 + 2.0 * dp * dc1dr1 / pow(c1, 3) );
    dR1du1 = deig1du1 * (dp - c1 * c1 * drho) / (c1 * c1);
    dR1dp1 = eig1 * (c1 - 2 * dp * dc1dp1) / pow(c1, 3);
    dR1dr2 = eig1;
    dR1du2 = deig1du2 * (dp - c1 * c1 * drho) / (c1 * c1);
    dR1dp2 = - eig1 / (c1 * c1);

    dR2dr1 = - du * eig2 * (c1 + r1 * dc1dr1) - (dp + cr * du) * deig2dr1;
    dR2du1 = - cr * eig2 - (dp + cr * du) * deig2du1;
    dR2dp1 = - eig2 * (1.0 + du * r1 * dc1dp1) - (dp + cr * du) * deig2dp1;
    dR2dr2 = - (dp + cr * du) * deig2dr2;
    dR2du2 = cr * eig2 - (dp + cr * du) * deig2du2;
    dR2dp2 = eig2 - (dp + cr * du) * deig2dp2;

//  dR3dr1 = c1 * du * eig3 + du * eig3 * r1 * dc1dr1 + (-dp + cr * du) * deig3dr1; 
//  dR3du1 = - cr * eig3 + (-dp + cr * du) * deig3du1;
//  dR3dp1 = - eig3 + du * eig3 * r1 * dc1dp1 + (-dp + cr * du) * deig3dp1;
//  dR3dr2 = (-dp + cr * du) * deig3dr2;
//  dR3du2 = - cr * eig3 + (-dp + cr * du) * deig3du2;
//  dR3dp2 = eig3 + (-dp + cr * du) * deig3dp2;

//  std::cout<<"dr3dr1"<<dR3dr1<<std::endl;
//  std::cout<<"dr3du1"<<dR3du1<<std::endl;
//  std::cout<<"dr3dp1"<<dR3dp1<<std::endl;
//  std::cout<<"dr3dr2"<<dR3dr2<<std::endl;
//  std::cout<<"dr3du2"<<dR3du2<<std::endl;
//  std::cout<<"dr3dp2"<<dR3dp2<<std::endl;

    dR3dr1 = eig3 * du * (c1 + r1 * dc1dr1) - (dp - cr * du) * deig3dr1; 
    dR3du1 = cr * eig3 - (dp - cr * du) * deig3du1;
    dR3dp1 = - eig3 - dp * deig3dp1 + du * r1 * eig3 * dc1dp1;
    dR3dr2 = - (dp - cr * du) * deig3dr2;
    dR3du2 = - cr * eig3 - (dp - cr * du) * deig3du2;
    dR3dp2 = eig3 - (dp - cr * du) * deig3dp2;

//  std::cout<<"dr3dr1"<<dR3dr1<<std::endl;
//  std::cout<<"dr3du1"<<dR3du1<<std::endl;
//  std::cout<<"dr3dp1"<<dR3dp1<<std::endl;
//  std::cout<<"dr3dr2"<<dR3dr2<<std::endl;
//  std::cout<<"dr3du2"<<dR3du2<<std::endl;
//  std::cout<<"dr3dp2"<<dR3dp2<<std::endl;
    // dp1/dt
    double dp1dt;
    double dp1dtdr1, dp1dtdu1, dp1dtdp1;
    double dp1dtdr2, dp1dtdu2, dp1dtdp2;
    if(u1 < c1)
    {
        dp1dt = 0;
        dp1dtdr1 = 0;
        dp1dtdu1 = 0;
        dp1dtdp1 = 0;
        dp1dtdr2 = 0;
        dp1dtdu2 = 0;
        dp1dtdp2 = 0;
    }
    else
    {
        dp1dt = (R2 + R3) / 2.0;
        dp1dtdr1 = (dR2dr1 + dR3dr1) / 2.0;
        dp1dtdu1 = (dR2du1 + dR3du1) / 2.0;
        dp1dtdp1 = (dR2dp1 + dR3dp1) / 2.0;
        dp1dtdr2 = (dR2dr2 + dR3dr2) / 2.0;
        dp1dtdu2 = (dR2du2 + dR3du2) / 2.0;
        dp1dtdp2 = (dR2dp2 + dR3dp2) / 2.0;
    }

    // drho1/dt
    double dr1dt;
    double dr1dtdr1, dr1dtdu1, dr1dtdp1;
    double dr1dtdr2, dr1dtdu2, dr1dtdp2;
    dr1dt = R1 + dp1dt / (c1 * c1);

    dr1dtdr1 = dR1dr1 + dp1dtdr1 / (c1 * c1) - 2.0 * dp1dt * dc1dr1 / pow(c1, 3);
    dr1dtdu1 = dR1du1 + dp1dtdu1 / (c1 * c1);
    dr1dtdp1 = dR1dp1 + dp1dtdp1 / (c1 * c1) - 2.0 * dp1dt * dc1dp1 / pow(c1, 3);
    dr1dtdr2 = dR1dr2 + dp1dtdr2 / (c1 * c1);
    dr1dtdu2 = dR1du2 + dp1dtdu2 / (c1 * c1);
    dr1dtdp2 = dR1dp2 + dp1dtdp2 / (c1 * c1);

    // du1/dt
    double du1dt;
    double du1dtdr1, du1dtdu1, du1dtdp1;
    double du1dtdr2, du1dtdu2, du1dtdp2;
    du1dt = (R2 - dp1dt) / (cr);

    du1dtdr1 = ( (dp1dt - R2) * r1 * dc1dr1
               + c1 * (dp1dt - R2 - r1 * dp1dtdr1 + r1 * dR2dr1) )
               / (cr * cr);
    du1dtdu1 = (dR2du1 - dp1dtdu1) / cr;
    du1dtdp1 = ( (dp1dt - R2) * dc1dp1 + c1 * (dR2dp1 - dp1dtdp1) ) / (cr * c1);
    du1dtdr2 = (dR2dr2 - dp1dtdr2) / cr;
    du1dtdu2 = (dR2du2 - dp1dtdu2) / cr;
    du1dtdp2 = (dR2dp2 - dp1dtdp2) / cr;

    // d(ru)1/dt
    double dru1dt;
    dru1dt = r1 * du1dt + u1 * dr1dt;
    double dru1dtdr1, dru1dtdu1, dru1dtdp1;
    double dru1dtdr2, dru1dtdu2, dru1dtdp2;
    dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
    dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
    dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
    dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
    dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
    dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

    // de1/dt
    double de1dt;
    de1dt = dp1dt * Cv / R + u1 * r1 * du1dt + uu * dr1dt / 2.0;
    double de1dtdr1, de1dtdu1, de1dtdp1;
    double de1dtdr2, de1dtdu2, de1dtdp2;

    de1dtdr1 = du1dt * u1 + dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1;
    de1dtdu1 = dr1dt * u1 + du1dt * r1 + dp1dtdu1 * Cv / R
              + uu * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1;
    de1dtdp1 = dp1dtdp1 * Cv / R + uu * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
    de1dtdr2 = dp1dtdr2 * Cv / R + u1 + r1 * u1 * du1dtdr2 + uu * dr1dtdr2 / 2.0;
    de1dtdu2 = dp1dtdu2 * Cv / R + u1 + r1 * u1 * du1dtdu2 + uu * dr1dtdu2 / 2.0;
    de1dtdp2 = dp1dtdp2 * Cv / R + u1 + r1 * u1 * du1dtdp2 + uu * dr1dtdp2 / 2.0;


    // BN
    BN[0] = dr1dt;
    std::cout<<"dr1dt outlet "<<BN[0]<<std::endl;
    std::cout<<"dp1dt outlet "<<dp1dt<<std::endl;
    std::cout<<"c2 outlet "<<c1 * c1<<std::endl;
    std::cout<<"R1"<<R1<<std::endl;
    std::cout<<"R2"<<R2<<std::endl;
    std::cout<<"R3"<<R3<<std::endl;
    BN[1] = dru1dt;
    BN[2] = de1dt;

    dbdwp[0] = dr1dtdr1;
    dbdwp[1] = dr1dtdu1;
    dbdwp[2] = dr1dtdp1;
    dbdwp[3] = dru1dtdr1;
    dbdwp[4] = dru1dtdu1;
    dbdwp[5] = dru1dtdp1;
    dbdwp[6] = de1dtdr1;
    dbdwp[7] = de1dtdu1;
    dbdwp[8] = de1dtdp1;

    std::cout.precision(17);
    // Get Transformation Matrix
    dWpdW(dwpdw, W, nx - 1);
    std::cout<<"dbodwo (1,1)"<<std::endl; 
    for(int row = 0; row < 1; row++)
    for(int col = 0; col < 1; col++)
    for(int k = 0; k < 3; k++)
    {
        std::cout<<"k = "<<k<<std::endl;
        std::cout<<dBodWo[row * 3 + col]<<std::endl;
        std::cout<<dbdwp[row * 3 + k] * dwpdw[k * 3 + col]<<std::endl;
        std::cout<<"dbdwp "<<dbdwp[row * 3 + k]<<" dwpdw "<<dwpdw[k * 3 + col]<<std::endl;
        dBodWo[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
        std::cout<<dBodWo[row * 3 + col]<<std::endl;
    }
    dbdwp[0] = dr1dtdr2;
    dbdwp[1] = dr1dtdu2;
    dbdwp[2] = dr1dtdp2;
    dbdwp[3] = dru1dtdr2;
    dbdwp[4] = dru1dtdu2;
    dbdwp[5] = dru1dtdp2;
    dbdwp[6] = de1dtdr2;
    dbdwp[7] = de1dtdu2;
    dbdwp[8] = de1dtdp2;

    std::cout<<"dbdwp"<<std::endl; 
    std::cout<<dbdwp[0]<<std::endl;
    std::cout<<dbdwp[1]<<std::endl;
    std::cout<<dbdwp[2]<<std::endl;
    std::cout<<dbdwp[3]<<std::endl;
    std::cout<<dbdwp[4]<<std::endl;
    std::cout<<dbdwp[5]<<std::endl;
    std::cout<<dbdwp[6]<<std::endl;
    std::cout<<dbdwp[7]<<std::endl;
    std::cout<<dbdwp[8]<<std::endl;

    // Get Transformation Matrix
    dWpdW(dwpdw, W, nx - 2);
    std::cout<<"dwpdw"<<std::endl; 
    std::cout<<dwpdw[0]<<std::endl;
    std::cout<<dwpdw[1]<<std::endl;
    std::cout<<dwpdw[2]<<std::endl;
    std::cout<<dwpdw[3]<<std::endl;
    std::cout<<dwpdw[4]<<std::endl;
    std::cout<<dwpdw[5]<<std::endl;
    std::cout<<dwpdw[6]<<std::endl;
    std::cout<<dwpdw[7]<<std::endl;
    std::cout<<dwpdw[8]<<std::endl;
    for(int row = 0; row < 3; row++)
    for(int col = 0; col < 3; col++)
    for(int k = 0; k < 3; k++)
        dBodWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
    std::cout<<"dBodWd"<<std::endl; 
    std::cout<<dBodWd[0]<<std::endl;
    std::cout<<dBodWd[1]<<std::endl;
    std::cout<<dBodWd[2]<<std::endl;
    std::cout<<dBodWd[3]<<std::endl;
    std::cout<<dBodWd[4]<<std::endl;
    std::cout<<dBodWd[5]<<std::endl;
    std::cout<<dBodWd[6]<<std::endl;
    std::cout<<dBodWd[7]<<std::endl;
    std::cout<<dBodWd[8]<<std::endl;

    // *********************
    // INLET JACOBIANS
    // *********************
    // Subsonic Inlet
    if(u[0] < c[0])
    {
        i1 = 0;
        i2 = 1;
        r1 = rho[i1];
        r2 = rho[i2];
        p1 = p[i1];
        p2 = p[i2];
        u1 = u[i1];
        u2 = u[i2];
        c1 = c[i1];
        c2 = c[i2];
        t1 = T[i1];

        // Shorthand
        drho = r2 - r1;
        dp = p2 - p1;
        du = u2 - u1;
        cr = r1 * c1;
        uu = u1 * u1;
        gamr = (gam - 1) / (gam + 1);
        fu = 1 - gamr * u1 * u1 / a2;

        // Speed of Sound
        dc1dr1 = - p1 * gam / (2.0 * cr * r1);
        dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
        dc1dp1 = gam / (2.0 * cr);
        dc2dp2 = gam / (2.0 * c2 * r2);

        // Eigenvalue
        eig1 = (u1 + u2) / 2.0;
        eig3 = eig1 - (c1 + c2) / 2.0;

        deig1du1 = 1.0 / 2.0;
        deig1du2 = 1.0 / 2.0;

        deig3dr1 = - dc1dr1 / 2.0;
        deig3du1 = deig1du1;
        deig3dp1 = - dc1dp1 / 2.0;
        deig3dr2 = - dc2dr2 / 2.0;
        deig3du2 = deig1du2;
        deig3dp2 = - dc2dp2 / 2.0;

        // Riemann Invariants
        R3 = - eig3 * (dp - cr * du);

        dR3dr1 = eig3 * (c1 * du + du * r1 * dc1dr1) - (dp - cr * du) * deig3dr1; 
        dR3du1 = - cr * eig3 - (dp - cr * du) * deig3du1;
        dR3dp1 = eig3 * (du * r1 * dc1dp1 + 1.0) - (dp - cr * du) * deig3dp1;
        dR3dr2 = -(dp - cr * du) * deig3dr2;
        dR3du2 = cr * eig3 - (dp - cr * du) * deig3du2;
        dR3dp2 = - eig3 - (dp - cr * du) * deig3dp2;


        // dp1
        double dp1du1, dp1du1du1;
//      dp1du1 = -2 * ptin * pow(fu, 1 / (gam - 1)) * gamr * u1 / a2 * (gam / (gam - 1));
//      dp1du1du1 = 2 * gamr * ptin * pow(fu, 1 / (gam - 1)) * gam
//                  * (a2 - a2 * gam + gamr * uu * (gam + 1))
//                  / (a2 * (a2 - gamr * uu) * pow(gam - 1, 2));
        // Same Values
        dp1du1 = -2.0 * gamr * ptin * u1 * pow(fu, 1.0 / (gam - 1.0)) * gam
                 / (a2 * (gam - 1.0));
        dp1du1du1 = 2 * gamr * ptin * pow(fu, gam/(gam - 1.0)) * gam
                    * (a2 - a2 * gam + gamr * uu * (gam + 1))
                    / pow((a2 - gamr * uu) * (gam - 1.0), 2);

        // du1
        du1dt = R3 / (dp1du1 - cr);
        du1dtdr1 = dR3dr1 / (dp1du1 - cr)
                   + R3 * (c1 + r1 * dc1dr1) / pow((dp1du1 - cr), 2);
        du1dtdu1 = dR3du1 / (dp1du1 - cr)
                   - R3 * dp1du1du1 / pow((dp1du1 - cr), 2);
        du1dtdp1 = (R3 * r1 * dc1dp1) / pow((dp1du1 - cr), 2)
                   + dR3dp1 / (dp1du1 - cr);
        du1dtdr2 = dR3dr2 / (dp1du1 - cr);
        du1dtdu2 = dR3du2 / (dp1du1 - cr);
        du1dtdp2 = dR3dp2 / (dp1du1 - cr);

        // dp1/dt
        dp1dt  = dp1du1 * du1dt;
        dp1dtdr1 = dp1du1 * du1dtdr1;
        dp1dtdu1 = du1dt * dp1du1du1 + dp1du1 * du1dtdu1;
        dp1dtdp1 = dp1du1 * du1dtdp1;
        dp1dtdr2 = dp1du1 * du1dtdr2;
        dp1dtdu2 = dp1du1 * du1dtdu2;
        dp1dtdp2 = dp1du1 * du1dtdp2;

        // dt1
        double dt1dp1, dt1dp1dp1;
        dt1dp1 = Ttin / ptin * (gam - 1.0) / gam * pow(p1 / ptin, - 1.0 / gam);
        dt1dp1dp1 = - Ttin * (gam - 1.0) / pow((gam * ptin), 2)
                    * pow(p1 / ptin, - (1.0 + gam) / gam);

        // dr1/dt
        double dr1dp1, dr1dp1dp1;
        dr1dp1 = 1.0 / (R * t1) - p1 * dt1dp1 / (R * t1 * t1);
        dr1dp1dp1 = 2.0 * p1 * dt1dp1 * dt1dp1 / (R * t1 * t1 * t1)
                    - 2.0 * dt1dp1 / (R * t1 * t1)
                    - p1 * dt1dp1dp1 / (R * t1 * t1);
        dr1dt = dr1dp1 * dp1dt;

        dr1dtdr1 = dr1dp1 * dp1dtdr1;
        dr1dtdu1 = dr1dp1 * dp1dtdu1;
        dr1dtdp1 = dr1dp1 * dp1dtdp1 + dp1dt * dr1dp1dp1;
        dr1dtdr2 = dr1dp1 * dp1dtdr2;
        dr1dtdu2 = dr1dp1 * dp1dtdu2;
        dr1dtdp2 = dr1dp1 * dp1dtdp2;

        // dru1/dt
        dru1dt = r1 * du1dt + u1 * dr1dt;

        dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
        dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
        dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
        dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
        dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
        dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

        // de1/dt
        de1dt = dp1dt * Cv / R + r1 * u1 * du1dt + uu * dr1dt / 2.0;
        de1dtdr1 = dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1 
                   + du1dt * u1;
        de1dtdu1 = dp1dtdu1 * Cv / R + uu * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1 
                   + du1dt * r1 + dr1dt * u1;
        de1dtdp1 = dp1dtdp1 / (gam - 1) + uu * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
        de1dtdr2 = dp1dtdr2 / (gam - 1) + uu * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
        de1dtdu2 = dp1dtdu2 / (gam - 1) + uu * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
        de1dtdp2 = dp1dtdp2 / (gam - 1) + uu * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

        // Same values
        //de1dtdr1 = Cv * dp1dt * dt1dp1 + du1dt * u1 + Cv * dt1dp1 * r1 * dp1dtdr1
        //           + (Cv * t1 + uu / 2) * dr1dtdr1 + u1 * r1 * du1dtdr1;
        //de1dtdu1 = dr1dt * u1 * du1dt * r1 + Cv * dt1dp1 * r1 * dp1dtdu1
        //           + (Cv * t1 + uu / 2) * dr1dtdu1 + u1 * r1 * du1dtdu1;
        //de1dtdp1 = Cv * dt1dp1 * r1 * dp1dtdp1 + (Cv * t1 + uu / 2) * dr1dtdp1
        //           + Cv * dp1dt * r1 * dt1dp1dp1 + u1 * r1 * du1dtdp1
        //           + Cv * dr1dt * dt1dp1;
        //de1dtdr2 = Cv * dt1dp1 * r1 * dp1dtdr2 + (Cv * t1 + uu / 2) * dr1dtdr2
        //           + u1 * r1 * du1dtdr2;
        //de1dtdu2 = Cv * dt1dp1 * r1 * dp1dtdu2 + (Cv * t1 + uu / 2) * dr1dtdu2
        //           + u1 * r1 * du1dtdu2;
        //de1dtdp2 = Cv * dt1dp1 * r1 * dp1dtdp2 + (Cv * t1 + uu / 2) * dr1dtdp2
        //           + u1 * r1 * du1dtdp2;

        // B1
        B1[0] = dr1dt;
        B1[1] = dru1dt;
        B1[2] = de1dt;

        dbdwp[0] = dr1dtdr1;
        dbdwp[1] = dr1dtdu1;
        dbdwp[2] = dr1dtdp1;
        dbdwp[3] = dru1dtdr1;
        dbdwp[4] = dru1dtdu1;
        dbdwp[5] = dru1dtdp1;
        dbdwp[6] = de1dtdr1;
        dbdwp[7] = de1dtdu1;
        dbdwp[8] = de1dtdp1;

        // Get Transformation Matrix
        dWpdW(dwpdw, W, 0);
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        for(int k = 0; k < 3; k++)
            dBidWi[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
        
        dbdwp[0] = dr1dtdr2;
        dbdwp[1] = dr1dtdu2;
        dbdwp[2] = dr1dtdp2;
        dbdwp[3] = dru1dtdr2;
        dbdwp[4] = dru1dtdu2;
        dbdwp[5] = dru1dtdp2;
        dbdwp[6] = de1dtdr2;
        dbdwp[7] = de1dtdu2;
        dbdwp[8] = de1dtdp2;

        // Get Transformation Matrix
        dWpdW(dwpdw, W, 1);
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        for(int k = 0; k < 3; k++)
            dBidWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
    }
    // Supersonic Inlet
    else
    {
        B1[0] = 0;
        B1[1] = 0;
        B1[2] = 0;
        for(int i = 0; i < 9; i++)
        {
            dBidWi[i] = 0;
            dBidWd[i] = 0;
            if(i % 4 == 0)
                dBidWi[i] = 1;
        }
    }
}

SparseMatrix<double> buildAMatrix(std::vector <double> Ap,
                                  std::vector <double> An,
                                  std::vector <double> dBidWi,
                                  std::vector <double> dBidWd,
                                  std::vector <double> dBodWd,
                                  std::vector <double> dBodWo,
                                  std::vector <double> dQdW,
                                  std::vector <double> dx,
                                  std::vector <double> dt,
                                  std::vector <double> S,
                                  double Min)
{
    SparseMatrix<double> matA(3 * nx, 3 * nx);
    int dwi, psii;
    int k, ri, ci;
    double val;
    // Input 4 lines where BC Jacobians occur
    // psi(1), psi(2), psi(n-1), psi(n)
    for(int row = 0; row < 3; row++)
    {
        for(int col = 0; col < 3; col++)
        {
            k = row * 3 + col;
            // dw0, psi0
            dwi = 0;
            psii = 0;
            ri = dwi * 3 + row;
            ci = psii * 3 + col;

            val = - dBidWi[k];
            matA.insert(ri, ci) = val;

            // dw1, psi0
            dwi = 1;
            psii = 0;
            ri = dwi * 3 + row;
            ci = psii * 3 + col;

            val = - dBidWd[k];
            matA.insert(ri, col) = val;

            // dw(nx-1), psi(nx-1)
            dwi = nx - 1;
            psii = nx - 1;
            ri = dwi * 3 + row;
            ci = psii * 3 + col;

            val = - dBodWo[k];
            matA.insert(ri, ci) = val;

            // dw(nx-2), psi(nx-1)
            dwi = nx - 2;
            psii = nx - 1;
            ri = dwi * 3 + row;
            ci = psii * 3 + col;

            val = - dBodWd[k];
            std::cout<<val<<std::endl;
            matA.insert(ri, ci) = val;
        }
    }
    for(int dwi = 0; dwi < nx; dwi++)
    {
        psii = dwi - 1;
        if(psii >= 1)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                ri = dwi * 3 + row;
                ci = psii * 3 + col;

                val = An[dwi * 9 + k] * S[psii];

                matA.insert(ri, ci) = val;
            }
        }

        psii = dwi;
        if(psii >= 1 && psii <= nx - 2)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                ri = dwi * 3 + row;
                ci = psii * 3 + col;

                val = Ap[dwi * 9 + k] * S[psii + 1];
                val -= An[dwi * 9 + k] * S[psii];
                if(row == 1) // Remember it is the transposed dQdW
                {
                    val -= dQdW[dwi * 3 + col];
                }

                matA.insert(ri, ci) = val;
            }
        }

        psii = dwi + 1;
        if(psii <= nx - 2)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                ri = dwi * 3 + row;
                ci = psii * 3 + col;

                val = - Ap[dwi * 9 + k] * S[psii + 1];

                matA.insert(ri, ci) = val;
            }
        }
    }
    if(Min > 1.0)
    {
        // Supersonic Inlet, don't solve for psi(0)
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            // dw(0), psi(1)
            dwi = 0;
            psii = 1;
            ri = dwi * 3 + row;
            ci = psii * 3 + col;
    
            matA.coeffRef(ri, ci) = 0;
        }
    }
    // Extrapolate Outlet for now... 
//  for(int row = 0; row < 3; row++)
//  for(int col = 0; col < 3; col++)
//  {
//      // dw(nx-1), psi(nx-1)
//      dwi = nx - 1;
//      psii = nx - 1;
//      ri = dwi * 3 + row;
//      ci = psii * 3 + col;

//      val = - dBodWo[k] * (S[nx - 1] + S[nx]) / 2;
//      matA.coeffRef(ri, ci) = matA.coeffRef(ri - 3, ci - 3);

//      matA.coeffRef(ri, ci) = 0;
//      if(row == col)
//          matA.coeffRef(ri, ci) = -1;

//      // dw(nx-2), psi(nx-1)
//      dwi = nx - 2;
//      psii = nx - 1;
//      ri = dwi * 3 + row;
//      ci = psii * 3 + col;

//      matA.coeffRef(ri, ci) = matA.coeffRef(ri, ci - 3);
//      matA.coeffRef(ri, ci) = 0;
//      if(row == col)
//          matA.coeffRef(ri, ci) = 1;
//  }
    return matA;
}

VectorXd buildbMatrix(std::vector <double> dIcdW)
{
    VectorXd matb(3 * nx);

    for(int i = 0; i < nx; i++)
    for(int k = 0; k < 3; k++)
        matb(i * 3 + k) = -dIcdW[i * 3 + k];
    
    return matb;
}

void evalddtdW(std::vector <double> &ddtdW,
                 std::vector <double> rho,
                 std::vector <double> u,
                 std::vector <double> p,
                 std::vector <double> c,
                 std::vector <double> dx,
                 std::vector <double> S)
{
    double dStepdr, dStepdru, dStepde;
    for(int i = 0; i < nx; i++)
    {
        dStepdr  = ( 2 * p[i] * gam 
                   + u[i] * rho[i] * ( u[i] * (1 - gam) * gam + 4 * c[i] ) )
                   / ( 4 * pow(u[i] + c[i], 2) * c[i] * rho[i] * rho[i] );
        dStepdru = ( u[i] * gam * (gam - 1) - 2 * c[i] )
                   / ( 2 * pow(u[i] + c[i], 2) * c[i] * rho[i] );
        dStepde  = (1 - gam) * c[i]
                   / ( 2 * p[i] * pow(u[i] + c[i], 2) );

        ddtdW[i * 3 + 0] = CFL * dx[i] * dStepdr;
        ddtdW[i * 3 + 1] = CFL * dx[i] * dStepdru;
        ddtdW[i * 3 + 2] = CFL * dx[i] * dStepde;
    }
}
