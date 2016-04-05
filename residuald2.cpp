#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"
#include"residuald1.h"
#include"BCDerivatives.h"

using namespace Eigen;

std::vector <MatrixXd> evalddFdWdW(double rho, double u, double p);

void evalddFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W);

void evalddScalarFluxdWdW_FD(
    std::vector <MatrixXd> &ddFluxdWdW1,
    int ii, int jj,
    std::vector <double> W);


std::vector <MatrixXd> evalddQdWdW(std::vector <double> W);

std::vector <SparseMatrix <double> > evalddRdWdW_FD(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <SparseMatrix <double> > ddRdWdW_FD(3 * nx);

    SparseMatrix <double> ddRidWdW_FD(3 * nx, 3 * nx);
    for(int Ri = 0; Ri < 3 * nx; Ri++)
        ddRdWdW_FD[Ri] = ddRidWdW_FD;

    double Resi0, Resi1, Resi2, Resi3, Resi4;
    std::vector <double> Resitemp(3 * nx, 0);
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 1e-3;
    double pertWi, pertWj;
    int Rik, Rikp;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int k = 0; k < 3; k++)
        {
            Rik = Ri * 3 + k;
            Rikp = (Ri + 1) * 3 + k;

            for(int Wi = (Ri - 1) * 3; Wi <= (Ri + 1) * 3 + 2; Wi++)
            {
                pertWi = W[Wi] * h;
                for(int Wj = (Ri - 1) * 3; Wj <= (Ri + 1) * 3 + 2; Wj++)
                {
                    pertWj = W[Wj] * h;
                    if(Wi != Wj)
                    {
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        // R1
                        Wd[Wi] = W[Wi] + pertWi;
                        Wd[Wj] = W[Wj] + pertWj;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi1 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R2
                        Wd[Wi] = W[Wi] + pertWi;
                        Wd[Wj] = W[Wj] - pertWj;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi2 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R3
                        Wd[Wi] = W[Wi] - pertWi;
                        Wd[Wj] = W[Wj] + pertWj;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi3 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R4
                        Wd[Wi] = W[Wi] - pertWi;
                        Wd[Wj] = W[Wj] - pertWj;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi4 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                     / (4 * pertWi * pertWj);

                        // Symmetry
//                      ddRdWdW_FD[Rik].insert(Wj, Wi) = ddRdWdW_FD[Rik].coeffRef(Wi, Wj);
                        // Reset Wd
                        Wd[Wi] = W[Wi];
                        Wd[Wj] = W[Wj];
                    }
                    else
                    {
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi0 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R1
                        Wd[Wi] = W[Wi] + 2.0 * pertWi;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi1 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R2
                        Wd[Wi] = W[Wi] + pertWi;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi2 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R3
                        Wd[Wi] = W[Wi] - pertWi;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi3 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        // R4
                        Wd[Wi] = W[Wi] - 2.0 * pertWi;

                        WtoQ(Wd, Q, S);
                        getFlux(Flux, Wd);

                        Resi4 = Flux[Rikp] * S[Ri + 1] - Flux[Rik] * S[Ri] - Q[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                                                     / (12 * pertWi * pertWj);
                        // Reset Wd
                        Wd[Wi] = W[Wi];
                        Wd[Wj] = W[Wj];
                    }
                }// Wj Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    //Inlet
    for(int Ri = 0; Ri < 1; Ri++)
    {
        for(int k = 0; k < 3; k++)
        {
            Rik = Ri * 3 + k;
            Rikp = (Ri + 1) * 3 + k;

            for(int Wi = Ri * 3; Wi <= (Ri + 1) * 3 + 2; Wi++)
            {
                pertWi = W[Wi] * h;
                for(int Wj = Ri * 3; Wj <= (Ri + 1) * 3 + 2; Wj++)
                {
                    pertWj = W[Wj] * h;
                    if(Wi != Wj)
                    {
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        // R1
                        Wd[Wi] = W[Wi] + pertWi;
                        Wd[Wj] = W[Wj] + pertWj;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi1 = Resitemp[Rik];

                        // R2
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] + pertWi;
                        Wd[Wj] = W[Wj] - pertWj;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi2 = Resitemp[Rik];

                        // R3
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] - pertWi;
                        Wd[Wj] = W[Wj] + pertWj;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi3 = Resitemp[Rik];

                        // R4
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] - pertWi;
                        Wd[Wj] = W[Wj] - pertWj;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi4 = Resitemp[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                     / (4 * pertWi * pertWj);
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                    }
                    else
                    {
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi0 = Resitemp[Rik];

                        // R1
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] + 2.0 * pertWi;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi1 = Resitemp[Rik];

                        // R2
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] + pertWi;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi2 = Resitemp[Rik];

                        // R3
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] - pertWi;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi3 = Resitemp[Rik];

                        // R4
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                        Wd[Wi] = W[Wi] - 2.0 * pertWi;

                        inletBC(Wd, Resitemp, 1, 1);
                        Resi4 = Resitemp[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) =
                            (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                            / (12 * pertWi * pertWj);
                        // Reset Wd
                        for(int m = 0; m < 3 * nx; m++)
                            Wd[m] = W[m];
                    }
                }// Wj Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    //Outlet
    int Ri = nx-1;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;
        Rikp = (Ri + 1) * 3 + Rk;

        for(int Wi = (Ri - 1) * 3; Wi <= Ri * 3 + 2; Wi++)
        {
            pertWi = W[Wi] * h;
            for(int Wj = (Ri - 1) * 3; Wj <= Ri * 3 + 2; Wj++)
            {
                pertWj = W[Wj] * h;
                if(Wi != Wj)
                {
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    // R1
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi1 = Resitemp[Rik];

                    // R2
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi2 = Resitemp[Rik];

                    // R3
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi3 = Resitemp[Rik];

                    // R4
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi4 = Resitemp[Rik];

                    ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                 / (4 * pertWi * pertWj);
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                }
                else
                {
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi0 = Resitemp[Rik];

                    // R1
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] + 2.0 * pertWi;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi1 = Resitemp[Rik];

                    // R2
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] + pertWi;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi2 = Resitemp[Rik];

                    // R3
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] - pertWi;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi3 = Resitemp[Rik];

                    // R4
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                    Wd[Wi] = W[Wi] - 2.0 * pertWi;

                    outletBC(Wd, Resitemp, 1, 1);
                    Resi4 = Resitemp[Rik];

                    ddRdWdW_FD[Rik].insert(Wi, Wj) =
                        (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset Wd
                    for(int m = 0; m < 3 * nx; m++) Wd[m] = W[m];
                } // If not diagonal
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
    return ddRdWdW_FD;
}

// Calculates Residual Hessian
std::vector < SparseMatrix <double> > evalddRdWdW(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <SparseMatrix <double> > ddRdWdW(3 * nx);
    SparseMatrix <double> dummySMat(3 * nx, 3 * nx);

    for(int Ri = 0; Ri < nx; Ri++)
    {
        for(int Rk = 0; Rk < 3; Rk++)
        {
            ddRdWdW[Ri * 3 + Rk] = dummySMat;
        }
    }
    // ********************************
    // BC Residual Hessian
    // ********************************
    // Inlet and Outlet Residual Hessian
    // ddRindWdW contains 3 matrices: dr1/dWdW, dr2/dWdW, dr3/dWdW
    // 6 state variables affect the inlet/outlet residuals
    std::vector <MatrixXd> ddRindWdW(3), ddRindWdW_FD(3);
    std::vector <MatrixXd> ddRoutdWdW(3), ddRoutdWdW_FD(3);
    MatrixXd dummyMat(6, 6);
    for(int Rk = 0; Rk < 3; Rk++)
    {
        ddRindWdW[Rk] = dummyMat;
        ddRoutdWdW[Rk] = dummyMat;
        ddRindWdW_FD[Rk] = dummyMat;
        ddRoutdWdW_FD[Rk] = dummyMat;
    }
    HessianBC(W, ddRindWdW, ddRoutdWdW);
//  HessianBCprim_FD(W, ddRindWdW_FD, ddRoutdWdW_FD);
//  HessianBC_FD(W, ddRindWdW, ddRoutdWdW);

//  std::cout<<std::endl;
//  std::cout<<"FD ddRdWdW_in 0"<<std::endl;
//  std::cout<<ddRindWdW_FD[0]<<std::endl;
//  std::cout<<"AN ddRdWdW_in 0"<<std::endl;
//  std::cout<<ddRindWdW[0]<<std::endl;
//  std::cout<<"FD - AN ddRdWdW_in 0"<<std::endl;
//  std::cout<<ddRindWdW_FD[0]+ddRindWdW[0]<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<"difference in norm"<<std::endl;
//  std::cout<<(ddRindWdW_FD[0].norm()-ddRindWdW[0].norm())/ddRindWdW[0].norm()<<std::endl;

//  std::cout<<std::endl;
//  std::cout<<"FD ddRdWdW_in 1"<<std::endl;
//  std::cout<<ddRindWdW_FD[1]<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<"AN ddRdWdW_in 1"<<std::endl;
//  std::cout<<ddRindWdW[1]<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<"FD - AN ddRdWdW_in 1"<<std::endl;
//  std::cout<<ddRindWdW_FD[1]+ddRindWdW[1]<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<"difference in norm"<<std::endl;
//  std::cout<<(ddRindWdW_FD[1].norm()-ddRindWdW[1].norm())/ddRindWdW[1].norm()<<std::endl;

//  std::cout<<std::endl;
//  std::cout<<"FD ddRdWdW_in 2"<<std::endl;
//  std::cout<<ddRindWdW_FD[2]<<std::endl;
//  std::cout<<"AN ddRdWdW_in 2"<<std::endl;
//  std::cout<<ddRindWdW[2]<<std::endl;
//  std::cout<<"FD - AN ddRdWdW_in 2"<<std::endl;
//  std::cout<<ddRindWdW_FD[2]+ddRindWdW[2]<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<"difference in norm"<<std::endl;
//  std::cout<<(ddRindWdW_FD[2].norm()-ddRindWdW[2].norm())/ddRindWdW[2].norm()<<std::endl;

    for(int Rk = 0; Rk < 3; Rk++)
    {
      for(int row = 0; row < 6; row++)
      {
        for(int col = 0; col < 6; col++)
        {
          // Inlet
          ddRdWdW[Rk].insert(row, col) = -ddRindWdW[Rk](row, col);
          // Outlet
          int rowi = (nx - 2) * 3 + row;
          int coli = (nx - 2) * 3 + col;
          ddRdWdW[(nx - 1) * 3 + Rk].insert(rowi, coli) = -ddRoutdWdW[Rk](row, col);
        }
      }
    }

    // ********************************
    // Domain Residual Hessian
    // Validated with FD
    // ********************************
    // Evaluate ddFluxdWdW
    // ddFluxdWdW1 = Flux(i) / W(i-1, i-1)
    // ddFluxdWdW2 = Flux(i) / W(i  , i  )
    // ddFluxdWdW3 = Flux(i) / W(i-1, i  )
    // ddFluxdWdW4 = TRANSPOSE(ddFluxdWdW3)
    std::vector <MatrixXd> ddFluxdWdW1(3 * (nx + 1)),
                           ddFluxdWdW2(3 * (nx + 1)),
                           ddFluxdWdW3(3 * (nx + 1));
    evalddFluxdWdW(ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3, W);

    // Evaluate ddQdWdW
    std::vector <MatrixXd> ddQdWdW(3 * nx);
    ddQdWdW = evalddQdWdW(W);

    int Rik, Rikp;
    int Wi1, Wi2;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int Rk = 0; Rk < 3; Rk++)
        {
            Rik = Ri * 3 + Rk;
            Rikp = (Ri + 1) * 3 + Rk;

            // R(i) = Flux(i + 1) * S(i + 1) - Flux(i) * S(i) - Q(i)

            // ***************
            // Diagonal Terms
            // ***************

            //ddR(i) / dW(i-1)dW(i-1)
            // = - ddFlux(i) / dW(i-1)dW(i-1) * S(i)
            Wi1 = Ri - 1;
            Wi2 = Ri - 1;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                            -ddFluxdWdW1[Rik].coeffRef(row, col) * S[Ri];
            }
//          if(Rik == 25*3)
//          {
//              std::cout<<ddRdWdW[Rik]<<std::endl;
//              std::cout<<ddFluxdWdW1[Rik]<<std::endl;
//          }

            //ddR(i) / dW(i)dW(i)
            // = ddFlux(i + 1) / dW(i)dW(i) * S(i + 1)
            //   - ddFlux(i) / dW(i)dW(i) * S(i)
            //   - ddQ(i) / dW(i)dW(i) * dS
            Wi1 = Ri;
            Wi2 = Ri;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                            ddFluxdWdW1[Rikp].coeffRef(row, col) * S[Ri + 1]
                            - ddFluxdWdW2[Rik].coeffRef(row, col) * S[Ri]
                            - ddQdWdW[Rik].coeffRef(row, col) * (S[Ri + 1] - S[Ri]);
            }

            //ddR(i) / dW(i+1)dW(i+1)
            // = ddFlux(i+1) / dW(i+1)dW(i+1) * S(i+1)
            Wi1 = Ri + 1;
            Wi2 = Ri + 1;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                            ddFluxdWdW2[Rikp].coeffRef(row, col) * S[Ri + 1];
            }

            // *******************
            // Off-Diagonal Terms
            // *******************

            // R(i) = Flux(i + 1) * S(i + 1) - Flux(i) * S(i) - Q(i)

            //ddR(i) / dW(i-1)dW(i)
            // = ddFlux(i) / dW(i-1)dW(i) * S(i)
            Wi1 = Ri - 1;
            Wi2 = Ri;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                    -ddFluxdWdW3[Rik].coeffRef(row, col) * S[Ri];
            }
            //ddR(i) / dW(i)dW(i-1)
            // = (ddFlux(i) / dW(i-1)dW(i) * S(i)).transpose()
            Wi1 = Ri;
            Wi2 = Ri - 1;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                    -ddFluxdWdW3[Rik].coeffRef(col, row) * S[Ri];
            }

            //ddR(i) / dW(i)dW(i+1)
            // = ddFlux(i) / dW(i)dW(i+1) * S(i+1)
            Wi1 = Ri;
            Wi2 = Ri + 1;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                    ddFluxdWdW3[Rikp].coeffRef(row, col) * S[Ri + 1];
            }
            //ddR(i) / dW(i+1)dW(i) Symmetry
            // = (ddFlux(i) / dW(i)dW(i+1) * S(i+1)).transpose()
            Wi1 = Ri + 1;
            Wi2 = Ri;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                    ddFluxdWdW3[Rikp].coeffRef(col, row) * S[Ri + 1];
            }

            //ddR(i) / dW(i-1)dW(i+1) = 0
            //ddR(i) / dW(i+1)dW(i-1) = 0
        }

    }
    return ddRdWdW;
}

std::vector <SparseMatrix <double> > evalddRdWdS_FD(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);

    SparseMatrix <double> ddRidWdS(3 * nx, nx + 1);
    for(int Ri = 0; Ri < 3 * nx; Ri++)
        ddRdWdS[Ri] = ddRidWdS;

    double Resi1, Resi2, Resi3, Resi4;
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Sd(nx + 1, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 1e-4;
    double pertW, pertS;
    int Rik, Rikp;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int k = 0; k < 3; k++)
        {
            Rik = Ri * 3 + k;
            Rikp = (Ri + 1) * 3 + k;

            for(int Wi = (Ri - 1) * 3; Wi <= (Ri + 1) * 3 + 2; Wi++)
            {
                pertW = W[Wi] * h;
                for(int Si = Ri; Si <= Ri+1; Si++)
                {
                    pertS = S[Si] * h;
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];
                    for(int m = 0; m < nx + 1; m++)
                        Sd[m] = S[m];
                    // R1
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] + pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi1 = Flux[Rikp] * Sd[Ri + 1] - Flux[Rik] * Sd[Ri] - Q[Rik];

                    // R2
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] - pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi2 = Flux[Rikp] * Sd[Ri + 1] - Flux[Rik] * Sd[Ri] - Q[Rik];

                    // R3
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] - pertW;
                    Sd[Si] = S[Si] + pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi3 = Flux[Rikp] * Sd[Ri + 1] - Flux[Rik] * Sd[Ri] - Q[Rik];

                    // R4
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] - pertW;
                    Sd[Si] = S[Si] - pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi4 = Flux[Rikp] * Sd[Ri + 1] - Flux[Rik] * Sd[Ri] - Q[Rik];

                    ddRdWdS[Rik].insert(Wi, Si) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                 / (4 * pertW * pertS);
                }// Si Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    return ddRdWdS;
}

std::vector <SparseMatrix <double> > evalddRdWdS(
    std::vector <double> W,
    std::vector <double> S)
{
    // Allocate ddRdWdS Sparse Matrix
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);
    SparseMatrix <double> ddRidWdS(3 * nx, nx + 1);
    ddRidWdS.reserve(9 * 2); // 9 State Variables and 2 Areas affect the Residual
    for(int Ri = 0; Ri < 3 * nx; Ri++)
        ddRdWdS[Ri] = ddRidWdS;

    // Get Jacobians and Fluxes
    std::vector <double> Ap(nx * 3 * 3, 0), An(nx * 3 * 3, 0);
    if(FluxScheme == 1) ScalarJac(W, Ap, An);
    // Evaluate dpdW
    std::vector <double> dpdW(3 * nx, 0);
    dpdW = evaldpdW(W, S);

    int Wi, Si;
    int Rik, Wik;
    double val;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;
        // Positive Jacobian of Left Incoming Flux
        Si = Ri;
        Wi = Ri - 1;
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            val = -Ap[Wi * 9 + Rk * 3 + Wk];
            ddRdWdS[Rik].insert(Wik, Si) = val;
        }
        // Negative Jacobian of Left Incoming Flux
        Si = Ri;
        Wi = Ri;
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            val = -An[Wi * 9 + Rk * 3 + Wk];
            if(Rk == 1)
            {
                val -= -dpdW[Wi * 3 + Wk];
            }
            ddRdWdS[Rik].insert(Wik, Si) = val;
        }
        // Positive Jacobian of Right Incoming Flux
        Si = Ri + 1;
        Wi = Ri;
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            val = Ap[Wi * 9 + Rk * 3 + Wk];
            if(Rk == 1)
            {
                val -= dpdW[Wi * 3 + Wk];
            }
            ddRdWdS[Rik].insert(Wik, Si) = val;
        }
        // Negative Jacobian of Right Incoming Flux
        Si = Ri + 1;
        Wi = Ri + 1;
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            val = An[Wi * 9 + Rk * 3 + Wk];
            ddRdWdS[Rik].insert(Wik, Si) = val;
        }
    }// Rk Loop
    }// Ri Loop
    return ddRdWdS;
}

std::vector <MatrixXd> evalddQdWdW(std::vector <double> W)
{
    std::vector <MatrixXd> ddQdWdW(3 * nx);
    MatrixXd ddQidWdW(3, 3);
    ddQidWdW.setZero();
    for(int i = 0; i < nx; i++)
    {
        ddQdWdW[i * 3 + 0] = ddQidWdW;
        ddQdWdW[i * 3 + 2] = ddQidWdW;
    }
    double rho, u;
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        ddQidWdW(0, 0) = u * u * (1.0 - gam) / rho;
        ddQidWdW(1, 1) = (1.0 - gam) / rho;
        ddQidWdW(0, 1) = u * (gam - 1.0) / rho;
        ddQidWdW(1, 0) = ddQidWdW(0, 1);
        // (:,3) = 0
        // (3,:) = 0
        ddQdWdW[i * 3 + 1] = ddQidWdW;
    }
    return ddQdWdW;
}
std::vector <MatrixXd> evalddFdWdW(double rho, double u, double p)
{
    std::vector <MatrixXd> ddFdWdW(3);
    MatrixXd ddfidWdW(3, 3);
    // F1
    ddfidWdW.setZero();
    ddFdWdW[0] = ddfidWdW;
    // F2
    ddfidWdW(0, 0) = -(u * u * (gam - 3.0) / rho);
    ddfidWdW(1, 1) = (3.0 - gam) / rho;
    ddfidWdW(0, 1) = (u * (gam - 3.0) / rho);
    ddfidWdW(1, 0) = ddfidWdW(0, 1);
    ddFdWdW[1] = ddfidWdW;
    // F3
    ddfidWdW(0, 0) = ( 2.0 * p * u * gam
                    + rho * pow(u, 3.0)
                    * (- 2.0 * gam * gam + 5.0 * gam - 3.0) )
                    / (pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 1) = 3.0 * u * (1.0 - gam) / rho;
    ddfidWdW(2, 2) = 0;

    ddfidWdW(0, 1) = (- 2.0 * p * gam
                    + rho * pow(u, 2.0)
                    * (5.0 * gam * gam - 11.0 * gam + 6.0) )
                    / (2.0 * pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 0) = ddfidWdW(0, 1);

    ddfidWdW(0, 2) = -u * gam / rho;
    ddfidWdW(2, 0) = ddfidWdW(0, 2);

    ddfidWdW(1, 2) = gam / rho;
    ddfidWdW(2, 1) = ddfidWdW(1, 2);

    ddFdWdW[2] = ddfidWdW;

    return ddFdWdW;
}

MatrixXd evalddlambdadWdW(std::vector <double> W, int i)
{
    double rho, u, p, c;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    p = (gam - 1) * ( W[i * 3 + 2] - rho * u * u / 2.0 );
    c = sqrt(p * gam / rho);

    MatrixXd ddlambdadWdWp(3, 3);

    ddlambdadWdWp(0, 0) = ( 6.0 * p * gam
                          + rho * u * (8.0 * c - u * gam * (gam - 1.0)) )
                          / (16.0 * c * pow(rho, 3.0));
    ddlambdadWdWp(1, 1) = -c * (gam - 1.0) / (4.0 * p);
    ddlambdadWdWp(2, 2) = c * (1.0 - gam) / (8.0 * p * p);

    ddlambdadWdWp(0, 1) = 0.25 * (c * u * (gam - 1.0) / p - 2.0 / rho);
    ddlambdadWdWp(0, 2) = - ( gam * gam * ( 2.0 * p + rho * u * u * (gam - 1.0) ) )
                          / ( 16.0 * pow(c * rho, 3) );

    ddlambdadWdWp(1, 0) = (-4.0 + u * gam * (gam - 1) / c)
                          / (8.0 * rho * rho);
    ddlambdadWdWp(1, 2) = c * u * (gam - 1.0) / (8.0 * p * p);

    ddlambdadWdWp(2, 0) = - gam * (gam - 1.0)
                          / (8.0 * c * rho * rho);
    ddlambdadWdWp(2, 1) = 0.0;

    return ddlambdadWdWp * dWpdW(W, i);
}

void evalddScalarFluxdWdW_FD(
    std::vector <MatrixXd> &ddFluxdWdW1,
    int ii, int jj,
    std::vector <double> W)
{
    int Fluxik;
    int Wi1, Wi2, Wi, Wj;
    double pertWi, pertWj;
    double Flux0, Flux1, Flux2, Flux3, Flux4;
    double h = 1e-3;
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Q(3 * nx, 0);

    MatrixXd dummy(3,3);
    dummy.setZero();
    for(int k = 0; k < 3 * (nx + 1); k++)
    {
        ddFluxdWdW1[k] = dummy;
    }
    for(int Fluxi = 1; Fluxi < nx; Fluxi++)
    {
        Wi1 = Fluxi + ii;
        Wi2 = Fluxi + jj;
        for(int Fluxk = 0; Fluxk < 3; Fluxk++)
        {
            Fluxik = Fluxi * 3 + Fluxk;
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                Wi = Wi1 * 3 + row;
                Wj = Wi2 * 3 + col;
                pertWi = W[Wi] * h;
                pertWj = W[Wj] * h;

                for(int m = 0; m < 3 * nx; m++)
                {
                    Wd[m] = W[m];
                }

                if(row != col || ii!=jj)
                {
                    // ddFluxdWdW1 = ddFlux(i)/dW(i-1)dW(i-1)
                    // R1
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    getFlux(Flux, Wd);

                    Flux1 = Flux[Fluxik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    getFlux(Flux, Wd);

                    Flux2 = Flux[Fluxik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    getFlux(Flux, Wd);

                    Flux3 = Flux[Fluxik];

                    // R4
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    getFlux(Flux, Wd);

                    Flux4 = Flux[Fluxik];

                    ddFluxdWdW1[Fluxik](row, col) =
                                (Flux1 - Flux2 - Flux3 + Flux4) / (4.0 * pertWi * pertWj);

                    // Symmetry
                    //ddFluxdWdW1[Fluxik](col, row) = ddFluxdWdW1[Fluxik](row, col);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }
                else
                {
                    getFlux(Flux, W);

                    Flux0 = Flux[Fluxik];

                    // R1
                    Wd[Wi] = W[Wi] + 2.0 * pertWi;

                    getFlux(Flux, Wd);

                    Flux1 = Flux[Fluxik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;

                    getFlux(Flux, Wd);

                    Flux2 = Flux[Fluxik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;

                    getFlux(Flux, Wd);

                    Flux3 = Flux[Fluxik];

                    // R4
                    Wd[Wi] = W[Wi] - 2.0 * pertWi;

                    getFlux(Flux, Wd);

                    Flux4 = Flux[Fluxik];

                    ddFluxdWdW1[Fluxik](row, col) =
                                    (-Flux1 + 16*Flux2 - 30*Flux0 + 16*Flux3 - Flux4)
                                                 / (12 * pertWi * pertWj);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }
            } // row col loop
        }// Fluxk
    }//Fluxi

}
void evalddScalarFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W)
{
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    std::vector <MatrixXd> ddFdWdW(3 * nx);
    std::vector <MatrixXd> ddFidWdW(3);

    std::vector <Vector3d> dlambdadW(nx);

    std::vector <MatrixXd> ddlambdadWdW(nx);
    for(int Wi = 0; Wi < nx; Wi++)
    {
        // Evaluate Convective Hessian
        ddFidWdW = evalddFdWdW(rho[Wi], u[Wi], p[Wi]);
        for(int k = 0; k < 3; k++)
        {
            ddFdWdW[Wi * 3 + k] = ddFidWdW[k];
        }
        // Evaluate dlambdadW
        dlambdadW[Wi] = Map<Vector3d>(evaldlambdadW(W, Wi).data());

        // Evaluate ddlambdadWdW
        ddlambdadWdW[Wi] = evalddlambdadWdW(W, Wi);
    }

    int Fluxkim, Fluxki;
    MatrixXd ddFluxidWdW(3, 3);
    for(int Fluxi = 1; Fluxi < nx; Fluxi++)
    {
        for(int Fluxk = 0; Fluxk < 3; Fluxk++)
        {
            Fluxkim = (Fluxi - 1) * 3 + Fluxk;
            Fluxki = Fluxi * 3 + Fluxk;
            // ddFluxdWdW1 = F(i+1/2) / W(i, i)
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxkim](row, col)
                    - 0.5 * Scalareps * ddlambdadWdW[Fluxi - 1](row, col)
                    * (W[Fluxki] - W[Fluxkim]);

                if(row == Fluxk)
                {
                    ddFluxidWdW(row, col) += 0.5 * Scalareps * dlambdadW[Fluxi - 1](col);
                }
                if(col == Fluxk)
                {
                    ddFluxidWdW(row, col) += 0.5 * Scalareps * dlambdadW[Fluxi - 1](row);
                }
            }
            ddFluxdWdW1[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW2 = F(i+1/2) / W(i + 1, i + 1)
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxki](row, col)
                    - 0.5 * Scalareps * ddlambdadWdW[Fluxi](row, col)
                    * (W[Fluxki] - W[Fluxkim]);
                if(row == Fluxk)
                {
                    ddFluxidWdW(row, col) -= 0.5 * Scalareps * dlambdadW[Fluxi](col);
                }
                if(col == Fluxk)
                {
                    ddFluxidWdW(row, col) -= 0.5 * Scalareps * dlambdadW[Fluxi](row);
                }
            }
            ddFluxdWdW2[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW3 = F(i+1/2) / W(i, i + 1)
            ddFluxidWdW.setZero();
            for(int row = 0; row < 3; row++)
            {
                ddFluxidWdW(row, Fluxk) -= 0.5 * Scalareps * dlambdadW[Fluxi - 1](row);
            }
            for(int col = 0; col < 3; col++)
            {
                // Important to have += since one of the entries is non-zero
                ddFluxidWdW(Fluxk, col) += 0.5 * Scalareps * dlambdadW[Fluxi](col);
            }
            ddFluxdWdW3[Fluxki] = ddFluxidWdW;
        }
    }

}


void evalddFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W)
{
    if(FluxScheme == 1)
    {
        evalddScalarFluxdWdW(ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3, W);
    }
}

std::vector <Matrix3d> ddWpdWdWp(
    std::vector <double> W,
    int i)
{
    double rho, u;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;

    std::vector <Matrix3d> M(3);
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
