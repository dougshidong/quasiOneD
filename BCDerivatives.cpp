#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"

using namespace Eigen;

std::vector <Matrix3d> ddWpdWdWp(
    std::vector <double> W,
    int i);

void HessianBC(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW,
    std::vector <MatrixXd> &ddRoutdWdW)
{
    // First Derivatives Required for Second Derivatives
    Matrix3d dRidWi, dRidWd, dRodWd, dRodWo;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        ddRindWdW[Rk].setZero();
        ddRoutdWdW[Rk].setZero();
    }

    std::vector <double> rho(nx), u(nx), e(nx), p(nx), c(nx), T(nx);
    WtoP(W, rho, u, e, p, c, T);

    // ************************
    // OUTLET HESSIAN
    // ************************

    double i1, i2;
    double r1, r2, p1, p2, u1, u2, c1, c2;
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

    // Shorthand
    double gamr, drho, dp, du, cr, uu;

    // ***********************************************************************
    // Speed of Sound
    double dc1dr1, dc2dr2, dc1dp1, dc2dp2;
    // First Derivative
    dc1dr1 = - p1 * gam / (2.0 * r1 * c1 * r1);
    dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
    dc1dp1 = gam / (2.0 * r1 * c1);
    dc2dp2 = gam / (2.0 * c2 * r2);

    // Second Derivative
    double ddc1dr1dr1, ddc1dp1dp1, ddc1dr1dp1, 
           ddc2dr2dr2, ddc2dp2dp2, ddc2dr2dp2; 
    ddc1dr1dr1 = 3.0 * c1 / (4.0 * r1 * r1);
    ddc1dp1dp1 = -c1 / (4.0 * p1 * p1);
    ddc1dr1dp1 = -gam / (4.0 * r1 * r1 * c1);

    ddc2dr2dr2 = 3.0 * c2 / (4.0 * r2 * r2);
    ddc2dp2dp2 = -c2 * (4.0 * p2 * p2);
    ddc2dr2dp2 = -gam / (4.0 * r2 * r2 * c2);

    // ***********************************************************************
    // Eigenvalue Derivatives
    double eig1, eig2, eig3;
    eig1 = (u1 + u2) / 2.0;
    eig2 = eig1 + (c1 + c2) / 2.0;
    eig3 = eig1 - (c1 + c2) / 2.0;

    // First Derivative
    double deig1du1, deig1du2;
    double deig2dr1, deig2du1, deig2dp1, deig2dr2, deig2du2, deig2dp2;
    double deig3dr1, deig3du1, deig3dp1, deig3dr2, deig3du2, deig3dp2;
    deig1du1 = 0.5;
    deig1du2 = 0.5;

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

    // Second Derivative
    double ddeig2dr1dr1, ddeig2dp1dp1, ddeig2dr1dp1;
    double ddeig2dr2dr2, ddeig2dp2dp2, ddeig2dr2dp2;

    double ddeig3dr1dr1, ddeig3dp1dp1, ddeig3dr1dp1;
    double ddeig3dr2dr2, ddeig3dp2dp2, ddeig3dr2dp2;

    // ddeig1 = 0

    // ddeig2
    ddeig2dr1dr1 = ddc1dr1dr1 / 2.0;
    ddeig2dp1dp1 = ddc1dp1dp1 / 2.0;
    ddeig2dr1dp1 = ddc1dr1dp1 / 2.0;

    ddeig2dr2dr2 = ddc2dr2dr2 / 2.0;
    ddeig2dp2dp2 = ddc2dp2dp2 / 2.0;
    ddeig2dr2dp2 = ddc2dr2dp2 / 2.0;

    // ddeig3
    ddeig3dr1dr1 = -ddc1dr1dr1 / 2.0;
    ddeig3dp1dp1 = -ddc1dp1dp1 / 2.0;
    ddeig3dr1dp1 = -ddc1dr1dp1 / 2.0;

    ddeig3dr2dr2 = -ddc2dr2dr2 / 2.0;
    ddeig3dp2dp2 = -ddc2dp2dp2 / 2.0;
    ddeig3dr2dp2 = -ddc2dr2dp2 / 2.0;

    // ***********************************************************************
    // Riemann invariants
    double R1, R2, R3;
    double dR1dr1, dR1du1, dR1dp1, dR1dr2, dR1du2, dR1dp2;
    double dR2dr1, dR2du1, dR2dp1, dR2dr2, dR2du2, dR2dp2;
    double dR3dr1, dR3du1, dR3dp1, dR3dr2, dR3du2, dR3dp2;
    R1 = - eig1 * ((r1 - r2) - (p1 - p2) / (c1 * c1));
    R2 = - eig2 * ((p1 - p2) + r1 * c1 * (u1 - u2));
    R3 = - eig3 * ((p1 - p2) - r1 * c1 * (u1 - u2));

    // First Derivative of Riemann Invariant
    dR1dr1 = -eig1 * (1.0 + 2.0 * (p1 - p2) * dc1dr1 / pow(c1, 3) );
    dR1du1 = deig1du1 * ((p1 - p2) - c1 * c1 * (r1 - r2)) / (c1 * c1);
    dR1dp1 = eig1 * (c1 - 2.0 * (p1 - p2) * dc1dp1) / pow(c1, 3);
    dR1dr2 = eig1;
    dR1du2 = deig1du2 * ((p1 - p2) - c1 * c1 * (r1 - r2)) / (c1 * c1);
    dR1dp2 = -eig1 / (c1 * c1);

    dR2dr1 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dr1
             - (u1 - u2) * eig2 * (c1 + r1 * dc1dr1); 
    dR2du1 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2du1
             - r1 * c1 * eig2;
    dR2dp1 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dp1
             - eig2 * (1.0 + (u1 - u2) * r1 * dc1dp1);
    dR2dr2 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dr2;
    dR2du2 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2du2
             + r1 * c1 * eig2;
    dR2dp2 = -((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dp2 + eig2;

    dR3dr1 = -((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dr1
             + eig3 * (u1 - u2) * (c1 + r1 * dc1dr1);
    dR3du1 = -((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3du1 + r1 * c1 * eig3;
    dR3dp1 = -(p1 - p2) * deig3dp1 + (u1 - u2) * r1 * eig3 * dc1dp1 - eig3;
    dR3dr2 = -((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dr2;
    dR3du2 = -((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3du2 - r1 * c1 * eig3;
    dR3dp2 = -((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dp2 + eig3;

    // Second Derivative of Riemann Invariant
    double ddR1dr1dr1, ddR1du1dr1, ddR1dp1dr1, ddR1dr2dr1, ddR1du2dr1, ddR1dp2dr1,
           ddR1dr1du1, ddR1du1du1, ddR1dp1du1, ddR1dr2du1, ddR1du2du1, ddR1dp2du1,
           ddR1dr1dp1, ddR1du1dp1, ddR1dp1dp1, ddR1dr2dp1, ddR1du2dp1, ddR1dp2dp1,
           ddR1dr1dr2, ddR1du1dr2, ddR1dp1dr2, ddR1dr2dr2, ddR1du2dr2, ddR1dp2dr2,
           ddR1dr1du2, ddR1du1du2, ddR1dp1du2, ddR1dr2du2, ddR1du2du2, ddR1dp2du2,
           ddR1dr1dp2, ddR1du1dp2, ddR1dp1dp2, ddR1dr2dp2, ddR1du2dp2, ddR1dp2dp2;
    double ddR2dr1dr1, ddR2du1dr1, ddR2dp1dr1, ddR2dr2dr1, ddR2du2dr1, ddR2dp2dr1,
           ddR2dr1du1, ddR2du1du1, ddR2dp1du1, ddR2dr2du1, ddR2du2du1, ddR2dp2du1,
           ddR2dr1dp1, ddR2du1dp1, ddR2dp1dp1, ddR2dr2dp1, ddR2du2dp1, ddR2dp2dp1,
           ddR2dr1dr2, ddR2du1dr2, ddR2dp1dr2, ddR2dr2dr2, ddR2du2dr2, ddR2dp2dr2,
           ddR2dr1du2, ddR2du1du2, ddR2dp1du2, ddR2dr2du2, ddR2du2du2, ddR2dp2du2,
           ddR2dr1dp2, ddR2du1dp2, ddR2dp1dp2, ddR2dr2dp2, ddR2du2dp2, ddR2dp2dp2;
    double ddR3dr1dr1, ddR3du1dr1, ddR3dp1dr1, ddR3dr2dr1, ddR3du2dr1, ddR3dp2dr1,
           ddR3dr1du1, ddR3du1du1, ddR3dp1du1, ddR3dr2du1, ddR3du2du1, ddR3dp2du1,
           ddR3dr1dp1, ddR3du1dp1, ddR3dp1dp1, ddR3dr2dp1, ddR3du2dp1, ddR3dp2dp1,
           ddR3dr1dr2, ddR3du1dr2, ddR3dp1dr2, ddR3dr2dr2, ddR3du2dr2, ddR3dp2dr2,
           ddR3dr1du2, ddR3du1du2, ddR3dp1du2, ddR3dr2du2, ddR3du2du2, ddR3dp2du2,
           ddR3dr1dp2, ddR3du1dp2, ddR3dp1dp2, ddR3dr2dp2, ddR3du2dp2, ddR3dp2dp2;

    // R1
    ddR1dr1dr1 = (-2.0 * eig1 * (p1 - p2) 
                 * (-3 * pow(dc1dr1, 2) 
                 + c1 * ddc1dr1dr1)) / pow(c1, 4);
    ddR1du1dr1 = -((1.0 + (2.0 * (p1 - p2) * dc1dr1) 
                 / pow(c1, 3)) * deig1du1);
    ddR1dp1dr1 = (-2.0 * eig1 * 
                 ((c1 + 3 * (-p1 + p2) * dc1dp1) 
                 * dc1dr1 + c1 * (p1 - p2) 
                 * ddc1dr1dp1)) / pow(c1, 4);
    ddR1dr2dr1 = 0.0;
    ddR1du2dr1 = -((1.0 + (2.0 * (p1 - p2) * dc1dr1)
                 / pow(c1, 3)) * deig1du2);
    ddR1dp2dr1 = (2.0 * eig1 * dc1dr1) / pow(c1, 3);

    ddR1dr1du1 = ddR1du1dr1;
    ddR1du1du1 = (p1 - p2) / pow(c1, 2);
    ddR1dp1du1 = ((c1 + 2.0 * (-p1 + p2) * dc1dp1) * deig1du1) / pow(c1, 3);
    ddR1dr2du1 = deig1du1;
    ddR1du2du1 = (p1 - p2) / pow(c1, 2);
    ddR1dp2du1 = -(deig1du1 / pow(c1, 2));

    ddR1dr1dp1 = ddR1dp1dr1;
    ddR1du1dp1 = ddR1dp1du1;
    ddR1dp1dp1 = (2.0 * eig1 * 
                 (-2.0 * c1 * dc1dp1 
                 + 3 * (p1 - p2) * pow(dc1dp1, 2) 
                 + c1 * (-p1 + p2) * ddc1dp1dp1))
                 / pow(c1, 4);
    ddR1dr2dp1 = 0.0;
    ddR1du2dp1 = ((c1 + 2.0 * (-p1 + p2) * dc1dp1) 
                 * deig1du2) / pow(c1, 3);
    ddR1dp2dp1 = (2.0 * eig1 * dc1dp1) / pow(c1, 3);

    ddR1dr1dr2 = 0.0;
    ddR1du1dr2 = deig1du1;
    ddR1dp1dr2 = 0.0;
    ddR1dr2dr2 = 0.0;
    ddR1du2dr2 = deig1du2;
    ddR1dp2dr2 = 0.0;

    ddR1dr1du2 = -((pow(c1, 3) + 2.0 * (p1 - p2) * dc1dr1) * deig1du2) 
                 / pow(c1, 3);
    ddR1du1du2 = (p1 - p2) / pow(c1, 2);
    ddR1dp1du2 = ((c1 + 2.0 * (-p1 + p2) * dc1dp1) 
                 * deig1du2) / pow(c1, 3);
    ddR1dr2du2 = ddR1du2dr2;
    ddR1du2du2 = 0.0;
    ddR1dp2du2 = -(deig1du2 / pow(c1, 2));

    ddR1dr1dp2 = (2.0 * eig1 * dc1dr1) / pow(c1, 3);
    ddR1du1dp2 = -(deig1du1 / pow(c1, 2));
    ddR1dp1dp2 = (2.0 * eig1 * dc1dp1) / pow(c1, 3);
    ddR1dr2dp2 = ddR1dp2dr2;
    ddR1du2dp2 = ddR1dp2du2;
    ddR1dp2dp2 = 0.0;

    // R2
    ddR2dr1dr1 = (-p1 + p2) * ddeig2dr1dr1 
                 - (u1 - u2) * (eig2 * r1 * ddc1dr1dr1 
                 + 2.0 * c1 * deig2dr1 
                 + 2.0 * dc1dr1 * (eig2 + r1 * deig2dr1) 
                 + c1 * r1 * ddeig2dr1dr1);
    ddR2du1dr1 = -(r1 * dc1dr1 * (eig2 + (u1 - u2) * deig2du1)) 
                 - c1 * (eig2 + r1 * deig2dr1 + (u1 - u2) * deig2du1);
    ddR2dp1dr1 = -deig2dr1 + (-p1 + p2) * ddeig2dr1dp1 
                 - (u1 - u2) * ((c1 + r1 * dc1dr1) * deig2dp1 
                 + dc1dp1 * (eig2 + r1 * deig2dr1) 
                 + r1 * (eig2 * ddc1dr1dp1 + c1 * ddeig2dr1dp1));
    ddR2dr2dr1 = -(u1 - u2) * (c1 + r1 * dc1dr1) * deig2dr2;
    ddR2du2dr1 = c1 * eig2 + eig2 * r1 * dc1dr1 
                 + c1 * r1 * deig2dr1 
                 - (u1 - u2) * (c1 + r1 * dc1dr1) * deig2du2;
    ddR2dp2dr1 = -((u1 - u2) * (c1 + r1 * dc1dr1) * deig2dp2) + deig2dr1;

    ddR2dr1du1 = ddR2du1dr1;
    ddR2du1du1 = -2.0 * c1 * r1 * deig2du1;
    ddR2dp1du1 = -deig2du1 - r1 * (c1 * deig2dp1 
                 + dc1dp1 * (eig2 + (u1 - u2) * deig2du1));
    ddR2dr2du1 = -(c1 * r1 * deig2dr2);
    ddR2du2du1 = c1 * r1 * deig2du1 
                 - c1 * r1 * deig2du2;
    ddR2dp2du1 = -(c1 * r1 * deig2dp2) + deig2du1;

    ddR2dr1dp1 = ddR2dp1dr1;
    ddR2du1dp1 = ddR2dp1du1;
    ddR2dp1dp1 = -2.0 * deig2dp1 + (-p1 + p2) * ddeig2dp1dp1 
                 -r1 * (u1 - u2) * (eig2 * ddc1dp1dp1 
                 + 2.0 * dc1dp1 * deig2dp1 + c1 * ddeig2dp1dp1);
    ddR2dr2dp1 = (-1 + r1 * (-u1 + u2) * dc1dp1) * deig2dr2;
    ddR2du2dp1 = c1 * r1 * deig2dp1 - deig2du2 
                 + r1 * dc1dp1 * (eig2 + (-u1 + u2) * deig2du2);
    ddR2dp2dp1 = deig2dp1 
                 + (-1 + r1 * (-u1 + u2) * dc1dp1) * deig2dp2;

    ddR2dr1dr2 = - (u1 - u2) * (c1 + r1 * dc1dr1) * deig2dr2;
    ddR2du1dr2 = ddR2dr2du1;
    ddR2dp1dr2 = ddR2dr2dp1;
    ddR2dr2dr2 = (-p1 + p2 + c1 * r1 * (-u1 + u2)) * ddeig2dr2dr2;
    ddR2du2dr2 = c1 * r1 * deig2dr2;
    ddR2dp2dr2 = deig2dr2 
                 + (-p1 + p2 + c1 * r1 * (-u1 + u2)) * ddeig2dr2dp2;

    ddR2dr1du2 = ddR2du2dr1;
    ddR2du1du2 = ddR2du2du1;
    ddR2dp1du2 = ddR2du2dp1;
    ddR2dr2du2 = ddR2du2dr2;
    ddR2du2du2 = 2.0 * c1 * r1 * deig2du2;
    ddR2dp2du2 = c1 * r1 * deig2dp2 + deig2du2;

    ddR2dr1dp2 = ddR2dp2dr1;
    ddR2du1dp2 = ddR2dp2du1;
    ddR2dp1dp2 = ddR2dp2dp1;
    ddR2dr2dp2 = ddR2dp2dr2;
    ddR2du2dp2 = ddR2dp2du2;
    ddR2dp2dp2 = 2.0 * deig2dp2 
                 + (-p1 + p2 + c1 * r1 * (-u1 + u2)) * ddeig2dp2dp2;

    // R3
    ddR3dr1dr1 = (-p1 + p2) * ddeig3dr1dr1 
                 + (u1 - u2) * (eig3 * r1 * ddc1dr1dr1 
                 + 2.0 * c1 * deig3dr1 
                 + 2.0 * dc1dr1 * (eig3 + r1 * deig3dr1) 
                 + c1 * r1 * ddeig3dr1dr1);
    ddR3du1dr1 = c1 * eig3 + eig3 * r1 * dc1dr1 
                 + c1 * r1 * deig3dr1 
                 + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du1; 
    ddR3dp1dr1 = -deig3dr1 + (-p1 + p2) * ddeig3dr1dp1 
                 + (u1 - u2) * ((c1 + r1 * dc1dr1) * deig3dp1 
                 + dc1dp1 * (eig3 + r1 * deig3dr1) 
                 + r1 * (eig3 * ddc1dr1dp1 + c1 * ddeig3dr1dp1));
    ddR3dr2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dr2;
    ddR3du2dr1 = -(c1 * eig3) - eig3 * r1 * dc1dr1
                 - c1 * r1 * deig3dr1 
                 + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du2;
    ddR3dp2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dp2 + deig3dr1;

    ddR3dr1du1 = ddR3du1dr1;
    ddR3du1du1 = 2.0 * c1 * r1 * deig3du1;
    ddR3dp1du1 = c1 * r1 * deig3dp1 - deig3du1 
                 + r1 * dc1dp1 * (eig3 + (u1 - u2) * deig3du1);
    ddR3dr2du1 = c1 * r1 * deig3dr2;
    ddR3du2du1 = -(c1 * r1 * deig3du1) + c1 * r1 * deig3du2;
    ddR3dp2du1 = c1 * r1 * deig3dp2 + deig3du1;

    ddR3dr1dp1 = ddR3dp1dr1;
    ddR3du1dp1 = ddR3dp1du1;
    ddR3dp1dp1 = -2.0 * deig3dp1 + (-p1 + p2) * ddeig3dp1dp1 
                 + r1 * (u1 - u2) * (eig3 * ddc1dp1dp1 
                 + 2.0 * dc1dp1 * deig3dp1 + c1 * ddeig3dp1dp1);
    ddR3dr2dp1 = (-1 + r1 * (u1 - u2) * dc1dp1) * deig3dr2;
    ddR3du2dp1 = -(c1 * r1 * deig3dp1) 
                 - deig3du2 - r1 * dc1dp1 * (eig3 + (-u1 + u2) * deig3du2);
    ddR3dp2dp1 = deig3dp1 + (-1 + r1 * (u1 - u2) * dc1dp1) * deig3dp2;

    ddR3dr1dr2 = ddR3dr2dr1;
    ddR3du1dr2 = ddR3dr2du1;
    ddR3dp1dr2 = ddR3dr2dp1;
    ddR3dr2dr2 = (-p1 + p2 + c1 * r1 * (u1 - u2)) * ddeig3dr2dr2;
    ddR3du2dr2 = -(c1 * r1 * deig3dr2);
    ddR3dp2dr2 = deig3dr2 
                 + (-p1 + p2 + c1 * r1 * (u1 - u2)) * ddeig3dr2dp2;

    ddR3dr1du2 = ddR3du2dr1;
    ddR3du1du2 = ddR3du2du1;
    ddR3dp1du2 = ddR3du2dp1;
    ddR3dr2du2 = ddR3du2dr2;
    ddR3du2du2 = -2.0 * c1 * r1 * deig3du2;
    ddR3dp2du2 = -(c1 * r1 * deig3dp2) + deig3du2;

    ddR3dr1dp2 = ddR3dp2dr1;
    ddR3du1dp2 = ddR3dp2du1;
    ddR3dp1dp2 = ddR3dp2dp1;
    ddR3dr2dp2 = ddR3dp2dr2;
    ddR3du2dp2 = ddR3dp2du2;
    ddR3dp2dp2 = 2.0 * deig3dp2 
                 + (-p1 + p2 + c1 * r1 * (u1 - u2)) * ddeig3dp2dp2;

    // ***********************************************************************
    // dp1/dt
    double dp1dt;
    double dp1dtdr1, dp1dtdu1, dp1dtdp1;
    double dp1dtdr2, dp1dtdu2, dp1dtdp2;
    double
    ddp1dtdr1dr1, ddp1dtdu1dr1, ddp1dtdp1dr1, ddp1dtdr2dr1, ddp1dtdu2dr1, ddp1dtdp2dr1,
    ddp1dtdr1du1, ddp1dtdu1du1, ddp1dtdp1du1, ddp1dtdr2du1, ddp1dtdu2du1, ddp1dtdp2du1,
    ddp1dtdr1dp1, ddp1dtdu1dp1, ddp1dtdp1dp1, ddp1dtdr2dp1, ddp1dtdu2dp1, ddp1dtdp2dp1,
    ddp1dtdr1dr2, ddp1dtdu1dr2, ddp1dtdp1dr2, ddp1dtdr2dr2, ddp1dtdu2dr2, ddp1dtdp2dr2,
    ddp1dtdr1du2, ddp1dtdu1du2, ddp1dtdp1du2, ddp1dtdr2du2, ddp1dtdu2du2, ddp1dtdp2du2,
    ddp1dtdr1dp2, ddp1dtdu1dp2, ddp1dtdp1dp2, ddp1dtdr2dp2, ddp1dtdu2dp2, ddp1dtdp2dp2;
    if(u1 < c1)
    {
        dp1dt = 0.0;

        dp1dtdr1 = 0.0;
        dp1dtdu1 = 0.0;
        dp1dtdp1 = 0.0;
        dp1dtdr2 = 0.0;
        dp1dtdu2 = 0.0;
        dp1dtdp2 = 0.0;

        ddp1dtdr1dr1 = 0.0;
        ddp1dtdu1dr1 = 0.0;
        ddp1dtdp1dr1 = 0.0;
        ddp1dtdr2dr1 = 0.0;
        ddp1dtdu2dr1 = 0.0;
        ddp1dtdp2dr1 = 0.0;

        ddp1dtdr1du1 = 0.0;
        ddp1dtdu1du1 = 0.0;
        ddp1dtdp1du1 = 0.0;
        ddp1dtdr2du1 = 0.0;
        ddp1dtdu2du1 = 0.0;
        ddp1dtdp2du1 = 0.0;

        ddp1dtdr1dp1 = 0.0;
        ddp1dtdu1dp1 = 0.0;
        ddp1dtdp1dp1 = 0.0;
        ddp1dtdr2dp1 = 0.0;
        ddp1dtdu2dp1 = 0.0;
        ddp1dtdp2dp1 = 0.0;

        ddp1dtdr1dr2 = 0.0;
        ddp1dtdu1dr2 = 0.0;
        ddp1dtdp1dr2 = 0.0;
        ddp1dtdr2dr2 = 0.0;
        ddp1dtdu2dr2 = 0.0;
        ddp1dtdp2dr2 = 0.0;

        ddp1dtdr1du2 = 0.0;
        ddp1dtdu1du2 = 0.0;
        ddp1dtdp1du2 = 0.0;
        ddp1dtdr2du2 = 0.0;
        ddp1dtdu2du2 = 0.0;
        ddp1dtdp2du2 = 0.0;

        ddp1dtdr1dp2 = 0.0;
        ddp1dtdu1dp2 = 0.0;
        ddp1dtdp1dp2 = 0.0;
        ddp1dtdr2dp2 = 0.0;
        ddp1dtdu2dp2 = 0.0;
        ddp1dtdp2dp2 = 0.0;
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

        ddp1dtdr1dr1 = (ddR2dr1dr1 + ddR3dr1dr1) / 2.0;
        ddp1dtdu1dr1 = (ddR2du1dr1 + ddR3du1dr1) / 2.0;
        ddp1dtdp1dr1 = (ddR2dp1dr1 + ddR3dp1dr1) / 2.0;
        ddp1dtdr2dr1 = (ddR2dr2dr1 + ddR3dr2dr1) / 2.0;
        ddp1dtdu2dr1 = (ddR2du2dr1 + ddR3du2dr1) / 2.0;
        ddp1dtdp2dr1 = (ddR2dp2dr1 + ddR3dp2dr1) / 2.0;

        ddp1dtdr1du1 = ddp1dtdu1dr1;
        ddp1dtdu1du1 = (ddR2du1du1 + ddR3du1du1) / 2.0;
        ddp1dtdp1du1 = (ddR2dp1du1 + ddR3dp1du1) / 2.0;
        ddp1dtdr2du1 = (ddR2dr2du1 + ddR3dr2du1) / 2.0;
        ddp1dtdu2du1 = (ddR2du2du1 + ddR3du2du1) / 2.0;
        ddp1dtdp2du1 = (ddR2dp2du1 + ddR3dp2du1) / 2.0;

        ddp1dtdr1dp1 = ddp1dtdp1dr1;
        ddp1dtdu1dp1 = ddp1dtdp1du1;
        ddp1dtdp1dp1 = (ddR2dp1dp1 + ddR3dp1dp1) / 2.0;
        ddp1dtdr2dp1 = (ddR2dr2dp1 + ddR3dr2dp1) / 2.0;
        ddp1dtdu2dp1 = (ddR2du2dp1 + ddR3du2dp1) / 2.0;
        ddp1dtdp2dp1 = (ddR2dp2dp1 + ddR3dp2dp1) / 2.0;

        ddp1dtdr1dr2 = ddp1dtdr2dr1;
        ddp1dtdu1dr2 = ddp1dtdr2du1;
        ddp1dtdp1dr2 = ddp1dtdr2dp1;
        ddp1dtdr2dr2 = (ddR2dr2dr2 + ddR3dr2dr2) / 2.0;
        ddp1dtdu2dr2 = (ddR2du2dr2 + ddR3du2dr2) / 2.0;
        ddp1dtdp2dr2 = (ddR2dp2dr2 + ddR3dp2dr2) / 2.0;

        ddp1dtdr1du2 = ddp1dtdu2dr1;
        ddp1dtdu1du2 = ddp1dtdu2du1;
        ddp1dtdp1du2 = ddp1dtdu2dp1;
        ddp1dtdr2du2 = ddp1dtdu2dr2;
        ddp1dtdu2du2 = (ddR2du2du2 + ddR3du2du2) / 2.0;
        ddp1dtdp2du2 = (ddR2dp2du2 + ddR3dp2du2) / 2.0;

        ddp1dtdr1dp2 = ddp1dtdp2dr1;
        ddp1dtdu1dp2 = ddp1dtdp2du1;
        ddp1dtdp1dp2 = ddp1dtdp2dp1;
        ddp1dtdr2dp2 = ddp1dtdp2dr2;
        ddp1dtdu2dp2 = ddp1dtdp2du2;
        ddp1dtdp2dp2 = (ddR2dp2dp2 + ddR3dp2dp2) / 2.0;
    }

    // ***********************************************************************
    // drho1/dt
    double dr1dt;
    double dr1dtdr1, dr1dtdu1, dr1dtdp1;
    double dr1dtdr2, dr1dtdu2, dr1dtdp2;
    dr1dt = R1 + dp1dt / (c1 * c1);

    // First Derivative
    dr1dtdr1 = dR1dr1 + dp1dtdr1 / (c1 * c1) - 2.0 * dp1dt * dc1dr1 / pow(c1, 3);
    dr1dtdu1 = dR1du1 + dp1dtdu1 / (c1 * c1);
    dr1dtdp1 = dR1dp1 + dp1dtdp1 / (c1 * c1) - 2.0 * dp1dt * dc1dp1 / pow(c1, 3);
    dr1dtdr2 = dR1dr2 + dp1dtdr2 / (c1 * c1);
    dr1dtdu2 = dR1du2 + dp1dtdu2 / (c1 * c1);
    dr1dtdp2 = dR1dp2 + dp1dtdp2 / (c1 * c1);
    // Second Derivative
    double
    ddr1dtdr1dr1, ddr1dtdu1dr1, ddr1dtdp1dr1, ddr1dtdr2dr1, ddr1dtdu2dr1, ddr1dtdp2dr1,
    ddr1dtdr1du1, ddr1dtdu1du1, ddr1dtdp1du1, ddr1dtdr2du1, ddr1dtdu2du1, ddr1dtdp2du1,
    ddr1dtdr1dp1, ddr1dtdu1dp1, ddr1dtdp1dp1, ddr1dtdr2dp1, ddr1dtdu2dp1, ddr1dtdp2dp1,
    ddr1dtdr1dr2, ddr1dtdu1dr2, ddr1dtdp1dr2, ddr1dtdr2dr2, ddr1dtdu2dr2, ddr1dtdp2dr2,
    ddr1dtdr1du2, ddr1dtdu1du2, ddr1dtdp1du2, ddr1dtdr2du2, ddr1dtdu2du2, ddr1dtdp2du2,
    ddr1dtdr1dp2, ddr1dtdu1dp2, ddr1dtdp1dp2, ddr1dtdr2dp2, ddr1dtdu2dp2, ddr1dtdp2dp2;

    ddr1dtdr1dr1 = (6.0 * dp1dt * pow(dc1dr1, 2) 
                   - 2.0 * c1 * dp1dt * ddc1dr1dr1 
                   - 4.0 * c1 * dc1dr1 * dp1dtdr1 
                   + pow(c1, 2) * ddp1dtdr1dr1) / pow(c1, 4) 
                   + ddR1dr1dr1;
    ddr1dtdu1dr1 = (-2.0 * dc1dr1 * dp1dtdu1 
                   + c1 * ddp1dtdr1du1) / pow(c1, 3) 
                   + ddR1dr1du1;
    ddr1dtdp1dr1 = (dc1dp1 * (6 * dp1dt * dc1dr1 
                   - 2.0 * c1 * dp1dtdr1) 
                   + c1 * (-2.0 * dc1dr1 * dp1dtdp1 
                   - 2.0 * dp1dt * ddc1dr1dp1 
                   + c1 * ddp1dtdp1dr1)) 
                   / pow(c1, 4) + ddR1dp1dr1;
    ddr1dtdr2dr1 = (-2.0 * dc1dr1 * dp1dtdr2 
                   + c1 * ddp1dtdr1dr2) / pow(c1, 3) 
                   + ddR1dr1dr2;
    ddr1dtdu2dr1 = (-2.0 * dc1dr1 * dp1dtdu2 
                   + c1 * ddp1dtdr1du2) / pow(c1, 3) 
                   + ddR1dr1du2;
    ddr1dtdp2dr1 = (-2.0 * dc1dr1 * dp1dtdp2 
                   + c1 * ddp1dtdp2dr1) / pow(c1, 3) 
                   + ddR1dp2dr1;

    ddr1dtdr1du1 = ddr1dtdu1dr1;
    ddr1dtdu1du1 = ddp1dtdu1du1 / pow(c1, 2) + ddR1du1du1;
    ddr1dtdp1du1 = (-2.0 * dc1dp1 * dp1dtdu1 
                   + c1 * ddp1dtdp1du1) / pow(c1, 3) 
                   + ddR1dp1du1;
    ddr1dtdr2du1 = ddp1dtdr2du1 / pow(c1, 2) + ddR1dr2du1;
    ddr1dtdu2du1 = ddp1dtdu1du2 / pow(c1, 2) + ddR1du1du2;
    ddr1dtdp2du1 = ddp1dtdp2du1 / pow(c1, 2) + ddR1dp2du1;

    ddr1dtdr1dp1 = ddr1dtdp1dr1;
    ddr1dtdu1dp1 = ddr1dtdp1du1;
    ddr1dtdp1dp1 = (6 * dp1dt * pow(dc1dp1, 2) 
                   - 2.0 * c1 * dp1dt * ddc1dp1dp1 
                   - 4 * c1 * dc1dp1 * dp1dtdp1 
                   + pow(c1, 2) * ddp1dtdp1dp1) 
                   / pow(c1, 4) + ddR1dp1dp1;
    ddr1dtdr2dp1 = (-2.0 * dc1dp1 * dp1dtdr2 
                   + c1 * ddp1dtdp1dr2) / pow(c1, 3) 
                   + ddR1dp1dr2;
    ddr1dtdu2dp1 = (-2.0 * dc1dp1 * dp1dtdu2 
                   + c1 * ddp1dtdp1du2) / pow(c1, 3) 
                   + ddR1dp1du2;
    ddr1dtdp2dp1 = (-2.0 * dc1dp1 * dp1dtdp2 
                   + c1 * ddp1dtdp1dp2) / pow(c1, 3) 
                   + ddR1dp1dp2;
                   
    ddr1dtdr1dr2 = ddr1dtdr2dr1;
    ddr1dtdu1dr2 = ddr1dtdr2du1;
    ddr1dtdp1dr2 = ddr1dtdr2dp1; 
    ddr1dtdr2dr2 = ddp1dtdr2dr2 / pow(c1, 2) + ddR1dr2dr2;
    ddr1dtdu2dr2 = ddp1dtdr2du2 / pow(c1, 2) + ddR1dr2du2;
    ddr1dtdp2dr2 = ddp1dtdp2dr2 / pow(c1, 2) + ddR1dp2dr2;

    ddr1dtdr1du2 = ddr1dtdu2dr1;
    ddr1dtdu1du2 = ddr1dtdu2du1;
    ddr1dtdp1du2 = ddr1dtdu2dp1;
    ddr1dtdr2du2 = ddr1dtdu2dr2;
    ddr1dtdu2du2 = ddp1dtdu2du2 / pow(c1, 2) + ddR1du2du2;
    ddr1dtdp2du2 = ddp1dtdp2du2 / pow(c1, 2) + ddR1dp2du2;


    ddr1dtdr1dp2 = ddr1dtdp2dr1;
    ddr1dtdu1dp2 = ddr1dtdp2du1;
    ddr1dtdp1dp2 = ddr1dtdp2dp1;
    ddr1dtdr2dp2 = ddr1dtdp2dr2;
    ddr1dtdu2dp2 = ddr1dtdp2du2;
    ddr1dtdp2dp2 = ddp1dtdp2dp2 / pow(c1, 2) + ddR1dp2dp2;

    // ***********************************************************************
    // du1/dt
    double du1dt;
    du1dt = (R2 - dp1dt) / (r1 * c1);

    // First Derivative
    double du1dtdr1, du1dtdu1, du1dtdp1;
    double du1dtdr2, du1dtdu2, du1dtdp2;
    du1dtdr1 = ( (dp1dt - R2) * r1 * dc1dr1
               + c1 * (dp1dt - R2 - r1 * dp1dtdr1 + r1 * dR2dr1) )
               / (r1 * c1 * r1 * c1);
    du1dtdu1 = (dR2du1 - dp1dtdu1) / (r1 * c1);
    du1dtdp1 = ( (dp1dt - R2) * dc1dp1 + c1 * (dR2dp1 - dp1dtdp1) ) / (r1 * c1 * c1);
    du1dtdr2 = (dR2dr2 - dp1dtdr2) / (r1 * c1);
    du1dtdu2 = (dR2du2 - dp1dtdu2) / (r1 * c1);
    du1dtdp2 = (dR2dp2 - dp1dtdp2) / (r1 * c1);

    // Second Derivative
    double
    ddu1dtdr1dr1, ddu1dtdu1dr1, ddu1dtdp1dr1, ddu1dtdr2dr1, ddu1dtdu2dr1, ddu1dtdp2dr1,
    ddu1dtdr1du1, ddu1dtdu1du1, ddu1dtdp1du1, ddu1dtdr2du1, ddu1dtdu2du1, ddu1dtdp2du1,
    ddu1dtdr1dp1, ddu1dtdu1dp1, ddu1dtdp1dp1, ddu1dtdr2dp1, ddu1dtdu2dp1, ddu1dtdp2dp1,
    ddu1dtdr1dr2, ddu1dtdu1dr2, ddu1dtdp1dr2, ddu1dtdr2dr2, ddu1dtdu2dr2, ddu1dtdp2dr2,
    ddu1dtdr1du2, ddu1dtdu1du2, ddu1dtdp1du2, ddu1dtdr2du2, ddu1dtdu2du2, ddu1dtdp2du2,
    ddu1dtdr1dp2, ddu1dtdu1dp2, ddu1dtdp1dp2, ddu1dtdr2dp2, ddu1dtdu2dp2, ddu1dtdp2dp2;

    ddu1dtdr1dr1 = (2.0 * c1 * c1 * (-dp1dt + R2) 
                   + r1 * (2.0 * r1 * (-dp1dt + R2) * dc1dr1 * dc1dr1 
                   + 2.0 * c1 * dc1dr1 * (-dp1dt 
                   + R2 + r1 * dp1dtdr1 - r1 * dR2dr1) 
                   + c1 * (r1 * (dp1dt - R2) * ddc1dr1dr1 
                   + c1 * (2.0 * dp1dtdr1 - r1 * ddp1dtdr1dr1 
                   - 2.0 * dR2dr1 + r1 * ddR2dr1dr1))))
                   / pow(c1 * r1, 3);

    ddu1dtdu1dr1 = ((c1 + r1 * dc1dr1) * dp1dtdu1 
                   - (c1 + r1 * dc1dr1) * dR2du1
                   + c1 * r1 * (-ddp1dtdr1du1 + ddR2dr1du1))
                   / pow(c1 * r1, 2.0);

    ddu1dtdp1dr1 = (dc1dp1 * (2.0 * r1 * (-dp1dt + R2) * dc1dr1 
                   + c1 * (-dp1dt + R2 + r1 * dp1dtdr1
                   - r1 * dR2dr1)) + c1 * ((c1 + r1 * dc1dr1) * dp1dtdp1
                   - (c1 + r1 * dc1dr1) * dR2dp1 
                   + r1 * ((dp1dt - R2) * ddc1dr1dp1 
                   + c1 * (-ddp1dtdp1dr1 + ddR2dp1dr1))))
                   / (pow(c1, 3) * pow(r1, 2));

    ddu1dtdr2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdr2 
                   - (c1 + r1 * dc1dr1) * dR2dr2 
                   + c1 * r1 * (-ddp1dtdr1dr2 + ddR2dr1dr2)) 
                   / pow(c1 * r1, 2);
//                 ((c1 + r1 * dc1dr1) * dp1dtdr2 
//                 - (c1 + r1 * dc1dr1) * dR2dr2 
//                 + c1 * r1 * (-ddp1dtdr1dr2 + ddR2dr1dr2))
//                 /(c1^2*r1^2)

    ddu1dtdu2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdu2 
                   - (c1 + r1 * dc1dr1) * dR2du2
                   + c1 * r1 * (-ddp1dtdr1du2 + ddR2dr1du2))
                   / pow(c1 * r1, 2);

    ddu1dtdp2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdp2
                   - (c1 + r1 * dc1dr1) * dR2dp2
                   + c1 * r1 * (-ddp1dtdp2dr1 + ddR2dp2dr1))
                   / pow(c1 * r1, 2);

    ddu1dtdr1du1 = ddu1dtdu1dr1;
    ddu1dtdu1du1 = (-ddp1dtdu1du1 + ddR2du1du1) / (c1 * r1);
    ddu1dtdp1du1 = (dc1dp1 * (dp1dtdu1 - dR2du1) 
                   + c1 * (-ddp1dtdp1du1 + ddR2dp1du1))
                   / (c1 * c1 * r1);
    ddu1dtdr2du1 = (-ddp1dtdr2du1 + ddR2dr2du1) / (c1 * r1);
    ddu1dtdu2du1 = (-ddp1dtdu1du2 + ddR2du1du2) / (c1 * r1);
    ddu1dtdp2du1 = (-ddp1dtdp2du1 + ddR2dp2du1) / (c1 * r1);

    ddu1dtdr1dp1 = ddu1dtdp1dr1;
    ddu1dtdu1dp1 = ddu1dtdp1du1;
    ddu1dtdp1dp1 = (2.0 * (-dp1dt + R2) * dc1dp1 * dc1dp1 
                   + 2.0 * c1 * dc1dp1 * (dp1dtdp1 - dR2dp1)
                   + c1 * ((dp1dt - R2) * ddc1dp1dp1 
                   + c1 * (-ddp1dtdp1dp1 + ddR2dp1dp1)))
                   / (pow(c1, 3) * r1);
    ddu1dtdr2dp1 = (dc1dp1 * (dp1dtdr2 - dR2dr2) 
                   + c1 * (-ddp1dtdp1dr2 + ddR2dp1dr2))
                   / (c1 * c1 * r1);
    ddu1dtdu2dp1 = (dc1dp1 * (dp1dtdu2 - dR2du2) 
                   + c1 * (-ddp1dtdp1du2 + ddR2dp1du2))
                   / (c1 * c1 * r1);
    ddu1dtdp2dp1 = (dc1dp1 * (dp1dtdp2 - dR2dp2)
                   + c1 * (-ddp1dtdp1dp2 + ddR2dp1dp2))
                   / (c1 * c1 * r1);

    ddu1dtdr1dr2 = ddu1dtdr2dr1;
    ddu1dtdu1dr2 = ddu1dtdr2du1;
    ddu1dtdp1dr2 = ddu1dtdr2dp1; 
    ddu1dtdr2dr2 = (-ddp1dtdr2dr2 + ddR2dr2dr2) / (c1 * r1);
    ddu1dtdu2dr2 = (-ddp1dtdr2du2 + ddR2dr2du2) / (c1 * r1);
    ddu1dtdp2dr2 = (-ddp1dtdp2dr2 + ddR2dp2dr2) / (c1 * r1);

    ddu1dtdr1du2 = ddu1dtdu2dr1;
    ddu1dtdu1du2 = ddu1dtdu2du1;
    ddu1dtdp1du2 = ddu1dtdu2dp1;
    ddu1dtdr2du2 = ddu1dtdu2dr2;
    ddu1dtdu2du2 = (-ddp1dtdu2du2 + ddR2du2du2) / (c1 * r1);
    ddu1dtdp2du2 = (-ddp1dtdp2du2 + ddR2dp2du2) / (c1 * r1);


    ddu1dtdr1dp2 = ddu1dtdp2dr1;
    ddu1dtdu1dp2 = ddu1dtdp2du1;
    ddu1dtdp1dp2 = ddu1dtdp2dp1;
    ddu1dtdr2dp2 = ddu1dtdp2dr2;
    ddu1dtdu2dp2 = ddu1dtdp2du2;
    ddu1dtdp2dp2 = (-ddp1dtdp2dp2 + ddR2dp2dp2) / (c1 * r1);

    // ***********************************************************************
    // d(ru)1/dt
//  double dru1dt;
//  dru1dt = r1 * du1dt + u1 * dr1dt;
    double dru1dtdr1, dru1dtdu1, dru1dtdp1;
    double dru1dtdr2, dru1dtdu2, dru1dtdp2;
    dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
    dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
    dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
    dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
    dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
    dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

    // Second Derivative
    double
    ddru1dtdr1dr1, ddru1dtdu1dr1, ddru1dtdp1dr1, ddru1dtdr2dr1, ddru1dtdu2dr1, ddru1dtdp2dr1,
    ddru1dtdr1du1, ddru1dtdu1du1, ddru1dtdp1du1, ddru1dtdr2du1, ddru1dtdu2du1, ddru1dtdp2du1,
    ddru1dtdr1dp1, ddru1dtdu1dp1, ddru1dtdp1dp1, ddru1dtdr2dp1, ddru1dtdu2dp1, ddru1dtdp2dp1,
    ddru1dtdr1dr2, ddru1dtdu1dr2, ddru1dtdp1dr2, ddru1dtdr2dr2, ddru1dtdu2dr2, ddru1dtdp2dr2,
    ddru1dtdr1du2, ddru1dtdu1du2, ddru1dtdp1du2, ddru1dtdr2du2, ddru1dtdu2du2, ddru1dtdp2du2,
    ddru1dtdr1dp2, ddru1dtdu1dp2, ddru1dtdp1dp2, ddru1dtdr2dp2, ddru1dtdu2dp2, ddru1dtdp2dp2;

    // Top Left
    ddru1dtdr2dr2 = u1 * ddr1dtdr2dr2 + r1 * ddu1dtdr2dr2;
    ddru1dtdu2du2 = u1 * ddr1dtdu2du2 + r1 * ddu1dtdu2du2;
    ddru1dtdp2dp2 = u1 * ddr1dtdp2dp2 + r1 * ddu1dtdp2dp2;

    ddru1dtdr2du2 = u1 * ddr1dtdr2du2 + r1 * ddu1dtdr2du2;
    ddru1dtdu2dr2 = u1 * ddr1dtdu2dr2 + r1 * ddu1dtdu2dr2;

    ddru1dtdr2dp2 = u1 * ddr1dtdr2dp2 + r1 * ddu1dtdr2dp2;
    ddru1dtdp2dr2 = u1 * ddr1dtdp2dr2 + r1 * ddu1dtdp2dr2;

    ddru1dtdu2dp2 = u1 * ddr1dtdu2dp2 + r1 * ddu1dtdu2dp2;
    ddru1dtdp2du2 = u1 * ddr1dtdp2du2 + r1 * ddu1dtdp2du2;

    // Bottom Right
    ddru1dtdr1dr1 = u1 * ddr1dtdr1dr1 + 2.0 * du1dtdr1 + r1 * ddu1dtdr1dr1;
    ddru1dtdu1du1 = u1 * ddr1dtdu1du1 + 2.0 * dr1dtdu1 + r1 * ddu1dtdu1du1;
    ddru1dtdp1dp1 = u1 * ddr1dtdp1dp1 + r1 * ddu1dtdp1dp1;

    ddru1dtdr1du1 = u1 * ddr1dtdr1du1 + dr1dtdr1 + du1dtdu1 + r1 * ddu1dtdr1du1;
    ddru1dtdu1dr1 = ddru1dtdr1du1;

    ddru1dtdr1dp1 = u1 * ddr1dtdp1dr1 + du1dtdp1 + r1 * ddu1dtdp1dr1;
    ddru1dtdp1dr1 = ddru1dtdr1dp1;
    
    ddru1dtdu1dp1 = u1 * ddr1dtdp1du1 + dr1dtdp1 + r1 * ddu1dtdp1du1;
    ddru1dtdp1du1 = ddru1dtdu1dp1;

    // Top Right

    ddru1dtdr2dr1 = u1 * ddr1dtdr1dr2 + du1dtdr2 + r1 * ddu1dtdr1dr2;
    ddru1dtdu2du1 = u1 * ddr1dtdu2du1 + dr1dtdu2 + r1 * ddu1dtdu2du1;
    ddru1dtdp2dp1 = u1 * ddr1dtdp2dp1 + r1 * ddu1dtdp2dp1;

    ddru1dtdr2du1 = u1 * ddr1dtdr2du1 + dr1dtdr2 + r1 * ddu1dtdr2du1;
    ddru1dtdu2dr1 = u1 * ddr1dtdr1du1 + du1dtdu2 + r1 * ddu1dtdr1du1;

    ddru1dtdr2dp1 = u1 * ddr1dtdr2dp1 + r1 * ddu1dtdr2dp1;
    ddru1dtdp2dr1 = u1 * ddr1dtdp2dr1 + du1dtdp2 + r1 * ddu1dtdp2dr1;

    ddru1dtdu2dp1 = u1 * ddr1dtdu2dp1 + r1 * ddu1dtdu2dp1;
    ddru1dtdp2du1 = u1 * ddr1dtdp2du1 + dr1dtdp2 + r1 * ddu1dtdp2du1;


    // Bottom Left

//  ddru1dtdr1dr2 = ddru1dtdr2dr1;
//  ddru1dtdu1du2 = ddru1dtdu2du1;
//  ddru1dtdp1dp2 = ddru1dtdp2dp1;

//  ddru1dtdr1du2 = ddru1dtdr1du2;
//  ddru1dtdu1dr2 = ddru1dtdu1dr2;

//  ddru1dtdp1du2 = ddru1dtdu2dp1;

//  ddru1dtdr1dp2 = ddru1dtdp2dr1;
//  ddru1dtdu1dp2 = ddru1dtdp2du1;
//  ddru1dtdp1dr2 = ddru1dtdr2dp1; 

    // ***********************************************************************
    // de1/dt
//  double de1dt;
//  de1dt = dp1dt * Cv / R + u1 * r1 * du1dt + (u1 * u1) * dr1dt / 2.0;
    double de1dtdr1, de1dtdu1, de1dtdp1;
    double de1dtdr2, de1dtdu2, de1dtdp2;

    de1dtdr1 = dp1dtdr1 * Cv / R + (u1 * u1) * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1
               + du1dt * u1;
    de1dtdu1 = dp1dtdu1 * Cv / R + (u1 * u1) * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1
               + du1dt * r1 + dr1dt * u1;
    de1dtdp1 = dp1dtdp1 / (gam - 1) + (u1 * u1) * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
    de1dtdr2 = dp1dtdr2 / (gam - 1) + (u1 * u1) * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
    de1dtdu2 = dp1dtdu2 / (gam - 1) + (u1 * u1) * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
    de1dtdp2 = dp1dtdp2 / (gam - 1) + (u1 * u1) * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

    // Second Derivative
    double
    dde1dtdr1dr1, dde1dtdu1dr1, dde1dtdp1dr1, dde1dtdr2dr1, dde1dtdu2dr1, dde1dtdp2dr1,
    dde1dtdr1du1, dde1dtdu1du1, dde1dtdp1du1, dde1dtdr2du1, dde1dtdu2du1, dde1dtdp2du1,
    dde1dtdr1dp1, dde1dtdu1dp1, dde1dtdp1dp1, dde1dtdr2dp1, dde1dtdu2dp1, dde1dtdp2dp1,
    dde1dtdr1dr2, dde1dtdu1dr2, dde1dtdp1dr2, dde1dtdr2dr2, dde1dtdu2dr2, dde1dtdp2dr2,
    dde1dtdr1du2, dde1dtdu1du2, dde1dtdp1du2, dde1dtdr2du2, dde1dtdu2du2, dde1dtdp2du2,
    dde1dtdr1dp2, dde1dtdu1dp2, dde1dtdp1dp2, dde1dtdr2dp2, dde1dtdu2dp2, dde1dtdp2dp2;

    dde1dtdr1dr1 = (Cv * ddp1dtdr1dr1) / R 
                   + (u1 * (u1 * ddr1dtdr1dr1 + 4.0 * du1dtdr1 
                   + 2.0 * r1 * ddu1dtdr1dr1)) / 2.0;
    dde1dtdu1dr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1 
                   + u1 * du1dtdu1 + (Cv * ddp1dtdr1du1) / R 
                   + (u1 * u1 * ddr1dtdr1du1) / 2.0 
                   + r1 * u1 * ddu1dtdr1du1;
    dde1dtdp1dr1 = u1 * du1dtdp1 
                   + (Cv * ddp1dtdp1dr1) / R 
                   + (u1 * u1 * ddr1dtdp1dr1) / 2.0 
                   + r1 * u1 * ddu1dtdp1dr1;
    dde1dtdr2dr1 = u1 * du1dtdr2 
                   + (Cv * ddp1dtdr1dr2) / R 
                   + (u1 * u1 * ddr1dtdr1dr2) / 2.0 
                   + r1 * u1 * ddu1dtdr1dr2;
    dde1dtdu2dr1 = u1 * du1dtdu2 
                   + (Cv * ddp1dtdr1du2) / R 
                   + (u1 * u1 * ddr1dtdr1du2) / 2.0 
                   + r1 * u1 * ddu1dtdr1du2;
    dde1dtdp2dr1 = u1 * du1dtdp2 
                   + (Cv * ddp1dtdp2dr1) / R 
                   + (u1 * u1 * ddr1dtdp2dr1) / 2.0 
                   + r1 * u1 * ddu1dtdp2dr1;

    dde1dtdr1du1 = dde1dtdu1dr1;
    dde1dtdu1du1 = dr1dt + (Cv * ddp1dtdu1du1) / R 
                   + 2.0 * u1 * dr1dtdu1 
                   + (u1 * u1 * ddr1dtdu1du1) / 2.0 
                   + 2.0 * r1 * du1dtdu1 
                   + r1 * u1 * ddu1dtdu1du1;
    dde1dtdp1du1 = u1 * dr1dtdp1 
                   + r1 * du1dtdp1 
                   + (Cv * ddp1dtdp1du1) / R 
                   + (u1 * u1 * ddr1dtdp1du1) / 2.0 
                   + r1 * u1 * ddu1dtdp1du1;
    dde1dtdr2du1 = u1 * dr1dtdr2 
                   + r1 * du1dtdr2 
                   + (Cv * ddp1dtdr2du1) / R 
                   + (u1 * u1 * ddr1dtdr2du1) / 2.0 
                   + r1 * u1 * ddu1dtdr2du1;
    dde1dtdu2du1 = u1 * dr1dtdu2 
                   + r1 * du1dtdu2 
                   + (Cv * ddp1dtdu1du2) / R 
                   + (u1 * u1 * ddr1dtdu1du2) / 2.0 
                   + r1 * u1 * ddu1dtdu1du2;
     dde1dtdp2du1 = u1 * dr1dtdp2 
                   + r1 * du1dtdp2 
                   + (Cv * ddp1dtdp2du1) / R 
                   + (u1 * u1 * ddr1dtdp2du1) / 2.0 
                   + r1 * u1 * ddu1dtdp2du1;

    dde1dtdr1dp1 = dde1dtdp1dr1;
    dde1dtdu1dp1 = dde1dtdp1du1;
    dde1dtdp1dp1 = (Cv * ddp1dtdp1dp1) / R 
                   + (u1 * u1 * ddr1dtdp1dp1) / 2.0 
                   + r1 * u1 * ddu1dtdp1dp1;
    dde1dtdr2dp1 = (Cv * ddp1dtdp1dr2) / R 
                   + (u1 * u1 * ddr1dtdp1dr2) / 2.0 
                   + r1 * u1 * ddu1dtdp1dr2;
    dde1dtdu2dp1 = (Cv * ddp1dtdp1du2) / R 
                   + (u1 * u1 * ddr1dtdp1du2) / 2.0 
                   + r1 * u1 * ddu1dtdp1du2;
    dde1dtdp2dp1 = (Cv * ddp1dtdp1dp2) / R 
                   + (u1 * u1 * ddr1dtdp1dp2) / 2.0 
                   + r1 * u1 * ddu1dtdp1dp2;

    dde1dtdr1dr2 = dde1dtdr2dr1;
    dde1dtdu1dr2 = dde1dtdr2du1;
    dde1dtdp1dr2 = dde1dtdr2dp1; 
    dde1dtdr2dr2 = (Cv * ddp1dtdr2dr2) / R
                   + (u1 * u1 * ddr1dtdr2dr2) / 2.0 
                   + r1 * u1 * ddu1dtdr2dr2;
    dde1dtdu2dr2 = (Cv * ddp1dtdr2du2) / R
                   + (u1 * u1 * ddr1dtdr2du2) / 2.0 
                   + r1 * u1 * ddu1dtdr2du2;
    dde1dtdp2dr2 = (Cv * ddp1dtdp2dr2) / R
                   + (u1 * u1 * ddr1dtdp2dr2) / 2.0 
                   + r1 * u1 * ddu1dtdp2dr2;

    dde1dtdr1du2 = dde1dtdu2dr1;
    dde1dtdu1du2 = dde1dtdu2du1;
    dde1dtdp1du2 = dde1dtdu2dp1;
    dde1dtdr2du2 = dde1dtdu2dr2;
    dde1dtdu2du2 = (Cv * ddp1dtdu2du2) / R 
                   + (u1 * u1 * ddr1dtdu2du2) / 2.0 
                   + r1 * u1 * ddu1dtdu2du2;
    dde1dtdp2du2 = (Cv * ddp1dtdp2du2) / R 
                   + (u1 * u1 * ddr1dtdp2du2) / 2.0 
                   + r1 * u1 * ddu1dtdp2du2;

    dde1dtdr1dp2 = dde1dtdp2dr1;
    dde1dtdu1dp2 = dde1dtdp2du1;
    dde1dtdp1dp2 = dde1dtdp2dp1;
    dde1dtdr2dp2 = dde1dtdp2dr2;
    dde1dtdu2dp2 = dde1dtdp2du2;
    dde1dtdp2dp2 = (Cv * ddp1dtdp2dp2) / R 
                   + (u1 * u1 * ddr1dtdp2dp2) / 2.0 
                   + r1 * u1 * ddu1dtdp2dp2;

// ***********************************************************************
    dRodWo(0,0) = dr1dtdr1;
    dRodWo(0,1) = dr1dtdu1;
    dRodWo(0,2) = dr1dtdp1;
    dRodWo(1,0) = dru1dtdr1;
    dRodWo(1,1) = dru1dtdu1;
    dRodWo(1,2) = dru1dtdp1;
    dRodWo(2,0) = de1dtdr1;
    dRodWo(2,1) = de1dtdu1;
    dRodWo(2,2) = de1dtdp1;

    dRodWd(0,0) = dr1dtdr2;
    dRodWd(0,1) = dr1dtdu2;
    dRodWd(0,2) = dr1dtdp2;
    dRodWd(1,0) = dru1dtdr2;
    dRodWd(1,1) = dru1dtdu2;
    dRodWd(1,2) = dru1dtdp2;
    dRodWd(2,0) = de1dtdr2;
    dRodWd(2,1) = de1dtdu2;
    dRodWd(2,2) = de1dtdp2;

//  Top Left
    ddRoutdWdW[0](0, 0) = ddr1dtdr2dr2;
    ddRoutdWdW[0](1, 1) = ddr1dtdu2du2;
    ddRoutdWdW[0](2, 2) = ddr1dtdp2dp2;

    ddRoutdWdW[0](0, 1) = ddr1dtdr2du2;
    ddRoutdWdW[0](1, 0) = ddr1dtdu2dr2;

    ddRoutdWdW[0](0, 2) = ddr1dtdr2dp2;
    ddRoutdWdW[0](2, 0) = ddr1dtdp2dr2;

    ddRoutdWdW[0](1, 2) = ddr1dtdu2dp2;
    ddRoutdWdW[0](2, 1) = ddr1dtdp2du2;

//  Bottom Right
    ddRoutdWdW[0](3, 3) = ddr1dtdr1dr1;
    ddRoutdWdW[0](4, 4) = ddr1dtdu1du1;
    ddRoutdWdW[0](5, 5) = ddr1dtdp1dp1;

    ddRoutdWdW[0](3, 4) = ddr1dtdr1du1;
    ddRoutdWdW[0](4, 3) = ddr1dtdu1dr1;
    
    ddRoutdWdW[0](3, 5) = ddr1dtdr1dp1;
    ddRoutdWdW[0](5, 3) = ddr1dtdp1dr1;

    ddRoutdWdW[0](4, 5) = ddr1dtdu1dp1;
    ddRoutdWdW[0](5, 4) = ddr1dtdp1du1;

//  Top Right
    ddRoutdWdW[0](0, 3) = ddr1dtdr2dr1;
    ddRoutdWdW[0](1, 4) = ddr1dtdu2du1;
    ddRoutdWdW[0](2, 5) = ddr1dtdp2dp1;

    ddRoutdWdW[0](0, 4) = ddr1dtdr2du1;
    ddRoutdWdW[0](1, 3) = ddr1dtdu2dr1;

    ddRoutdWdW[0](0, 5) = ddr1dtdr2dp1;
    ddRoutdWdW[0](2, 3) = ddr1dtdp2dr1;

    ddRoutdWdW[0](1, 5) = ddr1dtdu2dp1;
    ddRoutdWdW[0](2, 4) = ddr1dtdp2du1;


//  Bottom Left
    ddRoutdWdW[0](3, 0) = ddr1dtdr1dr2;
    ddRoutdWdW[0](4, 1) = ddr1dtdu1du2;
    ddRoutdWdW[0](5, 2) = ddr1dtdp1dp2;

    ddRoutdWdW[0](4, 0) = ddr1dtdu1dr2;
    ddRoutdWdW[0](3, 1) = ddr1dtdr1du2;

    ddRoutdWdW[0](5, 0) = ddr1dtdp1dr2;
    ddRoutdWdW[0](3, 2) = ddr1dtdr1dp2;

    ddRoutdWdW[0](5, 1) = ddr1dtdp1du2;
    ddRoutdWdW[0](4, 2) = ddr1dtdu1dp2;

    ddRoutdWdW[1](0, 0) = ddru1dtdr2dr2;
    ddRoutdWdW[1](0, 1) = ddru1dtdr2du2;
    ddRoutdWdW[1](0, 2) = ddru1dtdr2dp2;
    ddRoutdWdW[1](0, 3) = ddru1dtdr2dr1;
    ddRoutdWdW[1](0, 4) = ddru1dtdr2du1;
    ddRoutdWdW[1](0, 5) = ddru1dtdr2dp1;
    ddRoutdWdW[1](1, 0) = ddru1dtdu2dr2;
    ddRoutdWdW[1](1, 1) = ddru1dtdu2du2;
    ddRoutdWdW[1](1, 2) = ddru1dtdu2dp2;
    ddRoutdWdW[1](1, 3) = ddru1dtdu2dr1;
    ddRoutdWdW[1](1, 4) = ddru1dtdu2du1;
    ddRoutdWdW[1](1, 5) = ddru1dtdu2dp1;
    ddRoutdWdW[1](2, 0) = ddru1dtdp2dr2;
    ddRoutdWdW[1](2, 1) = ddru1dtdp2du2;
    ddRoutdWdW[1](2, 2) = ddru1dtdp2dp2;
    ddRoutdWdW[1](2, 3) = ddru1dtdp2dr1;
    ddRoutdWdW[1](2, 4) = ddru1dtdp2du1;
    ddRoutdWdW[1](2, 5) = ddru1dtdp2dp1;

    ddRoutdWdW[1](3, 0) = ddru1dtdr1dr2;
    ddRoutdWdW[1](3, 1) = ddru1dtdr1du2;
    ddRoutdWdW[1](3, 2) = ddru1dtdr1dp2;
    ddRoutdWdW[1](3, 3) = ddru1dtdr1dr1;
    ddRoutdWdW[1](3, 4) = ddru1dtdr1du1;
    ddRoutdWdW[1](3, 5) = ddru1dtdr1dp1;
    ddRoutdWdW[1](4, 0) = ddru1dtdu1dr2;
    ddRoutdWdW[1](4, 1) = ddru1dtdu1du2;
    ddRoutdWdW[1](4, 2) = ddru1dtdu1dp2;
    ddRoutdWdW[1](4, 3) = ddru1dtdu1dr1;
    ddRoutdWdW[1](4, 4) = ddru1dtdu1du1;
    ddRoutdWdW[1](4, 5) = ddru1dtdu1dp1;
    ddRoutdWdW[1](5, 0) = ddru1dtdp1dr2;
    ddRoutdWdW[1](5, 1) = ddru1dtdp1du2;
    ddRoutdWdW[1](5, 2) = ddru1dtdp1dp2;
    ddRoutdWdW[1](5, 3) = ddru1dtdp1dr1;
    ddRoutdWdW[1](5, 4) = ddru1dtdp1du1;
    ddRoutdWdW[1](5, 5) = ddru1dtdp1dp1;

    ddRoutdWdW[2](0, 0) = dde1dtdr2dr2;
    ddRoutdWdW[2](0, 1) = dde1dtdr2du2;
    ddRoutdWdW[2](0, 2) = dde1dtdr2dp2;
    ddRoutdWdW[2](0, 3) = dde1dtdr2dr1;
    ddRoutdWdW[2](0, 4) = dde1dtdr2du1;
    ddRoutdWdW[2](0, 5) = dde1dtdr2dp1;
    ddRoutdWdW[2](1, 0) = dde1dtdu2dr2;
    ddRoutdWdW[2](1, 1) = dde1dtdu2du2;
    ddRoutdWdW[2](1, 2) = dde1dtdu2dp2;
    ddRoutdWdW[2](1, 3) = dde1dtdu2dr1;
    ddRoutdWdW[2](1, 4) = dde1dtdu2du1;
    ddRoutdWdW[2](1, 5) = dde1dtdu2dp1;
    ddRoutdWdW[2](2, 0) = dde1dtdp2dr2;
    ddRoutdWdW[2](2, 1) = dde1dtdp2du2;
    ddRoutdWdW[2](2, 2) = dde1dtdp2dp2;
    ddRoutdWdW[2](2, 3) = dde1dtdp2dr1;
    ddRoutdWdW[2](2, 4) = dde1dtdp2du1;
    ddRoutdWdW[2](2, 5) = dde1dtdp2dp1;

    ddRoutdWdW[2](3, 0) = dde1dtdr1dr2;
    ddRoutdWdW[2](3, 1) = dde1dtdr1du2;
    ddRoutdWdW[2](3, 2) = dde1dtdr1dp2;
    ddRoutdWdW[2](3, 3) = dde1dtdr1dr1;
    ddRoutdWdW[2](3, 4) = dde1dtdr1du1;
    ddRoutdWdW[2](3, 5) = dde1dtdr1dp1;
    ddRoutdWdW[2](4, 0) = dde1dtdu1dr2;
    ddRoutdWdW[2](4, 1) = dde1dtdu1du2;
    ddRoutdWdW[2](4, 2) = dde1dtdu1dp2;
    ddRoutdWdW[2](4, 3) = dde1dtdu1dr1;
    ddRoutdWdW[2](4, 4) = dde1dtdu1du1;
    ddRoutdWdW[2](4, 5) = dde1dtdu1dp1;
    ddRoutdWdW[2](5, 0) = dde1dtdp1dr2;
    ddRoutdWdW[2](5, 1) = dde1dtdp1du2;
    ddRoutdWdW[2](5, 2) = dde1dtdp1dp2;
    ddRoutdWdW[2](5, 3) = dde1dtdp1dr1;
    ddRoutdWdW[2](5, 4) = dde1dtdp1du1;
    ddRoutdWdW[2](5, 5) = dde1dtdp1dp1;


//  std::cout<<"ddRoutdWpdWp[2]"<<std::endl;
//  std::cout<<ddRoutdWdW[2]<<std::endl;
    // Get Transformation Matrices
//  MatrixXd dwpdw = dWpdW(W, nx - 2);
//  std::vector <Matrix3d> ddwpdwdwp = ddWpdWdWp(W, nx - 2);

//  MatrixXd temp(3, 3);
//  for(int Ri = 0; Ri < 3; Ri++)
//  {
//      temp.setZero();
//      temp += dwpdw.transpose() * ddRoutdWdW[Ri].topLeftCorner(3, 3);
//      for(int Wpi = 0; Wpi < 3; Wpi++)
//      {
//          temp += dRodWd(Ri, Wpi) * ddwpdwdwp[Wpi];
//      }
//      ddRoutdWdW[Ri].topLeftCorner(3, 3) = temp * dwpdw;
//  }

//  // Get Transformation Matrices
//  dwpdw = dWpdW(W, nx - 1);
//  ddwpdwdwp = ddWpdWdWp(W, nx - 1);

//  for(int Ri = 0; Ri < 3; Ri++)
//  {
//      temp.setZero();
//      temp += dwpdw.transpose() * ddRoutdWdW[Ri].bottomRightCorner(3, 3);
//      for(int Wpi = 0; Wpi < 3; Wpi++)
//      {
//          temp += dRodWo(Ri, Wpi) * ddwpdwdwp[Wpi];
//      }
//      ddRoutdWdW[Ri].bottomRightCorner(3, 3) = temp * dwpdw;
//  }

//  // Get Transformation Matrices
//  dwpdw = dWpdW(W, nx - 1);
//  MatrixXd dwpdw2 = dWpdW(W, nx - 2);

//  for(int Ri = 0; Ri < 3; Ri++)
//  {
//      ddRoutdWdW[Ri].topRightCorner(3, 3) = 
//          dwpdw.transpose() * ddRoutdWdW[Ri].topRightCorner(3, 3) * dwpdw2;
//  }
//  for(int Ri = 0; Ri < 3; Ri++)
//  {
//      ddRoutdWdW[Ri].bottomLeftCorner(3, 3) = 
//          dwpdw2.transpose() * ddRoutdWdW[Ri].bottomLeftCorner(3, 3) * dwpdw;
//  }

    // Subsonic Inlet
    if(u[0] < c[0])
    {
        // Values at time-step N
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

        // Shorthand
        drho = r2 - r1;
        dp = p2 - p1;
        du = u2 - u1;
        cr = r1 * c1;
        uu = u1 * u1;
        gamr = (gam - 1.0) / (gam + 1.0);

        // Speed of Sound
        dc1dr1 = - p1 * gam / (2.0 * cr * r1);
        dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
        dc1dp1 = gam / (2.0 * cr);
        dc2dp2 = gam / (2.0 * c2 * r2);

        // Eigenvalue
        eig1 = (u1 + u2) / 2.0;
        eig3 = eig1 - (c1 + c2) / 2.0;

        deig1du1 = 0.5;
        deig1du2 = 0.5;

        deig3dr1 = - dc1dr1 / 2.0;
        deig3du1 = deig1du1;
        deig3dp1 = - dc1dp1 / 2.0;
        deig3dr2 = - dc2dr2 / 2.0;
        deig3du2 = deig1du2;
        deig3dp2 = - dc2dp2 / 2.0;

        // Riemann Invariants
        R3 = - eig3 * (dp - cr * du);

        dR3dr1 = -eig3 * (-c1 * du - du * r1 * dc1dr1) 
                 - (dp - cr * du) * deig3dr1;
        dR3du1 = -cr * eig3 - (dp - cr * du) * deig3du1;
        dR3dp1 = eig3 * (1.0 + du * r1 * dc1dp1) 
                 - (dp - cr * du) * deig3dp1;
        dR3dr2 = -(dp - cr * du) * deig3dr2;
        dR3du2 = cr * eig3 - (dp - cr * du) * deig3du2;
        dR3dp2 = -eig3 - (dp - cr * du) * deig3dp2;

        ddR3dr1dr1 = (p1 - p2) * ddeig3dr1dr1 
                     + (u1 - u2) * (eig3 * r1 * ddc1dr1dr1 
                     + 2.0 * c1 * deig3dr1 
                     + 2.0 * dc1dr1 * (eig3 + r1 * deig3dr1) 
                     + c1 * r1 * ddeig3dr1dr1);
        ddR3du1dr1 = c1 * eig3 + eig3 * r1 * dc1dr1 
                     + c1 * r1 * deig3dr1 
                     + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du1;
        ddR3dp1dr1 = deig3dr1 + (p1 - p2) * ddeig3dr1dp1 
                     + (u1 - u2) * ((c1 + r1 * dc1dr1) * deig3dp1 
                     + dc1dp1 * (eig3 + r1 * deig3dr1) 
                     + r1 * (eig3 * ddc1dr1dp1 + c1 * ddeig3dr1dp1));
        ddR3dr2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dr2;
        ddR3du2dr1 = -(c1 * eig3) - eig3 * r1 * dc1dr1 
                     - c1 * r1 * deig3dr1 
                     + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du2;
        ddR3dp2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dp2 - deig3dr1;

        ddR3dr1du1 = ddR3du1dr1;
        ddR3du1du1 = 2.0 * c1 * r1 * deig3du1;
        ddR3dp1du1 = c1 * r1 * deig3dp1 + deig3du1 
                     + r1 * dc1dp1 * (eig3 + (u1 - u2) * deig3du1);
        ddR3dr2du1 = c1 * r1 * deig3dr2;
        ddR3du2du1 = -(c1 * r1 * deig3du1) + c1 * r1 * deig3du2;
        ddR3dp2du1 = c1 * r1 * deig3dp2 - deig3du1;

        ddR3dr1dp1 = ddR3dp1dr1;
        ddR3du1dp1 = ddR3dp1du1;
        ddR3dp1dp1 = 2.0 * deig3dp1 + (p1 - p2) * ddeig3dp1dp1 
                     + r1 * (u1 - u2) * (eig3 * ddc1dp1dp1 
                     + 2.0 * dc1dp1 * deig3dp1 + c1 * ddeig3dp1dp1);
        ddR3dr2dp1 = (1.0 + r1 * (u1 - u2) * dc1dp1) * deig3dr2;
        ddR3du2dp1 = -(c1 * r1 * deig3dp1) + deig3du2 
                     - r1 * dc1dp1 * (eig3 + (-u1 + u2) * deig3du2);
        ddR3dp2dp1 = -deig3dp1 + (1.0 + r1 * (u1 - u2) * dc1dp1) * deig3dp2;

        ddR3dr1dr2 = ddR3dr2dr1;
        ddR3du1dr2 = ddR3dr2du1;
        ddR3dp1dr2 = ddR3dr2dp1;
        ddR3dr2dr2 = (p1 - p2 + c1 * r1 * (u1 - u2)) * ddeig3dr2dr2;
        ddR3du2dr2 = -(c1 * r1 * deig3dr2);
        ddR3dp2dr2 = -deig3dr2 + (p1 - p2 + c1 * r1 * (u1 - u2)) * ddeig3dr2dp2;

        ddR3dr1du2 = ddR3du2dr1;
        ddR3du1du2 = ddR3du2du1;
        ddR3dp1du2 = ddR3du2dp1;
        ddR3dr2du2 = ddR3du2dr2;
        ddR3du2du2 = -2.0 * c1 * r1 * deig3du2;
        ddR3dp2du2 = -(c1 * r1 * deig3dp2) - deig3du2;

        ddR3dr1dp2 = ddR3dp1dr1;
        ddR3du1dp2 = ddR3dp1du1;
        ddR3dp1dp2 = ddR3dp1dp1;
        ddR3dr2dp2 = ddR3dp1dr2;
        ddR3du2dp2 = ddR3dp1du1;
        ddR3dp2dp2 = -2.0 * deig3dp2 
                     + (p1 - p2 + c1 * r1 * (u1 - u2)) * ddeig3dp2dp2;


        // dp1
        double dp1du1, ddp1du1du1, dddp1du1du1du1;
        // Same Values
        dp1du1 = (ptin * u1 * pow(1.0 - (gamr * u1 * u1) / a2, gam/(-1.0 + gam)) 
                 * gam * (2.0 * gamr)) / ((-a2 + gamr * u1 * u1) * (-1.0 + gam));

        ddp1du1du1 = (2.0 * gamr * ptin * 
                      pow(1.0 - (gamr * pow(u1, 2)) / a2
                      , gam/(-1 + gam)) * gam * 
                      (a2 - a2 * gam + gamr * pow(u1, 2) * (1 + gam)))
                      /(pow(a2 - gamr * pow(u1,2), 2) * pow(-1 + gam, 2));
                     
        dddp1du1du1du1 = (4.0 * pow(gamr, 2) * ptin * u1 * 
                         pow(1 - (gamr * pow(u1, 2)) / a2
                         , gam / (-1 + gam)) * gam *
                         (-3.0 * a2 * (-1.0 + gam) 
                         + gamr * pow(u1, 2) * (1.0 + gam)))
                         /(pow(-a2 + gamr * pow(u1, 2), 3) * 
                         pow(-1.0 + gam, 3));
                         
        // du1
        du1dt = R3 / (dp1du1 - cr);
        du1dtdr1 = dR3dr1 / (dp1du1 - cr)
                   - R3 * (-c1 - r1 * dc1dr1) / pow((dp1du1 - cr), 2);
        du1dtdu1 = dR3du1 / (dp1du1 - cr)
                   - R3 * ddp1du1du1 / pow((dp1du1 - cr), 2);
        du1dtdp1 = dR3dp1 / (dp1du1 - cr)
                   + (R3 * r1 * dc1dp1) / pow((dp1du1 - cr), 2);
        du1dtdr2 = dR3dr2 / (dp1du1 - cr);
        du1dtdu2 = dR3du2 / (dp1du1 - cr);
        du1dtdp2 = dR3dp2 / (dp1du1 - cr);

        ddu1dtdr1dr1 = (2.0 * R3 * pow(c1 + r1 * dc1dr1, 2) 
                       + (dp1du1 - c1 * r1) * R3 * 
                       (2.0 * dc1dr1 + r1 * ddc1dr1dr1) 
                       + 2.0 * (dp1du1 - c1 * r1) * 
                       (c1 + r1 * dc1dr1) * dR3dr1 
                       + pow(dp1du1 - c1 * r1, 2) * ddR3dr1dr1) 
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdu1dr1 = (ddp1du1du1 * (-2.0 * R3 * (c1 + r1 * dc1dr1) 
                       + (-dp1du1 + c1 * r1) * dR3dr1) 
                       + (dp1du1 - c1 * r1) * 
                       ((c1 + r1 * dc1dr1) * dR3du1 
                       + (dp1du1 - c1 * r1) * ddR3dr1du1))
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdp1dr1 = (dc1dp1 * (R3 * 
                       (dp1du1 + c1 * r1 + 2.0 * pow(r1, 2) * dc1dr1) 
                       + r1 * (dp1du1 - c1 * r1) * dR3dr1) 
                       + (dp1du1 - c1 * r1) * 
                       ((c1 + r1 * dc1dr1) * dR3dp1 
                       + r1 * R3 * ddc1dr1dp1 
                       + (dp1du1 - c1 * r1) * ddR3dp1dr1)) 
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdr2dr1 = ((c1 + r1 * dc1dr1) * dR3dr2 
                       + (dp1du1 - c1 * r1) * ddR3dr1dr2) 
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdu2dr1 = ((c1 + r1 * dc1dr1) * dR3du2 
                       + (dp1du1 - c1 * r1) * ddR3dr1du2) 
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdp2dr1 = ((c1 + r1 * dc1dr1) * dR3dp2 
                       + (dp1du1 - c1 * r1) * ddR3dp2dr1) 
                       / pow(dp1du1 - c1 * r1, 2);

        ddu1dtdr1du1 = ddu1dtdu1dr1;
        ddu1dtdu1du1 = (2.0 * R3 * pow(ddp1du1du1, 2) 
                       - 2.0 * (dp1du1 - c1 * r1) * ddp1du1du1 * dR3du1 
                       + (dp1du1 - c1 * r1) * (-(R3 * dddp1du1du1du1) 
                       + (dp1du1 - c1 * r1) * ddR3du1du1))
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdp1du1 = (r1 * dc1dp1 * (-2.0 * R3 * ddp1du1du1 
                       + (dp1du1 - c1 * r1) * dR3du1) 
                       + (dp1du1 - c1 * r1) * (-(ddp1du1du1 * dR3dp1) 
                       + (dp1du1 - c1 * r1) * ddR3dp1du1)) 
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdr2du1 = (-(ddp1du1du1 * dR3dr2) 
                       + (dp1du1 - c1 * r1) * ddR3dr2du1) 
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdu2du1 = (-(ddp1du1du1 * dR3du2) 
                       + (dp1du1 - c1 * r1) * ddR3du1du2) 
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdp2du1 = (-(ddp1du1du1 * dR3dp2) 
                       + (dp1du1 - c1 * r1) * ddR3dp2du1) 
                       / pow(dp1du1 - c1 * r1, 2);

        ddu1dtdr1dp1 = ddu1dtdp1dr1;
        ddu1dtdu1dp1 = ddu1dtdp1du1;
        ddu1dtdp1dp1 = (r1 * R3 * (2.0 * r1 * pow(dc1dp1, 2) 
                       + (dp1du1 - c1 * r1) * ddc1dp1dp1) 
                       + 2.0 * r1 * (dp1du1 - c1 * r1) * dc1dp1 * dR3dp1
                       + pow(dp1du1 - c1 * r1, 2) * ddR3dp1dp1)
                       / pow(dp1du1 - c1 * r1, 3);
        ddu1dtdr2dp1 = (r1 * dc1dp1 * dR3dr2 
                       + (dp1du1 - c1 * r1) * ddR3dp1dr2)
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdu2dp1 = (r1 * dc1dp1 * dR3du2 
                       + (dp1du1 - c1 * r1) * ddR3dp1du2) 
                       / pow(dp1du1 - c1 * r1, 2);
        ddu1dtdp2dp1 = (r1 * dc1dp1 * dR3dp2 
                       + (dp1du1 - c1 * r1) * ddR3dp1dp2) 
                       / pow(dp1du1 - c1 * r1, 2);

        ddu1dtdr1dr2 = ddu1dtdr2dr1;
        ddu1dtdu1dr2 = ddu1dtdr2du1;
        ddu1dtdp1dr2 = ddu1dtdr2dp1; 
        ddu1dtdr2dr2 = ddR3dr2dr2 / (dp1du1 - c1 * r1);
        ddu1dtdu2dr2 = ddR3dr2du2 / (dp1du1 - c1 * r1);
        ddu1dtdp2dr2 = ddR3dp2dr2 / (dp1du1 - c1 * r1);

        ddu1dtdr1du2 = ddu1dtdu2dr1;
        ddu1dtdu1du2 = ddu1dtdu2du1;
        ddu1dtdp1du2 = ddu1dtdu2dp1;
        ddu1dtdr2du2 = ddu1dtdu2dr2;
        ddu1dtdu2du2 = ddR3du2du2 / (dp1du1 - c1 * r1);
        ddu1dtdp2du2 = ddR3dp2du2 / (dp1du1 - c1 * r1);


        ddu1dtdr1dp2 = ddu1dtdp2dr1;
        ddu1dtdu1dp2 = ddu1dtdp2du1;
        ddu1dtdp1dp2 = ddu1dtdp2dp1;
        ddu1dtdr2dp2 = ddu1dtdp2dr2;
        ddu1dtdu2dp2 = ddu1dtdp2du2;
        ddu1dtdp2dp2 = ddR3dp2dp2 / (dp1du1 - c1 * r1);

        // Primitive values at time-step n+1
        double unp1, pnp1, rnp1, tnp1, funp1;
        unp1 = u1 + du1dt;
        pnp1 = ptin * pow(1 - gamr * pow(unp1, 2) / a2, gam / (gam - 1.0));
        tnp1 = Ttin * ( 1 - gamr * unp1 * unp1 / a2 );
        rnp1 = pnp1 / (R * tnp1);
        funp1 = 1.0 - gamr * unp1 * unp1 / a2;
        double dpnp1dunp1;
        dpnp1dunp1 = -2.0 * gamr * ptin * unp1 * pow(funp1, 1.0 / (gam - 1.0)) * gam
                     / (a2 * (gam - 1.0));

        // dp1
        dp1dt  = pnp1 - p1;
        dp1dtdr1 = dpnp1dunp1 * du1dtdr1;
        dp1dtdu1 = dpnp1dunp1 * (du1dtdu1 + 1.0);
        dp1dtdp1 = dpnp1dunp1 * du1dtdp1 - 1.0;
        dp1dtdr2 = dpnp1dunp1 * du1dtdr2;
        dp1dtdu2 = dpnp1dunp1 * du1dtdu2;
        dp1dtdp2 = dpnp1dunp1 * du1dtdp2;

        ddp1dtdr1dr1 = dpnp1dunp1 * ddu1dtdr1dr1;
        ddp1dtdu1dr1 = dpnp1dunp1 * ddu1dtdu1dr1;
        ddp1dtdp1dr1 = dpnp1dunp1 * ddu1dtdp1dr1;
        ddp1dtdr2dr1 = dpnp1dunp1 * ddu1dtdr2dr1;
        ddp1dtdu2dr1 = dpnp1dunp1 * ddu1dtdu2dr1;
        ddp1dtdp2dr1 = dpnp1dunp1 * ddu1dtdp2dr1;

        ddp1dtdr1du1 = ddp1dtdu1dr1;
        ddp1dtdu1du1 = dpnp1dunp1 * ddu1dtdu1du1;
        ddp1dtdp1du1 = dpnp1dunp1 * ddu1dtdp1du1;
        ddp1dtdr2du1 = dpnp1dunp1 * ddu1dtdr2du1;
        ddp1dtdu2du1 = dpnp1dunp1 * ddu1dtdu2du1;
        ddp1dtdp2du1 = dpnp1dunp1 * ddu1dtdp2du1;

        ddp1dtdr1dp1 = ddp1dtdp1dr1;
        ddp1dtdu1dp1 = ddp1dtdp1du1;
        ddp1dtdp1dp1 = dpnp1dunp1 * ddu1dtdp1dp1;
        ddp1dtdr2dp1 = dpnp1dunp1 * ddu1dtdr2dp1;
        ddp1dtdu2dp1 = dpnp1dunp1 * ddu1dtdu2dp1;
        ddp1dtdp2dp1 = dpnp1dunp1 * ddu1dtdp2dp1;

        ddp1dtdr1dr2 = ddp1dtdr2dr1;
        ddp1dtdu1dr2 = ddp1dtdr2du1;
        ddp1dtdp1dr2 = ddp1dtdr2dp1;
        ddp1dtdr2dr2 = dpnp1dunp1 * ddu1dtdr2dr2;
        ddp1dtdu2dr2 = dpnp1dunp1 * ddu1dtdu2dr2;
        ddp1dtdp2dr2 = dpnp1dunp1 * ddu1dtdp2dr2;

        ddp1dtdr1du2 = ddp1dtdu2dr1;
        ddp1dtdu1du2 = ddp1dtdu2du1;
        ddp1dtdp1du2 = ddp1dtdu2dp1;
        ddp1dtdr2du2 = ddp1dtdu2dr2;
        ddp1dtdu2du2 = dpnp1dunp1 * ddu1dtdu2du2;
        ddp1dtdp2du2 = dpnp1dunp1 * ddu1dtdp2du2;

        ddp1dtdr1dp2 = ddp1dtdp2dr1;
        ddp1dtdu1dp2 = ddp1dtdp2du1;
        ddp1dtdp1dp2 = ddp1dtdp2dp1;
        ddp1dtdr2dp2 = ddp1dtdp2dr2;
        ddp1dtdu2dp2 = ddp1dtdp2du2;
        ddp1dtdp2dp2 = dpnp1dunp1 * ddu1dtdp2dp2;

        // dr1
        // Total derivative from rho_n+1 to p_n+1 and u_n+1
        double drnp1dpnp1, drnp1dtnp1, dtnp1dpnp1;
        drnp1dpnp1 = 1 / (R * tnp1);
        drnp1dtnp1 = -pnp1 / (R * tnp1 * tnp1);
        dtnp1dpnp1 = Ttin / ptin * (gam - 1.0) / gam * pow(pnp1 / ptin, - 1.0 / gam);
        double Drnp1Dpnp1 = drnp1dpnp1 + drnp1dtnp1 * dtnp1dpnp1;
        double drnp1dunp1 = Drnp1Dpnp1 * dpnp1dunp1;

        dr1dt = rnp1 - r1;

        dr1dtdr1 = drnp1dunp1 * du1dtdr1 - 1;
        dr1dtdu1 = drnp1dunp1 * (du1dtdu1 + 1);
        dr1dtdp1 = drnp1dunp1 * du1dtdp1;

        dr1dtdr2 = drnp1dunp1 * du1dtdr2;
        dr1dtdu2 = drnp1dunp1 * du1dtdu2;
        dr1dtdp2 = drnp1dunp1 * du1dtdp2;

        ddr1dtdr1dr1 = drnp1dunp1 * ddu1dtdr1dr1;
        ddr1dtdu1dr1 = drnp1dunp1 * ddu1dtdu1dr1;
        ddr1dtdp1dr1 = drnp1dunp1 * ddu1dtdp1dr1;
        ddr1dtdr2dr1 = drnp1dunp1 * ddu1dtdr2dr1;
        ddr1dtdu2dr1 = drnp1dunp1 * ddu1dtdu2dr1;
        ddr1dtdp2dr1 = drnp1dunp1 * ddu1dtdp2dr1;

        ddr1dtdr1du1 = ddp1dtdu1dr1;
        ddr1dtdu1du1 = drnp1dunp1 * ddu1dtdu1du1;
        ddr1dtdp1du1 = drnp1dunp1 * ddu1dtdp1du1;
        ddr1dtdr2du1 = drnp1dunp1 * ddu1dtdr2du1;
        ddr1dtdu2du1 = drnp1dunp1 * ddu1dtdu2du1;
        ddr1dtdp2du1 = drnp1dunp1 * ddu1dtdp2du1;

        ddr1dtdr1dp1 = ddp1dtdp1dr1;
        ddr1dtdu1dp1 = ddp1dtdp1du1;
        ddr1dtdp1dp1 = drnp1dunp1 * ddu1dtdp1dp1;
        ddr1dtdr2dp1 = drnp1dunp1 * ddu1dtdr2dp1;
        ddr1dtdu2dp1 = drnp1dunp1 * ddu1dtdu2dp1;
        ddr1dtdp2dp1 = drnp1dunp1 * ddu1dtdp2dp1;

        ddr1dtdr1dr2 = ddp1dtdr2dr1;
        ddr1dtdu1dr2 = ddp1dtdr2du1;
        ddr1dtdp1dr2 = ddp1dtdr2dp1;
        ddr1dtdr2dr2 = drnp1dunp1 * ddu1dtdr2dr2;
        ddr1dtdu2dr2 = drnp1dunp1 * ddu1dtdu2dr2;
        ddr1dtdp2dr2 = drnp1dunp1 * ddu1dtdp2dr2;

        ddr1dtdr1du2 = ddp1dtdu2dr1;
        ddr1dtdu1du2 = ddp1dtdu2du1;
        ddr1dtdp1du2 = ddp1dtdu2dp1;
        ddr1dtdr2du2 = ddp1dtdu2dr2;
        ddr1dtdu2du2 = drnp1dunp1 * ddu1dtdu2du2;
        ddr1dtdp2du2 = drnp1dunp1 * ddu1dtdp2du2;

        ddr1dtdr1dp2 = ddp1dtdp2dr1;
        ddr1dtdu1dp2 = ddp1dtdp2du1;
        ddr1dtdp1dp2 = ddp1dtdp2dp1;
        ddr1dtdr2dp2 = ddp1dtdp2dr2;
        ddr1dtdu2dp2 = ddp1dtdp2du2;
        ddr1dtdp2dp2 = drnp1dunp1 * ddu1dtdp2dp2;

        // dru1/dt
//      dru1dt = r1 * du1dt + u1 * dr1dt;

        dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
        dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
        dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
        dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
        dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
        dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

        ddru1dtdr1dr1 = u1 * ddr1dtdr1dr1 + 2.0 * du1dtdr1 + r1 * ddu1dtdr1dr1;
        ddru1dtdu1dr1 = u1 * ddr1dtdr1du1 + du1dtdr1 + du1dtdu1 + r1 * ddu1dtdr1du1;
        ddru1dtdp1dr1 = u1 * ddr1dtdp1dr1 + du1dtdp1 + r1 * ddu1dtdp1dr1;
        ddru1dtdr2dr1 = u1 * ddr1dtdr1dr2 + du1dtdr2 + r1 * ddu1dtdr1dr2;
        ddru1dtdu2dr1 = u1 * ddr1dtdr1du1 + du1dtdu2 + r1 * ddu1dtdr1du1;
        ddru1dtdp2dr1 = u1 * ddr1dtdp2dr1 + du1dtdp2 + r1 * ddu1dtdp2dr1;

        ddru1dtdr1du1 = ddru1dtdu1dr1;
        ddru1dtdu1du1 = u1 * ddr1dtdu1du1 + 2.0 * dr1dtdu1 + r1 * ddu1dtdu1du1;
        ddru1dtdp1du1 = u1 * ddr1dtdp1du1 + dr1dtdp1 + r1 * ddu1dtdp1du1;
        ddru1dtdr2du1 = u1 * ddr1dtdr2du1 + dr1dtdr2 + r1 * ddu1dtdr2du1;
        ddru1dtdu2du1 = u1 * ddr1dtdu2du1 + dr1dtdu2 + r1 * ddu1dtdu2du1;
        ddru1dtdp2du1 = u1 * ddr1dtdp2du1 + dr1dtdp2 + r1 * ddu1dtdp2du1;


        ddru1dtdr1dp1 = ddru1dtdp1dr1;
        ddru1dtdu1dp1 = ddru1dtdp1du1;
        ddru1dtdp1dp1 = u1 * ddr1dtdp1dp1 + r1 * ddu1dtdp1dp1;
        ddru1dtdr2dp1 = u1 * ddr1dtdr2dp1 + r1 * ddu1dtdr2dp1;
        ddru1dtdu2dp1 = u1 * ddr1dtdu2dp1 + r1 * ddu1dtdu2dp1;
        ddru1dtdp2dp1 = u1 * ddr1dtdp2dp1 + r1 * ddu1dtdp2dp1;

        ddru1dtdr1dr2 = ddru1dtdr2dr1;
        ddru1dtdu1dr2 = ddru1dtdr2du1;
        ddru1dtdp1dr2 = ddru1dtdr2dp1; 
        ddru1dtdr2dr2 = u1 * ddr1dtdr2dr2 + r1 * ddu1dtdr2dr2;
        ddru1dtdu2dr2 = u1 * ddr1dtdu2dr2 + r1 * ddu1dtdu2dr2;
        ddru1dtdp2dr2 = u1 * ddr1dtdp2dr2 + r1 * ddu1dtdp2dr2;

        ddru1dtdr1du2 = ddru1dtdu2dr1;
        ddru1dtdu1du2 = ddru1dtdu2du1;
        ddru1dtdp1du2 = ddru1dtdu2dp1;
        ddru1dtdr2du2 = ddru1dtdu2dr2;
        ddru1dtdu2du2 = u1 * ddr1dtdr2du2 + r1 * ddu1dtdr2du2;
        ddru1dtdp2du2 = u1 * ddr1dtdu2du2 + r1 * ddu1dtdu2du2;

        ddru1dtdr1dp2 = ddru1dtdp2dr1;
        ddru1dtdu1dp2 = ddru1dtdp2du1;
        ddru1dtdp1dp2 = ddru1dtdp2dp1;
        ddru1dtdr2dp2 = ddru1dtdp2dr2;
        ddru1dtdu2dp2 = ddru1dtdp2du2;
        ddru1dtdp2dp2 = u1 * ddr1dtdp2dp2 + r1 * ddu1dtdp2dp2;
        // de1/dt
//      de1dt = dp1dt / (gam - 1.0) + r1 * u1 * du1dt + uu * dr1dt / 2.0;

        de1dtdr1 = dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1
                   + du1dt * u1;
        de1dtdu1 = dp1dtdu1 * Cv / R + uu * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1
                   + du1dt * r1 + dr1dt * u1;
        de1dtdp1 = dp1dtdp1 / (gam - 1.0) + uu * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
        de1dtdr2 = dp1dtdr2 / (gam - 1.0) + uu * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
        de1dtdu2 = dp1dtdu2 / (gam - 1.0) + uu * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
        de1dtdp2 = dp1dtdp2 / (gam - 1.0) + uu * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

        dde1dtdr1dr1 = (Cv * ddp1dtdr1dr1) / R 
                       + (u1 * (u1 * ddr1dtdr1dr1 + 4.0 * du1dtdr1 
                       + 2.0 * r1 * ddu1dtdr1dr1)) / 2.0;
        dde1dtdu1dr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1 
                       + u1 * du1dtdu1 + (Cv * ddp1dtdr1du1) / R 
                       + (u1 * u1 * ddr1dtdr1du1) / 2.0 
                       + r1 * u1 * ddu1dtdr1du1;
        dde1dtdp1dr1 = u1 * du1dtdp1 
                       + (Cv * ddp1dtdp1dr1) / R 
                       + (u1 * u1 * ddr1dtdp1dr1) / 2.0 
                       + r1 * u1 * ddu1dtdp1dr1;
        dde1dtdr2dr1 = u1 * du1dtdr2 
                       + (Cv * ddp1dtdr1dr2) / R 
                       + (u1 * u1 * ddr1dtdr1dr2) / 2.0 
                       + r1 * u1 * ddu1dtdr1dr2;
        dde1dtdu2dr1 = u1 * du1dtdu2 
                       + (Cv * ddp1dtdr1du2) / R 
                       + (u1 * u1 * ddr1dtdr1du2) / 2.0 
                       + r1 * u1 * ddu1dtdr1du2;
        dde1dtdp2dr1 = u1 * du1dtdp2 
                       + (Cv * ddp1dtdp2dr1) / R 
                       + (u1 * u1 * ddr1dtdp2dr1) / 2.0 
                       + r1 * u1 * ddu1dtdp2dr1;

        dde1dtdr1du1 = dde1dtdu1dr1;
        dde1dtdu1du1 = dr1dt + (Cv * ddp1dtdu1du1) / R 
                       + 2.0 * u1 * dr1dtdu1 
                       + (u1 * u1 * ddr1dtdu1du1) / 2.0 
                       + 2.0 * r1 * du1dtdu1 
                       + r1 * u1 * ddu1dtdu1du1;
        dde1dtdp1du1 = u1 * dr1dtdp1 
                       + r1 * du1dtdp1 
                       + (Cv * ddp1dtdp1du1) / R 
                       + (u1 * u1 * ddr1dtdp1du1) / 2.0 
                       + r1 * u1 * ddu1dtdp1du1;
        dde1dtdr2du1 = u1 * dr1dtdr2 
                       + r1 * du1dtdr2 
                       + (Cv * ddp1dtdr2du1) / R 
                       + (u1 * u1 * ddr1dtdr2du1) / 2.0 
                       + r1 * u1 * ddu1dtdr2du1;
        dde1dtdu2du1 = u1 * dr1dtdu2 
                       + r1 * du1dtdu2 
                       + (Cv * ddp1dtdu1du2) / R 
                       + (u1 * u1 * ddr1dtdu1du2) / 2.0 
                       + r1 * u1 * ddu1dtdu1du2;
         dde1dtdp2du1 = u1 * dr1dtdp2 
                       + r1 * du1dtdp2 
                       + (Cv * ddp1dtdp2du1) / R 
                       + (u1 * u1 * ddr1dtdp2du1) / 2.0 
                       + r1 * u1 * ddu1dtdp2du1;

        dde1dtdr1dp1 = dde1dtdp1dr1;
        dde1dtdu1dp1 = dde1dtdp1du1;
        dde1dtdp1dp1 = (Cv * ddp1dtdp1dp1) / R 
                       + (u1 * u1 * ddr1dtdp1dp1) / 2.0 
                       + r1 * u1 * ddu1dtdp1dp1;
        dde1dtdr2dp1 = (Cv * ddp1dtdp1dr2) / R 
                       + (u1 * u1 * ddr1dtdp1dr2) / 2.0 
                       + r1 * u1 * ddu1dtdp1dr2;
        dde1dtdu2dp1 = (Cv * ddp1dtdp1du2) / R 
                       + (u1 * u1 * ddr1dtdp1du2) / 2.0 
                       + r1 * u1 * ddu1dtdp1du2;
        dde1dtdp2dp1 = (Cv * ddp1dtdp1dp2) / R 
                       + (u1 * u1 * ddr1dtdp1dp2) / 2.0 
                       + r1 * u1 * ddu1dtdp1dp2;

        dde1dtdr1dr2 = dde1dtdr2dr1;
        dde1dtdu1dr2 = dde1dtdr2du1;
        dde1dtdp1dr2 = dde1dtdr2dp1; 
        dde1dtdr2dr2 = (Cv * ddp1dtdr2dr2) / R
                       + (u1 * u1 * ddr1dtdr2dr2) / 2.0 
                       + r1 * u1 * ddu1dtdr2dr2;
        dde1dtdu2dr2 = (Cv * ddp1dtdr2du2) / R
                       + (u1 * u1 * ddr1dtdr2du2) / 2.0 
                       + r1 * u1 * ddu1dtdr2du2;
        dde1dtdp2dr2 = (Cv * ddp1dtdp2dr2) / R
                       + (u1 * u1 * ddr1dtdp2dr2) / 2.0 
                       + r1 * u1 * ddu1dtdp2dr2;

        dde1dtdr1du2 = dde1dtdu2dr1;
        dde1dtdu1du2 = dde1dtdu2du1;
        dde1dtdp1du2 = dde1dtdu2dp1;
        dde1dtdr2du2 = dde1dtdu2dr2;
        dde1dtdu2du2 = (Cv * ddp1dtdu2du2) / R 
                       + (u1 * u1 * ddr1dtdu2du2) / 2.0 
                       + r1 * u1 * ddu1dtdu2du2;
        dde1dtdp2du2 = (Cv * ddp1dtdp2du2) / R 
                       + (u1 * u1 * ddr1dtdp2du2) / 2.0 
                       + r1 * u1 * ddu1dtdp2du2;

        dde1dtdr1dp2 = dde1dtdp2dr1;
        dde1dtdu1dp2 = dde1dtdp2du1;
        dde1dtdp1dp2 = dde1dtdp2dp1;
        dde1dtdr2dp2 = dde1dtdp2dr2;
        dde1dtdu2dp2 = dde1dtdp2du2;
        dde1dtdp2dp2 = (Cv * ddp1dtdp2dp2) / R 
                       + (u1 * u1 * ddr1dtdp2dp2) / 2.0 
                       + r1 * u1 * ddu1dtdp2dp2;
// ***********************************************************************
        dRidWi(0,0) = dr1dtdr1;
        dRidWi(0,1) = dr1dtdu1;
        dRidWi(0,2) = dr1dtdp1;
        dRidWi(1,0) = du1dtdr1;
        dRidWi(1,1) = du1dtdu1;
        dRidWi(1,2) = du1dtdp1;
        dRidWi(2,0) = dp1dtdr1;
        dRidWi(2,1) = dp1dtdu1;
        dRidWi(2,2) = dp1dtdp1;

        dRidWd(0,0) = dr1dtdr2;
        dRidWd(0,1) = dr1dtdu2;
        dRidWd(0,2) = dr1dtdp2;
        dRidWd(1,0) = du1dtdr2;
        dRidWd(1,1) = du1dtdu2;
        dRidWd(1,2) = du1dtdp2;
        dRidWd(2,0) = dp1dtdr2;
        dRidWd(2,1) = dp1dtdu2;
        dRidWd(2,2) = dp1dtdp2;

        ddRindWdW[0](0, 0) = ddr1dtdr2dr2;
        ddRindWdW[0](1, 1) = ddr1dtdu2du2;
        ddRindWdW[0](2, 2) = ddr1dtdp2dp2;
        ddRindWdW[0](0, 1) = ddr1dtdr2du2;
        ddRindWdW[0](1, 0) = ddr1dtdu2dr2;
        ddRindWdW[0](0, 2) = ddr1dtdr2dp2;
        ddRindWdW[0](2, 0) = ddr1dtdp2dr2;
        ddRindWdW[0](1, 2) = ddr1dtdu2dp2;
        ddRindWdW[0](2, 1) = ddr1dtdp2du2;

        ddRindWdW[0](3, 3) = ddr1dtdr1dr1;
        ddRindWdW[0](4, 4) = ddr1dtdu1du1;
        ddRindWdW[0](5, 5) = ddr1dtdp1dp1;
        ddRindWdW[0](3, 4) = ddr1dtdr1du1;
        ddRindWdW[0](4, 3) = ddr1dtdu1dr1;
        ddRindWdW[0](3, 5) = ddr1dtdr1dp1;
        ddRindWdW[0](5, 3) = ddr1dtdp1dr1;
        ddRindWdW[0](4, 5) = ddr1dtdu1dp1;
        ddRindWdW[0](5, 4) = ddr1dtdp1du1;

        ddRindWdW[1](0, 0) = ddu1dtdr2dr2;
        ddRindWdW[1](1, 1) = ddu1dtdu2du2;
        ddRindWdW[1](2, 2) = ddu1dtdp2dp2;
        ddRindWdW[1](0, 1) = ddu1dtdr2du2;
        ddRindWdW[1](1, 0) = ddu1dtdu2dr2;
        ddRindWdW[1](0, 2) = ddu1dtdr2dp2;
        ddRindWdW[1](2, 0) = ddu1dtdp2dr2;
        ddRindWdW[1](1, 2) = ddu1dtdu2dp2;
        ddRindWdW[1](2, 1) = ddu1dtdp2du2;

        ddRindWdW[1](3, 3) = ddu1dtdr1dr1;
        ddRindWdW[1](4, 4) = ddu1dtdu1du1;
        ddRindWdW[1](5, 5) = ddu1dtdp1dp1;
        ddRindWdW[1](3, 4) = ddu1dtdr1du1;
        ddRindWdW[1](4, 3) = ddu1dtdu1dr1;
        ddRindWdW[1](3, 5) = ddu1dtdr1dp1;
        ddRindWdW[1](5, 3) = ddu1dtdp1dr1;
        ddRindWdW[1](4, 5) = ddu1dtdu1dp1;
        ddRindWdW[1](5, 4) = ddu1dtdp1du1;

        ddRindWdW[2](0, 0) = ddp1dtdr2dr2;
        ddRindWdW[2](1, 1) = ddp1dtdu2du2;
        ddRindWdW[2](2, 2) = ddp1dtdp2dp2;
        ddRindWdW[2](0, 1) = ddp1dtdr2du2;
        ddRindWdW[2](1, 0) = ddp1dtdu2dr2;
        ddRindWdW[2](0, 2) = ddp1dtdr2dp2;
        ddRindWdW[2](2, 0) = ddp1dtdp2dr2;
        ddRindWdW[2](1, 2) = ddp1dtdu2dp2;
        ddRindWdW[2](2, 1) = ddp1dtdp2du2;

        ddRindWdW[2](3, 3) = ddp1dtdr1dr1;
        ddRindWdW[2](4, 4) = ddp1dtdu1du1;
        ddRindWdW[2](5, 5) = ddp1dtdp1dp1;
        ddRindWdW[2](3, 4) = ddp1dtdr1du1;
        ddRindWdW[2](4, 3) = ddp1dtdu1dr1;
        ddRindWdW[2](3, 5) = ddp1dtdr1dp1;
        ddRindWdW[2](5, 3) = ddp1dtdp1dr1;
        ddRindWdW[2](4, 5) = ddp1dtdu1dp1;
        ddRindWdW[2](5, 4) = ddp1dtdp1du1;

        // Get Transformation Matrices
        MatrixXd dwpdw = dWpdW(W, 0);
        std::vector <Matrix3d> ddwpdwdwp = ddWpdWdWp(W, 0);

        MatrixXd temp(3, 3);
        for(int Ri = 0; Ri < 3; Ri++)
        {
            temp.setZero();
            temp += dwpdw.transpose() * ddRindWdW[Ri].topLeftCorner(3, 3);
            for(int Wpi = 0; Wpi < 3; Wpi++)
            {
                temp += dRodWd(Ri, Wpi) * ddwpdwdwp[Wpi];
            }
            ddRindWdW[Ri].topLeftCorner(3, 3) = dwpdw.transpose() * temp;
        }

        // Get Transformation Matrices
        dwpdw = dWpdW(W, 1);
        ddwpdwdwp = ddWpdWdWp(W, 2);

        for(int Ri = 0; Ri < 3; Ri++)
        {
            temp.setZero();
            temp += dwpdw.transpose() * ddRindWdW[Ri].bottomRightCorner(3, 3);
            for(int Wpi = 0; Wpi < 3; Wpi++)
            {
                temp += dRodWo(Ri, Wpi) * ddwpdwdwp[Wpi];
            }
            ddRindWdW[Ri].bottomRightCorner(3, 3) = dwpdw.transpose() * temp;
        }
        }

    // Supersonic Inlet
    else
    {
        for(int Ri = 0; Ri < 3; Ri++)
        {
            ddRindWdW[Ri].setZero();
            ddRindWdW[Ri].setZero();
        }
    }
}

// Hessian of Residual at the Boundaries
void HessianBC_FD(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW,
    std::vector <MatrixXd> &ddRoutdWdW)
{
    std::vector <double> Resi(3 * nx);
    double Resi0, Resi1, Resi2, Resi3, Resi4;
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 1e-3;
    double pertWi, pertWj;
    // Inlet
    int Ri = 0;
    int Rik;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;

        for(int Wi = 0; Wi <= (Ri + 1) * 3 + 2; Wi++)
        {
            pertWi = W[Wi] * h;
            for(int Wj = 0; Wj <= (Ri + 1) * 3 + 2; Wj++)
            {
                pertWj = W[Wj] * h;
                if(Wi != Wj)
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];
                    // R1
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi4 = Resi[Rik];

                    ddRindWdW[Rk](Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                        / (4 * pertWi * pertWj);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }
                else
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi0 = Resi[Rik];

                    // R1
                    Wd[Wi] = W[Wi] + 2.0 * pertWi;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wd[Wi] = W[Wi] - 2.0 * pertWi;

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[m] = W[m];
                    Resi4 = Resi[Rik];

                    ddRindWdW[Rk](Wi, Wj) = (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
    // Outlet
    Ri = nx - 1;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;

        for(int Wi = (Ri - 1) * 3; Wi <= Ri * 3 + 2; Wi++)
        {
            pertWi = W[Wi] * h;
            for(int Wj = (Ri - 1) * 3; Wj <= Ri * 3 + 2; Wj++)
            {
                pertWj = W[Wj] * h;
                if(Wi != Wj)
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];
                    // R1
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] + pertWj;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wd[Wi] = W[Wi] - pertWi;
                    Wd[Wj] = W[Wj] - pertWj;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi4 = Resi[Rik];

                    ddRoutdWdW[Rk](Wi - (Ri - 1) * 3, Wj - (Ri - 1) * 3)
                        = (Resi1 - Resi2 - Resi3 + Resi4) / (4 * pertWi * pertWj);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }
                else
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi0 = Resi[Rik];

                    // R1
                    Wd[Wi] = W[Wi] + 2.0 * pertWi;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wd[Wi] = W[Wi] + pertWi;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wd[Wi] = W[Wi] - pertWi;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wd[Wi] = W[Wi] - 2.0 * pertWi;

                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wd[3 * (nx - 1) + m] = W[3 * (nx - 1) + m];
                    Resi4 = Resi[Rik];

                    ddRoutdWdW[Rk](Wi - (Ri - 1) * 3, Wj - (Ri - 1) * 3) =
                        (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset Wd
                    Wd[Wi] = W[Wi];
                    Wd[Wj] = W[Wj];
                }// If Not Diagonal
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
//  std::cout<<"FD ddRoutdWpdWp"<<std::endl;
//  std::cout<<ddRoutdWdW[0]<<std::endl;
//  if(Min > 1.0)
//  {
//      for(int Rk = 0; Rk < 3; Rk++)
//      {
//          // Supersonic Inlet
//          for(int row = 0; row < 6; row++)
//          for(int col = 0; col < 6; col++)
//          {
//              // ddRindWindWin
//              ddRindWdW[Rk](row, col) = 0;
//              //if(row == col)
//              //    ddRindWdW[Rk](row, col) = 1;
//          }
//      }
//  }
}


// Jacobian at the Boundaries
void BCJac(
    std::vector <double> W,
    std::vector <double> dt,
    std::vector <double> dx,
    std::vector <double> &dBidWi,
    std::vector <double> &dBidWd,
    std::vector <double> &dBodWd,
    std::vector <double> &dBodWo)
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
    double r1, r2, p1, p2, u1, u2, c1, c2;
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

    deig1du1 = 0.5;
    deig1du2 = 0.5;

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
    dR1dp1 = eig1 * (c1 - 2.0 * dp * dc1dp1) / pow(c1, 3);
    dR1dr2 = eig1;
    dR1du2 = deig1du2 * (dp - c1 * c1 * drho) / (c1 * c1);
    dR1dp2 = - eig1 / (c1 * c1);

    dR2dr1 = - du * eig2 * (c1 + r1 * dc1dr1) - (dp + cr * du) * deig2dr1;
    dR2du1 = - cr * eig2 - (dp + cr * du) * deig2du1;
    dR2dp1 = - eig2 * (1.0 + du * r1 * dc1dp1) - (dp + cr * du) * deig2dp1;
    dR2dr2 = - (dp + cr * du) * deig2dr2;
    dR2du2 = cr * eig2 - (dp + cr * du) * deig2du2;
    dR2dp2 = eig2 - (dp + cr * du) * deig2dp2;

    dR3dr1 = eig3 * du * (c1 + r1 * dc1dr1) - (dp - cr * du) * deig3dr1;
    dR3du1 = cr * eig3 - (dp - cr * du) * deig3du1;
    dR3dp1 = - eig3 - dp * deig3dp1 + du * r1 * eig3 * dc1dp1;
    dR3dr2 = - (dp - cr * du) * deig3dr2;
    dR3du2 = - cr * eig3 - (dp - cr * du) * deig3du2;
    dR3dp2 = eig3 - (dp - cr * du) * deig3dp2;

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
//  double dru1dt;
//  dru1dt = r1 * du1dt + u1 * dr1dt;
    double dru1dtdr1, dru1dtdu1, dru1dtdp1;
    double dru1dtdr2, dru1dtdu2, dru1dtdp2;
    dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
    dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
    dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
    dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
    dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
    dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

    // de1/dt
//  double de1dt;
//  de1dt = dp1dt * Cv / R + u1 * r1 * du1dt + uu * dr1dt / 2.0;
    double de1dtdr1, de1dtdu1, de1dtdp1;
    double de1dtdr2, de1dtdu2, de1dtdp2;

    de1dtdr1 = dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1
               + du1dt * u1;
    de1dtdu1 = dp1dtdu1 * Cv / R + uu * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1
               + du1dt * r1 + dr1dt * u1;
    de1dtdp1 = dp1dtdp1 / (gam - 1) + uu * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
    de1dtdr2 = dp1dtdr2 / (gam - 1) + uu * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
    de1dtdu2 = dp1dtdu2 / (gam - 1) + uu * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
    de1dtdp2 = dp1dtdp2 / (gam - 1) + uu * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

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

    for(int row = 0; row < 3; row++)
    for(int col = 0; col < 3; col++)
    for(int k = 0; k < 3; k++)
        dBodWo[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];

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
    dWpdW(dwpdw, W, nx - 2);

    for(int row = 0; row < 3; row++)
    for(int col = 0; col < 3; col++)
    for(int k = 0; k < 3; k++)
        dBodWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];


    // *********************
    // INLET JACOBIANS
    // *********************
    // Subsonic Inlet
    if(u[0] < c[0])
    {
        // Values at time-step N
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

        // Shorthand
        drho = r2 - r1;
        dp = p2 - p1;
        du = u2 - u1;
        cr = r1 * c1;
        uu = u1 * u1;
        gamr = (gam - 1.0) / (gam + 1.0);
        fu = 1.0 - gamr * u1 * u1 / a2;

        // Speed of Sound
        dc1dr1 = - p1 * gam / (2.0 * cr * r1);
        dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
        dc1dp1 = gam / (2.0 * cr);
        dc2dp2 = gam / (2.0 * c2 * r2);

        // Eigenvalue
        eig1 = (u1 + u2) / 2.0;
        eig3 = eig1 - (c1 + c2) / 2.0;

        deig1du1 = 0.5;
        deig1du2 = 0.5;

        deig3dr1 = - dc1dr1 / 2.0;
        deig3du1 = deig1du1;
        deig3dp1 = - dc1dp1 / 2.0;
        deig3dr2 = - dc2dr2 / 2.0;
        deig3du2 = deig1du2;
        deig3dp2 = - dc2dp2 / 2.0;

        // Riemann Invariants
        R3 = - eig3 * (dp - cr * du);

        dR3dr1 = -eig3 * (-c1 * du - du * r1 * dc1dr1) - (dp - cr * du) * deig3dr1;
        dR3du1 = -cr * eig3 - (dp - cr * du) * deig3du1;
        dR3dp1 = eig3 * (1.0 + du * r1 * dc1dp1) - (dp - cr * du) * deig3dp1;
        dR3dr2 = -(dp - cr * du) * deig3dr2;
        dR3du2 = cr * eig3 - (dp - cr * du) * deig3du2;
        dR3dp2 = -eig3 - (dp - cr * du) * deig3dp2;
        // dp1
        double dp1du1, dp1du1du1;
        // Same Values
        dp1du1 = -2.0 * gamr * ptin * u1 * pow(fu, 1.0 / (gam - 1.0)) * gam
                 / (a2 * (gam - 1.0));

        dp1du1du1 = 2.0 * gamr * ptin * pow(fu, gam / (gam - 1.0)) * gam
                    * (a2 - a2 * gam + gamr * uu * (gam + 1))
                    / pow((a2 - gamr * uu) * (gam - 1.0), 2);

        // du1
        du1dt = R3 / (dp1du1 - cr);
        du1dtdr1 = dR3dr1 / (dp1du1 - cr)
                   - R3 * (-c1 - r1 * dc1dr1) / pow((dp1du1 - cr), 2);
        du1dtdu1 = dR3du1 / (dp1du1 - cr)
                   - R3 * dp1du1du1 / pow((dp1du1 - cr), 2);
        du1dtdp1 = dR3dp1 / (dp1du1 - cr)
                   + (R3 * r1 * dc1dp1) / pow((dp1du1 - cr), 2);
        du1dtdr2 = dR3dr2 / (dp1du1 - cr);
        du1dtdu2 = dR3du2 / (dp1du1 - cr);
        du1dtdp2 = dR3dp2 / (dp1du1 - cr);

        // Primitive values at time-step n+1
        double unp1, pnp1, rnp1, tnp1, funp1;
        unp1 = u1 + du1dt;
        pnp1 = ptin * pow(1 - gamr * pow(unp1, 2) / a2, gam / (gam - 1.0));
        tnp1 = Ttin * ( 1 - gamr * unp1 * unp1 / a2 );
        rnp1 = pnp1 / (R * tnp1);
        funp1 = 1.0 - gamr * unp1 * unp1 / a2;
        double dpnp1dunp1;
        dpnp1dunp1 = -2.0 * gamr * ptin * unp1 * pow(funp1, 1.0 / (gam - 1.0)) * gam
                     / (a2 * (gam - 1.0));

        // dp1
        dp1dt  = pnp1 - p1;
        dp1dtdr1 = dpnp1dunp1 * du1dtdr1;
        dp1dtdu1 = dpnp1dunp1 * (1.0 + du1dtdu1);// d(unp1)/du1 = du1du1 + du1dtdu1
        dp1dtdp1 = dpnp1dunp1 * du1dtdp1 - 1.0; // -1.0 due to dp1dp1
        dp1dtdr2 = dpnp1dunp1 * du1dtdr2;
        dp1dtdu2 = dpnp1dunp1 * du1dtdu2;
        dp1dtdp2 = dpnp1dunp1 * du1dtdp2;

        // dr1
        // Total derivative from rho_n+1 to p_n+1 and u_n+1
        double drnp1dpnp1, drnp1dtnp1, dtnp1dpnp1;
        drnp1dpnp1 = 1.0 / (R * tnp1);
        drnp1dtnp1 = -pnp1 / (R * tnp1 * tnp1);
        dtnp1dpnp1 = Ttin / ptin * (gam - 1.0) / gam * pow(pnp1 / ptin, - 1.0 / gam);
        double Drnp1Dpnp1 = drnp1dpnp1 + drnp1dtnp1 * dtnp1dpnp1;
        double drnp1dunp1 = Drnp1Dpnp1 * dpnp1dunp1;

        dr1dt = rnp1 - r1;

        dr1dtdr1 = drnp1dunp1 * du1dtdr1 - 1.0;
        dr1dtdu1 = drnp1dunp1 * (du1dtdu1 + 1.0);
        dr1dtdp1 = drnp1dunp1 * du1dtdp1;

        dr1dtdr2 = drnp1dunp1 * du1dtdr2;
        dr1dtdu2 = drnp1dunp1 * du1dtdu2;
        dr1dtdp2 = drnp1dunp1 * du1dtdp2;

        // dru1/dt
//      dru1dt = r1 * du1dt + u1 * dr1dt;

        dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
        dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
        dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
        dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
        dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
        dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

        // de1/dt
//      de1dt = dp1dt / (gam - 1.0) + r1 * u1 * du1dt + uu * dr1dt / 2.0;

        de1dtdr1 = dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1
                   + du1dt * u1;
        de1dtdu1 = dp1dtdu1 * Cv / R + uu * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1
                   + du1dt * r1 + dr1dt * u1;
        de1dtdp1 = dp1dtdp1 / (gam - 1.0) + uu * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
        de1dtdr2 = dp1dtdr2 / (gam - 1.0) + uu * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
        de1dtdu2 = dp1dtdu2 / (gam - 1.0) + uu * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
        de1dtdp2 = dp1dtdp2 / (gam - 1.0) + uu * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

        // Assign dR1/dWp1
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

        // Assign dR1/dWp2
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
        for(int i = 0; i < 9; i++)
        {
            dBidWi[i] = 0;
            dBidWd[i] = 0;
            if(i % 4 == 0)
                dBidWi[i] = 1;
        }
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

void HessianBCprim_FD(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW,
    std::vector <MatrixXd> &ddRoutdWdW)
{
    std::vector <double> Resi(3 * nx);
    double Resi0, Resi1, Resi2, Resi3, Resi4;
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Wp(3 * nx, 0);
    std::vector <double> Wpd(3 * nx, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 5e-4;
    double pertWi, pertWj;
    for(int i = 0; i < nx; i++)
    {
        Wp[i * 3 + 0] = W[i * 3 + 0];
        Wp[i * 3 + 1] = W[i * 3 + 1] / W[i * 3 + 0];
        Wp[i * 3 + 2] = (gam - 1.0) * ( W[i * 3 + 2] 
                        - (pow(W[i * 3 + 1], 2.0) / W[i * 3 + 0]) / 2.0 );
    }

    // Inlet
    int Ri = 0;
    int Rik;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;

        for(int Wi = 0; Wi <= (Ri + 1) * 3 + 2; Wi++)
        {
            pertWi = Wp[Wi] * h;
            for(int Wj = 0; Wj <= (Ri + 1) * 3 + 2; Wj++)
            {
                pertWj = Wp[Wj] * h;
                if(Wi != Wj)
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wpd[m] = Wp[m];
                    // R1
                    Wpd[Wi] = Wp[Wi] + pertWi;
                    Wpd[Wj] = Wp[Wj] + pertWj;
                    
                    PtoW(Wd, Wpd);

                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wpd[Wi] = Wp[Wi] + pertWi;
                    Wpd[Wj] = Wp[Wj] - pertWj;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wpd[Wi] = Wp[Wi] - pertWi;
                    Wpd[Wj] = Wp[Wj] + pertWj;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wpd[Wi] = Wp[Wi] - pertWi;
                    Wpd[Wj] = Wp[Wj] - pertWj;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi4 = Resi[Rik];

                    ddRindWdW[Rk](Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                        / (4 * pertWi * pertWj);
                    // Reset Wpd
                    Wpd[Wi] = Wp[Wi];
                    Wpd[Wj] = Wp[Wj];
                }
                else
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wpd[m] = Wp[m];

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi0 = Resi[Rik];

                    // R1
                    Wpd[Wi] = Wp[Wi] + 2.0 * pertWi;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wpd[Wi] = Wp[Wi] + pertWi;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wpd[Wi] = Wp[Wi] - pertWi;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wpd[Wi] = Wp[Wi] - 2.0 * pertWi;

                    PtoW(Wd, Wpd);
                    inletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[m] = Wp[m];
                    Resi4 = Resi[Rik];

                    ddRindWdW[Rk](Wi, Wj) = (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset Wd
                    Wpd[Wi] = Wp[Wi];
                    Wpd[Wj] = Wp[Wj];
                }
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
    // Outlet
    Ri = nx - 1;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;

        for(int Wi = (Ri - 1) * 3; Wi <= Ri * 3 + 2; Wi++)
        {
            pertWi = Wp[Wi] * h;
            for(int Wj = (Ri - 1) * 3; Wj <= Ri * 3 + 2; Wj++)
            {
                pertWj = Wp[Wj] * h;
                if(Wi != Wj)
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wpd[m] = Wp[m];
                    // R1
                    Wpd[Wi] = Wp[Wi] + pertWi;
                    Wpd[Wj] = Wp[Wj] + pertWj;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wpd[Wi] = Wp[Wi] + pertWi;
                    Wpd[Wj] = Wp[Wj] - pertWj;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wpd[Wi] = Wp[Wi] - pertWi;
                    Wpd[Wj] = Wp[Wj] + pertWj;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wpd[Wi] = Wp[Wi] - pertWi;
                    Wpd[Wj] = Wp[Wj] - pertWj;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi4 = Resi[Rik];

                    ddRoutdWdW[Rk](Wi - (Ri - 1) * 3, Wj - (Ri - 1) * 3)
                        = (Resi1 - Resi2 - Resi3 + Resi4) / (4 * pertWi * pertWj);
                    // Reset Wd
                    Wpd[Wi] = Wp[Wi];
                    Wpd[Wj] = Wp[Wj];
                }
                else
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wpd[m] = Wp[m];

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi0 = Resi[Rik];

                    // R1
                    Wpd[Wi] = Wp[Wi] + 2.0 * pertWi;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi1 = Resi[Rik];

                    // R2
                    Wpd[Wi] = Wp[Wi] + pertWi;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi2 = Resi[Rik];

                    // R3
                    Wpd[Wi] = Wp[Wi] - pertWi;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi3 = Resi[Rik];

                    // R4
                    Wpd[Wi] = Wp[Wi] - 2.0 * pertWi;

                    PtoW(Wd, Wpd);
                    outletBC(Wd, Resi, 1, 1);
                    for(int m = 0; m < 3; m++) Wpd[3 * (nx - 1) + m] = Wp[3 * (nx - 1) + m];
                    Resi4 = Resi[Rik];

                    ddRoutdWdW[Rk](Wi - (Ri - 1) * 3, Wj - (Ri - 1) * 3) =
                        (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset Wd
                    Wpd[Wi] = Wp[Wi];
                    Wpd[Wj] = Wp[Wj];
                }// If Not Diagonal
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
//  std::cout<<"FD ddRoutdWpdWp"<<std::endl;
//  std::cout<<ddRoutdWdW[0]<<std::endl;
//  if(Min > 1.0)
//  {
//      for(int Rk = 0; Rk < 3; Rk++)
//      {
//          // Supersonic Inlet
//          for(int row = 0; row < 6; row++)
//          for(int col = 0; col < 6; col++)
//          {
//              // ddRindWindWin
//              ddRindWdW[Rk](row, col) = 0;
//              //if(row == col)
//              //    ddRindWdW[Rk](row, col) = 1;
//          }
//      }
//  }
}
