#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"

using namespace Eigen;

void HessianBC(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW,
    std::vector <MatrixXd> &ddRoutdWdW)
{
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
    double gamr, fu, drho, dp, du, cr, uu;
    drho = r1 - r2;
    dp = p1 - p2;
    du = u1 - u2;
    cr = r1 * c1;
    uu = u1 * u1;

    // ***********************************************************************
    // Speed of Sound
    double dc1dr1, dc2dr2, dc1dp1, dc2dp2;
    // First Derivative
    dc1dr1 = - p1 * gam / (2.0 * cr * r1);
    dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
    dc1dp1 = gam / (2.0 * cr);
    dc2dp2 = gam / (2.0 * c2 * r2);

    // Second Derivative
    ddc1dr1dr1 = 3.0 * c1 / (4.0 * r1 * r1);
    ddc1dp1dp1 = -c1 / (4.0 * p1 * p1);
    ddc1dr1dp1 = -gam / (4.0 * cr * r1);
    ddc1dp1dr1 = ddc1dr1dp1;

    ddc2dr2dr2 = 3.0 * c2 / (4.0 * r2 * r2);
    ddc2dp2dp2 = -c2 * (4 * p2 * p2);
    ddc2dr2dp2 = -gam / (4 * c2 * r2 * r2);
    ddc2dp2dr2 = ddc2dr2dp2;

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

    // ddeig1 = 0

    // ddeig2
    ddeig2dr1dr1 = ddc1dr1dr1 / 2.0;
    ddeig2dp1dp1 = ddc1dp1dp1 / 2.0;
    ddeig2dr1dp1 = ddc1dr1dp1 / 2.0;
    ddeig2dp1dr1 = ddeig2dr1dp1;

    ddeig2dr2dr2 = ddc2dr2dr2 / 2.0;
    ddeig2dp2dp2 = ddc2dp2dp2 / 2.0;
    ddeig2dr2dp2 = ddc2dr2dp2 / 2.0;
    ddeig2dp2dr2 = ddeig2dr2dp2;

    // ddeig3
    ddeig3dr1dr1 = -ddc1dr1dr1 / 2.0;
    ddeig3dp1dp1 = -ddc1dp1dp1 / 2.0;
    ddeig3dr1dp1 = -ddc1dr1dp1 / 2.0;
    ddeig3dp1dr1 = ddeig3dr1dp1;

    ddeig3dr2dr2 = -ddc2dr2dr2 / 2.0;
    ddeig3dp2dp2 = -ddc2dp2dp2 / 2.0;
    ddeig3dr2dp2 = -ddc2dr2dp2 / 2.0;
    ddeig3dp2dr2 = ddeig3dr2dp2;


    // ***********************************************************************
    // Riemann invariants
    double R1, R2, R3;
    double dR1dr1, dR1du1, dR1dp1, dR1dr2, dR1du2, dR1dp2;
    double dR2dr1, dR2du1, dR2dp1, dR2dr2, dR2du2, dR2dp2;
    double dR3dr1, dR3du1, dR3dp1, dR3dr2, dR3du2, dR3dp2;
    R1 = - eig1 * (drho - dp / (c1 * c1));
    R2 = - eig2 * (dp + cr * du);
    R3 = - eig3 * (dp - cr * du);

    // First Derivative of Riemann Invariant
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

    dR3dr1 = eig3 * du * (c1 + r1 * dc1dr1) - (dp - cr * du) * deig3dr1;
    dR3du1 = cr * eig3 - (dp - cr * du) * deig3du1;
    dR3dp1 = - eig3 - dp * deig3dp1 + du * r1 * eig3 * dc1dp1;
    dR3dr2 = - (dp - cr * du) * deig3dr2;
    dR3du2 = - cr * eig3 - (dp - cr * du) * deig3du2;
    dR3dp2 = eig3 - (dp - cr * du) * deig3dp2;

    // Second Derivative of Riemann Invariant
    double ddR1dr1dr1, ddR1du1dr1, ddR1dp1dr1, ddR1dr2dr1, ddR1du2dr1, ddR1dp2dr1,
           ddR1dr1du1, ddR1du1du1, ddR1dp1du1, ddR1dr2du1, ddR1du2du1, ddR1dp2du1,
           ddR1dr1dp1, ddR1du1dp1, ddR1dp1dp1, ddR1dr2dp1, ddR1du2dp1, ddR1dp2dp1,
           ddR1dr1dr2, ddR1du1dr2, ddR1dp1dr2, ddR1dr2dr2, ddR1du2dr2, ddR1dp2dr2,
           ddR1dr1du2, ddR1du1du2, ddR1dp1du2, ddR1dr2du2, ddR1du2du2, ddR1dp2du2,
           ddR1dr1dp2, ddR1du1dp2, ddR1dp1dp2, rdR1dr2dp2, ddR1du2dp2, ddR1dp2dp2;
    double ddR2dr1dr1, ddR2du1dr1, ddR2dp1dr1, ddR2dr2dr1, ddR2du2dr1, ddR2dp2dr1,
           ddR2dr1du1, ddR2du1du1, ddR2dp1du1, ddR2dr2du1, ddR2du2du1, ddR2dp2du1,
           ddR2dr1dp1, ddR2du1dp1, ddR2dp1dp1, ddR2dr2dp1, ddR2du2dp1, ddR2dp2dp1,
           ddR2dr1dr2, ddR2du1dr2, ddR2dp1dr2, ddR2dr2dr2, ddR2du2dr2, ddR2dp2dr2,
           ddR2dr1du2, ddR2du1du2, ddR2dp1du2, ddR2dr2du2, ddR2du2du2, ddR2dp2du2,
           ddR2dr1dp2, ddR2du1dp2, ddR2dp1dp2, rdR2dr2dp2, ddR2du2dp2, ddR2dp2dp2;
    double ddR3dr1dr1, ddR3du1dr1, ddR3dp1dr1, ddR3dr2dr1, ddR3du2dr1, ddR3dp2dr1,
           ddR3dr1du1, ddR3du1du1, ddR3dp1du1, ddR3dr2du1, ddR3du2du1, ddR3dp2du1,
           ddR3dr1dp1, ddR3du1dp1, ddR3dp1dp1, ddR3dr2dp1, ddR3du2dp1, ddR3dp2dp1,
           ddR3dr1dr2, ddR3du1dr2, ddR3dp1dr2, ddR3dr2dr2, ddR3du2dr2, ddR3dp2dr2,
           ddR3dr1du2, ddR3du1du2, ddR3dp1du2, ddR3dr2du2, ddR3du2du2, ddR3dp2du2,
           ddR3dr1dp2, ddR3du1dp2, ddR3dp1dp2, rdR3dr2dp2, ddR3du2dp2, ddR3dp2dp2;

    // R1
    ddR1dr1dr1 = -2.0 * eig1 * dp *
                 (-3.0 * dc1dr1 * dc1dr1 + c1 * ddc1dr1dr1)
                 / pow(c1, 4);
    ddR1du1dr1 = -(pow(c1, 3) + 2.0 * dp * dc1dr1) * deig1du1
                 / pow(c1, 3);
    ddR1dp1dr1 = -2.0 * eig1 *
                 ((c1 + 3 * (-p1 + p2) * dc1dp1) * dc1dr1
                 + c1 * (p1 - p2) * ddc1dp1dr1)
                 / pow(c1, 4);
    ddR1dr2dr1 = 0.0;
    ddR1du2dr1 = -(pow(c1, 3) + 2.0 * (p1 - p2) * dc1dr1) * deig1du2
                 / pow(c1, 3);
    ddR1dp2dr1 = 2.0 * eig1 * dc1dr1 / pow(c1, 3);

    ddR1dr1du1 = ddR1du1dr1;
    ddR1du1du1 = 0.0;
    ddR1dp1du1 = (c1 + 2.0 * (-p1 + p2) * dc1dp1) * deig1du1
                 / pow(c1, 3);
    ddR1dr2du1 = deig1du1;
    ddR1du2du1 = 0.0;
    ddR1dp2du1 = -deig1du1 / pow(c1, 2);

    ddR1dr1dp1 = ddR1dp1dr1;
    ddR1du1dp1 = ddR1dp1du1;
    ddR1dp1dp1 = 2.0 * eig1 *
                 (-2.0 * c1 * dc1dp1 + 3.0 * (p1 - p2) * dc1cp1 * dc1dp1
                 + c1 * (p2 - p1) * ddc1dp1dp1)
                 / pow(c1, 4);
    ddR1dr2dp1 = 0.0;
    ddR1du2dp1 = (c1 + 2.0 * (p2 - p1) * dc1dp1) * deig1du2 / pow(c1, 3);
    ddR1dp2dp1 = 2.0 * eig1 * dc1dp1 / pow(c1, 3);

    ddR1dr1dr2 = ddR1dr2dr1;
    ddR1du1dr2 = ddR1dr2du1;
    ddR1dp1dr2 = ddR1dr2dp1;
    ddR1dr2dr2 = 0.0;
    ddR1du2dr2 = deig1du2;
    ddR1dp2dr2 = 0.0;

    ddR1dr1du2 = ddR1du2dr1;
    ddR1du1du2 = ddR1du2du1;
    ddR1dp1du2 = ddR1du2dp1;
    ddR1dr2du2 = ddR1du2dr2;
    ddR1du2du2 = 0.0;
    ddR1dp2du2 = -deig1du2 / pow(c1, 2);

    ddR1dr1dp2 = ddR1dp1dr1;
    ddR1du1dp2 = ddR1dp1du1;
    ddR1dp1dp2 = ddR1dp1dp1;
    ddR1dr2dp2 = ddR1dp1dr2;
    ddR1du2dp2 = ddR1dp1du1;
    ddR1dp2dp2 = 0.0;

    // R2
    ddR2dr1dr1 = (p2 - p1) * ddeig2dr1dr1
                 - (u1 - u2) *
                 (
                  eig2 * r1 * ddc1dr1dr1
                  + 2.0 * c1 * deig2dr1
                  + 2.0 * dc1dr1 * (eig2 + r1 * deig2dr1)
                  + c1 * r1 * ddeig2dr1dr1
                 );
    ddR2du1dr1 = -r1 * dc1dr1 * (eig2 + (u1 - u2) * deig2du1)
                 -c1 * (eig2 + (u1 - u2) * deig2du2 + r1 * deig2dr1)
                 +(-p1 + p2 + c1 * (u2 - u1) * r1) * ddeig2du1dr1;
    ddR2dp1dr1 = -deig2dr1 + (p2 - p1) * ddeig2dp1dr1
                 - (u1 - u2) *
                 (
                  (c1 + r1 * dc1dr1) * deig2dp1
                  + dc1dp1 * (eig2 + r1 * deig2dr1)
                  + r1 * (eig2 * ddc1dp1dr1 + c1 * ddeig2dp1dr1)
                 );
    ddR2dr2dr1 = -(u1 - u2) * (c1 + r1 * dc1dr1) * deig2dr2
                 + (p2 - p1 + c1 * (u2 - u1) * r1) * deig2dr1dr2;
    ddR2du2dr1 = c1 * eig2 + eig2 * r1 * dc1dr1
                 - (u1 - u2) * (c1 + r1 * dc1dr1) * deig2du2
                 + c1 * r1 * deig2dr1
                 + (p2 - p1 + c1 * (u2 - u1) * r1) * ddeig2du2dr1;
    ddR2dp2dr1 = -(u1 - u2) * (c1 + r1 * dc1dr1) * deig2dp2;

    ddR2dr1du1 = ddR2du1dr1;
    ddR2du1du1 = -2.0 * c1 * r1 * deig2du1;
    ddR2dp1du1 = -deig2du1 + r1 * (c1 * deig2dp1
                 + dc1dp1 * (eig2 + (u1 - u2) * deig2du1));
    ddR2dr2du1 = -c1 * r1 * deig2dr2;
    ddR2du2du1 = c1 * r1 * deig2du1 - c1 * r1 * deig2du2;
    ddR2dp2du1 = -c1 * r1 * deig2dp2 + deig2du1;

    ddR2dr1dp1 = ddR2dp1dr1;
    ddR2du1dp1 = ddR2dp1du1;
    ddR2dp1dp1 = -2.0 * deig2dp1 + (p2 - p1) * ddeig2dp1dp1
                 -(u1 - u2) * r1 *
                 (
                  eig2 * ddc1dp1dp1
                  + 2.0 * dc1dp1 * deig2dp1
                  + c1 * ddeig2dp1dp1
                 );
    ddR2dr2dp1 = (-1.0 + (u2 - u1) * r1 * dc1dp1) * deig2dr2;
    ddR2du2dp1 = c1 * r1 * deig2dp1 - deig2du2
                 + r1 * dc1dp1 * eig2;
    ddR2dp2dp1 = deig2dp1 + (-1.0 + (u2 - u1) * r1 * dc1dp1) * deig2dp2

    ddR2dr1dr2 = ddR2dr2dr1;
    ddR2du1dr2 = ddR2dr2du1;
    ddR2dp1dr2 = ddR2dr2dp1;
    ddR2dr2dr2 = (p2 - p1 + c1 * (u2 - u1) * r1) * ddeig2dr2dr2;
    ddR2du2dr2 = c1 * r1 * deig2dr2;
    ddR2dp2dr2 = deig2dr2 + (p2 - p1 + c1 * (u2 - u1) * r1) * ddeig2dp2dr2;

    ddR2dr1du2 = ddR2du2dr1;
    ddR2du1du2 = ddR2du2du1;
    ddR2dp1du2 = ddR2du2dp1;
    ddR2dr2du2 = ddR2du2dr2;
    ddR2du2du2 = 2.0 * c1 * r1 * deig2du2;
    ddR2dp2du2 = c1 * r1 * deig2dp2 + deig2du2;

    ddR2dr1dp2 = ddR2dp1dr1;
    ddR2du1dp2 = ddR2dp1du1;
    ddR2dp1dp2 = ddR2dp1dp1;
    ddR2dr2dp2 = ddR2dp1dr2;
    ddR2du2dp2 = ddR2dp1du1;
    ddR2dp2dp2 = 2.0 * deig2dp2 + (p2 - p1 + c1 * (u2 - u1) * r1) * ddeig2dp2dp2;

    // R3
    ddR3dr1dr1 = (p2 - p1) * ddeig3dr1dr1
                 + (u1 - u2) *
                 (
                  eig3 * r1 * ddc1dr1dr1
                  + 2.0 * c1 * deig3dr1
                  + 2.0 * dc1dr1 * (eig3 + r1 * deig3dr1)
                  + c1 * r1 * ddeig3dr1dr1
                 );
    ddR3du1dr1 = c1 * eig3 + eig3 * r1 * dc1dr1
                 + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du1
                 + c1 * r1 * deig3dr1;
    ddR3dp1dr1 = -deig3dr1 + (p2 - p1) * ddeig3dp1dr1
                 + (u1 - u2) *
                 (
                  (c1 + r1 * dc1dr1) * deig3dp1
                  + dc1dp1 * (eig3 + r1 * deig3dr1)
                  + r1 * (eig3 * ddc1dp1dr1 + c1 * ddeig3dp1dr1)
                 );
    ddR3dr2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dr2;
    ddR3du2dr1 = -c1 * eig3 - eig3 * r1 * dc1dr1
                 + (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du2
                 - dc1dr1 * deig3dr1;
    ddR3dp2dr1 = (u1 - u2) * (c1 + r1 * dc1dr1) * deig3dp2 + deig3dr1;

    ddR3dr1du1 = ddR3du1dr1;
    ddR3du1du1 = 2.0 * c1 * r1 * deig3du1;
    ddR3dp1du1 = c1 * r1 * deig3dp1 - deig3du1
                 + r1 * dc1dp1 * (eig3 + (u1 - u2) * deig3du1);
    ddR3dr2du1 = c1 * r1 * deig3dr2;
    ddR3du2du1 = -c1 * r1 * deig3du1 + c1 * r1 * deig3du2;
    ddR3dp2du1 = c1 * r1 * deig3dp2 + deig3du1;

    ddR3dr1dp1 = ddR3dp1dr1;
    ddR3du1dp1 = ddR3dp1du1;
    ddR3dp1dp1 = -2.0 * deig3dp1 + (p2 - p1) * ddeig3dp1dp1
                 + (u1 - u2) * r1 *
                 (
                  eig3 * ddc1dp1dp1
                  + 2.0 * dc1dp1 * deig3dp1
                  + c1 * ddeig3dp1dp1
                 );
    ddR3dr2dp1 = (-1.0 + (u1 - u2) * r1 * dc1dp1) * deig3dr2;
    ddR3du2dp1 = -c1 * r1 * deig3dp1 - deig3du2
                 - r1 * dc1dp1 * (eig3  * (u2 - u1) * deig3du2);
    ddR3dp2dp1 = deig3dp1 + (-1.0 + (u1 - u2) * r1 * dc1dp1) * deig3dp2;

    ddR3dr1dr2 = ddR3dr2dr1;
    ddR3du1dr2 = ddR3dr2du1;
    ddR3dp1dr2 = ddR3dr2dp1;
    ddR3dr2dr2 = (p2 - p1 + c1 * (u1 - u2) * r1) * ddeig3dr2dr2;
    ddR3du2dr2 = -c1 * r1 * deig3dr2;
    ddR3dp2dr2 = deig3dr2 + (p2 - p1 + c1 * (u1 -u2) * r1) * deig3dp2dr2;

    ddR3dr1du2 = ddR3du2dr1;
    ddR3du1du2 = ddR3du2du1;
    ddR3dp1du2 = ddR3du2dp1;
    ddR3dr2du2 = ddR3du2dr2;
    ddR3du2du2 = -2.0 * c1 * r1 * deig3du2;
    ddR3dp2du2 = -c1 * r1 * deig3dp2 + deig3du2;

    ddR3dr1dp2 = ddR3dp1dr1;
    ddR3du1dp2 = ddR3dp1du1;
    ddR3dp1dp2 = ddR3dp1dp1;
    ddR3dr2dp2 = ddR3dp1dr2;
    ddR3du2dp2 = ddR3dp1du1;
    ddR3dp2dp2 = 2.0 * deig3dp2 + (p2 - p1 + c1 * (u1 - u2) * r1) * deig3dp2dp2;

    // ***********************************************************************
    // dp1/dt
    double dp1dt;
    double dp1dtdr1, dp1dtdu1, dp1dtdp1;
    double dp1dtdr2, dp1dtdu2, dp1dtdp2;
    double
    ddp1dtdr1dr1, ddp1dtdu1dr1, ddp1dtdp1dr1, ddp1dtdr2dr1, ddp1dtdu2dr1, ddp1dtdp2dt1,
    ddp1dtdr1du1, ddp1dtdu1du1, ddp1dtdp1du1, ddp1dtdr2du1, ddp1dtdu2du1, ddp1dtdp2du1,
    ddp1dtdr1dp1, ddp1dtdu1dp1, ddp1dtdp1dp1, ddp1dtdr2dp1, ddp1dtdu2dp1, ddp1dtdp2dp1,
    ddp1dtdr1dr2, ddp1dtdu1dr2, ddp1dtdp1dr2, ddp1dtdr2dr2, ddp1dtdu2dr2, ddp1dtdp2dr2,
    ddp1dtdr1du2, ddp1dtdu1du2, ddp1dtdp1du2, ddp1dtdr2du2, ddp1dtdu2du2, ddp1dtdp2du2,
    ddp1dtdr1dp2, ddp1dtdu1dp2, ddp1dtdp1dp2, ddp1dtdr2dp2, ddp1dtdu2dp2, ddp1dtdp2dp2;
    if(u1 < c1)
    {
        dp1dt = 0;

        ddp1dtdr1dr1 = 0.0;
        ddp1dtdu1dr1 = 0.0;
        ddp1dtdp1dr1 = 0.0;
        ddp1dtdr2dr1 = 0.0;
        ddp1dtdu2dr1 = 0.0;
        ddp1dtdp2dt1 = 0.0;

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

        ddp1dtdr1dr1 = (ddR2dr1dr1 + ddR3dr1dr1) / 2.0;
        ddp1dtdu1dr1 = (ddR2du1dr1 + ddR3du1dr1) / 2.0;
        ddp1dtdp1dr1 = (ddR2dp1dr1 + ddR3dp1dr1) / 2.0;
        ddp1dtdr2dr1 = (ddR2dr2dr1 + ddR3dr2dr1) / 2.0;
        ddp1dtdu2dr1 = (ddR2du2dr1 + ddR3du2dr1) / 2.0;
        ddp1dtdp2dt1 = (ddR2dp2dr1 + ddR3dp2dr1) / 2.0;

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
        ddp1dtdr2dr2 = (ddR2dr2dr1 + ddR3dr2dr1) / 2.0;
        ddp1dtdu2dr2 = (ddR2du2dr1 + ddR3du2dr1) / 2.0;
        ddp1dtdp2dr2 = (ddR2dp2dr1 + ddR3dp2dr1) / 2.0;

        ddp1dtdr1du2 = ddp1dtdr2dr1;
        ddp1dtdu1du2 = ddp1dtdr2du1;
        ddp1dtdp1du2 = ddp1dtdr2dp1;
        ddp1dtdr2du2 = ddp1dtdu2dr2;
        ddp1dtdu2du2 = (ddR2du2dr1 + ddR3du2dr1) / 2.0;
        ddp1dtdp2du2 = (ddR2dp2dr1 + ddR3dp2dr1) / 2.0;

        ddp1dtdr1dp2 = ddp1dtdr2dr1;
        ddp1dtdu1dp2 = ddp1dtdr2du1;
        ddp1dtdp1dp2 = ddp1dtdr2dp1;
        ddp1dtdr2dp2 = ddp1dtdp2dr2;
        ddp1dtdu2dp2 = ddp1dtdp2du2;
        ddp1dtdp2dp2 = (ddR2dp2dr1 + ddR3dp2dr1) / 2.0;
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
    ddr1dtdr1dr1, ddr1dtdu1dr1, ddr1dtdp1dr1, ddr1dtdr2dr1, ddr1dtdu2dr1, ddr1dtdp2dt1,
    ddr1dtdr1du1, ddr1dtdu1du1, ddr1dtdp1du1, ddr1dtdr2du1, ddr1dtdu2du1, ddr1dtdp2du1,
    ddr1dtdr1dp1, ddr1dtdu1dp1, ddr1dtdp1dp1, ddr1dtdr2dp1, ddr1dtdu2dp1, ddr1dtdp2dp1,
    ddr1dtdr1dr2, ddr1dtdu1dr2, ddr1dtdp1dr2, ddr1dtdr2dr2, ddr1dtdu2dr2, ddr1dtdp2dr2,
    ddr1dtdr1du2, ddr1dtdu1du2, ddr1dtdp1du2, ddr1dtdr2du2, ddr1dtdu2du2, ddr1dtdp2du2,
    ddr1dtdr1dp2, ddr1dtdu1dp2, ddr1dtdp1dp2, ddr1dtdr2dp2, ddr1dtdu2dp2, ddr1dtdp2dp2;

    ddr1dtdr1dr1 = (6.0 * dp1dt * dc1dr1 * dc1dr1
                   - 2.0 * c1 * dp1dt * ddc1dr1dr1
                   - 4.0 * c1 * dc1dr1 * dp1dtdr1
                   + c1 * c1 * ddp1dtdr1dr1) / pow(c1, 4)
                   + ddR1dr1dr1;
    ddr1dtdu1dr1 = (-2.0 * dc1dr1 * dp1dtdu1 
                   + c1 * ddp1dtdu1dr1) / pow(c1, 3)
                   + ddR1du1dr1;
    ddr1dtdp1dr1 = (dc1dp1 * ( 6.0 * dp1dt * dc1dr1 - 2.0 * c1 * dp1dtdr1)
                   + c1 * (-2.0 * dc1dr1 * dp1dtdp1 - 2.0 * dp1dt * ddc1dp1dr1
                           + c1 * ddp1dtdp1dr1)) / pow(c1, 4)
                   + ddR1dp1dr1;
    ddr1dtdr2dr1 = (-2.0 * dc1dr1 * dp1dtdr2
                   + c1 * ddp1dtdr1dr2) / pow(c1, 3)
                   + ddR1dr1dr2;
    ddr1dtdu2dr1 = (-2.0 * dc1dr1 * dp1dtdu2
                   + c1 * ddp1dtdu2dr1) / pow(c1, 3)
                   + ddR1du2dr1;
    ddr1dtdp2dt1 = (-2.0 * dc1dr1 * dp1dtdp2
                   + c1 * ddp1dtdp2dr1) / pow(c1, 3)
                   + ddR1dp2dr1;

    ddr1dtdr1du1 = ddr1dtdu1dr1;
    ddr1dtdu1du1 = ddp1dtdu1du1 / pow(c1, 2) + ddR1du1du1;
    ddr1dtdp1du1 = (-2.0 * dc1dp1 * dp1dtdu1
                   - 2.0 * dp1dt * ddc1dp1du1
                   + c1 * ddp1dtdp1du1) / pow(c1, 3)
                   + ddR1dp1du1;
    ddr1dtdr2du1 = ddp1du1dr2 / pow(c1, 2) + ddR1du1dr2;
    ddr1dtdu2du1 = ddp1du1du2 / pow(c1, 2) + ddR1du1du2;
    ddr1dtdp2du1 = ddp1dp2du1 / pow(c1, 2) + ddR1dp2du1;

    ddr1dtdr1dp1 = ddr1dtdp1dr1;
    ddr1dtdu1dp1 = ddr1dtdp1du1;
    ddr1dtdp1dp1 = (6.0 * dp1dt * dc1dp1 * dc1dp1
                   - 2.0 * c1 * dp1dt * ddc1dp1dp1
                   - 4.0 * c1 * dc1dp1 * dp1dtdp1
                   + c1 * c1 * ddp1dtdp1dp1) / pow(c1, 4)
                   + ddR1dp1dp1;
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
    ddr1dtdu2dr2 = ddp1dtdu2dr2 / pow(c1, 2) + ddR1du2dr2;
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
    du1dt = (R2 - dp1dt) / (cr);

    // First Derivative
    double du1dtdr1, du1dtdu1, du1dtdp1;
    double du1dtdr2, du1dtdu2, du1dtdp2;
    du1dtdr1 = ( (dp1dt - R2) * r1 * dc1dr1
               + c1 * (dp1dt - R2 - r1 * dp1dtdr1 + r1 * dR2dr1) )
               / (cr * cr);
    du1dtdu1 = (dR2du1 - dp1dtdu1) / cr;
    du1dtdp1 = ( (dp1dt - R2) * dc1dp1 + c1 * (dR2dp1 - dp1dtdp1) ) / (cr * c1);
    du1dtdr2 = (dR2dr2 - dp1dtdr2) / cr;
    du1dtdu2 = (dR2du2 - dp1dtdu2) / cr;
    du1dtdp2 = (dR2dp2 - dp1dtdp2) / cr;

    // Second Derivative
    double
    ddu1dtdr1dr1, ddu1dtdu1dr1, ddu1dtdp1dr1, ddu1dtdr2dr1, ddu1dtdu2dr1, ddu1dtdp2dt1,
    ddu1dtdr1du1, ddu1dtdu1du1, ddu1dtdp1du1, ddu1dtdr2du1, ddu1dtdu2du1, ddu1dtdp2du1,
    ddu1dtdr1dp1, ddu1dtdu1dp1, ddu1dtdp1dp1, ddu1dtdr2dp1, ddu1dtdu2dp1, ddu1dtdp2dp1,
    ddu1dtdr1dr2, ddu1dtdu1dr2, ddu1dtdp1dr2, ddu1dtdr2dr2, ddu1dtdu2dr2, ddu1dtdp2dr2,
    ddu1dtdr1du2, ddu1dtdu1du2, ddu1dtdp1du2, ddu1dtdr2du2, ddu1dtdu2du2, ddu1dtdp2du2,
    ddu1dtdr1dp2, ddu1dtdu1dp2, ddu1dtdp1dp2, ddu1dtdr2dp2, ddu1dtdu2dp2, ddu1dtdp2dp2;

    ddu1dtdr1dr1 = (2 * c1 * c1 * (-dp1dt + R2) 
                   + r1 * (2 * r1 * (-dp1dt + R2) * dc1dr1 * dc1dr1 
                   + 2 * c1 * dc1dr1 * (-dp1dt 
                   + R2 + r1 * dp1dtdr1 - r1 * dR2dr1) 
                   + c1 * (r1 * (dp1dt - R2) * ddc1dr1dr1 
                   + c1 * (2 * dp1dtdr1 - r1 * ddp1dtdr1dr1 
                   - 2 * dR2dr1 + r1 * ddR2dr1dr1))))
                   / pow(c1 * r1, 3);

    ddu1dtdu1dr1 = ((c1 + r1 * dc1dr1) * dp1dtdu1 
                   - (c1 + r1 * dc1dr1) * dR2du1
                   + c1 * r1 * (-ddp1dtdr1du1 + dR2dr1du1))
                   / pow(c1 * r1, 2);

    ddu1dtdp1dr1 = (dc1dp1 * (2 * r1 * (-dp1dt + R2) * dc1dr1 
                   + c1 * (-dp1dt + R2 + r1 * dp1dtdr1
                   - r1 * dR2dr1)) + c1 * ((c1 + r1 * dc1dr1) * dp1dtdp1
                   - (c1 + r1 * dc1dr1) * dR2dp1 
                   + r1 * ((dp1dt - R2) * ddc1dp1dr1 
                   + c1 * (-ddp1dtdp1dr1 + ddR2dp1dr1))))
                   /(pow(c1, 3) * pow(r1, 2));

    ddu1dtdr2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdr2 
                   - (c1 + r1 * dc1dr1) * dR1dr2 
                   + c1 * r1 * (-ddp1dtdr1dr2 + ddR2dr1dr2)) 
                   / pow(c1 * r1, 2);

    ddu1dtdu2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdu2 
                   - (c1 + r1 * dc1dr1) * dR2du2
                   + c1 * r1 * (-ddp1dtdr1du2 + ddR2dr1du2))
                   / pow(c1 * r1, 2);

    ddu1dtdp2dr1 = ((c1 + r1 * dc1dr1) * dp1dtdp2
                   - (c1 + r1 * dc1dr1) * dR2dp2
                   + c1 * r1 * (-ddp1dtdp2dr1 + ddR2dp2dr1))
                   /pow(c1 * r1, 2);

    ddu1dtdr1du1 = ddu1dtdu1dr1;
    ddu1dtdu1du1 = (-ddp1dtdu1du1 + ddR2du1du1) / (c1 * r1)
    ddu1dtdp1du1 = (dc1dp1 * (dp1dtdu1 - dR2du1) 
                   + (dp1dt - R2) * ddc1dp1du1
                   + c1 * (-ddp1dtdp1du1 + ddR2dp1du1))
                   / (c1 * c1 * r1);
    ddu1dtdr2du1 = (-ddp1dtdr2du1 + ddR2dr2du1) / (c1 * r1);
    ddu1dtdu2du1 = (-ddp1dtdu1du2 + ddR2du1du2) / (c1 * r1);
    ddu1dtdp2du1 = (-ddp1dtdp2du1 + ddR2dp2du1) / (c1 * r1);

    ddu1dtdr1dp1 = ddu1dtdp1dr1;
    ddu1dtdu1dp1 = ddu1dtdp1du1;
    ddu1dtdp1dp1 = (2 * (-dp1dt + R2) * dc1dp1 * dc1dp1 
                   + 2 * c1 * dc1dp1 * (dp1dtdp1 - dR2dp1)
                   + c1 * ((dp1dt - R2) * ddc1dp1dp1 
                   + c1 * (-ddp1dtdp1dp1 + ddR2dp1dp1)))
                   / (pow(c1, 3) * 1);
    ddu1dtdr2dp1 = (dc1dp1 * (dp1dtdr2 - dR2dr2) 
                   + c1 * (-ddp1dtdp1dr2 + ddR2dp1dr2))
                   /(c1 * c1 * r1);
    ddu1dtdu2dp1 = (dc1dp1 * (dp1dtdu2 - dR2du2) 
                   + c1 * (-ddp1dtdp1du2 + ddR2dp1du2))
                   /(c1 * c1 * r1);
    ddu1dtdp2dp1 = (dc1dp1 * (dp1dtdp2 - dR2dp2)
                   + c1 * (-ddp1dtdp1dp2 + ddR2dp1dp2))
                   /(c1 * c1 * r1);

    ddu1dtdr1dr2 = ddu1dtdr2dr1;
    ddu1dtdu1dr2 = ddu1dtdr2du1;
    ddu1dtdp1dr2 = ddu1dtdr2dp1; 
    ddu1dtdr2dr2 = (-ddp1dtdr2dr2 + ddR2dr2dr2)/(c1*r1)
    ddu1dtdu2dr2 = (-ddp1dtdr2du2 + ddR2dr2du2)/(c1*r1)
    ddu1dtdp2dr2 = (-ddp1dtdp2dr2 + ddR2dp2dr2)/(c1*r1)

    ddu1dtdr1du2 = ddu1dtdu2dr1;
    ddu1dtdu1du2 = ddu1dtdu2du1;
    ddu1dtdp1du2 = ddu1dtdu2dp1;
    ddu1dtdr2du2 = ddu1dtdu2dr2;
    ddu1dtdu2du2 = (-ddp1dtdu2du2 + ddR2du2du2)/(c1*r1)
    ddu1dtdp2du2 = (-ddp1dtdp2du2 + ddR2dp2du2)/(c1*r1)


    ddu1dtdr1dp2 = ddu1dtdp2dr1;
    ddu1dtdu1dp2 = ddu1dtdp2du1;
    ddu1dtdp1dp2 = ddu1dtdp2dp1;
    ddu1dtdr2dp2 = ddu1dtdp2dr2;
    ddu1dtdu2dp2 = ddu1dtdp2du2;
    ddu1dtdp2dp2 = (-ddp1dtdp2dp2 + ddR2dp2dp2)/(c1*r1)

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
    ddru1dtdr1dr1, ddru1dtdu1dr1, ddru1dtdp1dr1, ddru1dtdr2dr1, ddru1dtdu2dr1, ddru1dtdp2dt1,
    ddru1dtdr1du1, ddru1dtdu1du1, ddru1dtdp1du1, ddru1dtdr2du1, ddru1dtdu2du1, ddru1dtdp2du1,
    ddru1dtdr1dp1, ddru1dtdu1dp1, ddru1dtdp1dp1, ddru1dtdr2dp1, ddru1dtdu2dp1, ddru1dtdp2dp1,
    ddru1dtdr1dr2, ddru1dtdu1dr2, ddru1dtdp1dr2, ddru1dtdr2dr2, ddru1dtdu2dr2, ddru1dtdp2dr2,
    ddru1dtdr1du2, ddru1dtdu1du2, ddru1dtdp1du2, ddru1dtdr2du2, ddru1dtdu2du2, ddru1dtdp2du2,
    ddru1dtdr1dp2, ddru1dtdu1dp2, ddru1dtdp1dp2, ddru1dtdr2dp2, ddru1dtdu2dp2, ddru1dtdp2dp2;

    ddru1dtdr1dr1 = u1 * ddr1dtdr1dr1 + 2 * du1dtdr1 + r1 * ddu1dtdr1dr1;
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

    // ***********************************************************************
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

    // Second Derivative
    double
    dde1dtdr1dr1, dde1dtdu1dr1, dde1dtdp1dr1, dde1dtdr2dr1, dde1dtdu2dr1, dde1dtdp2dt1,
    dde1dtdr1du1, dde1dtdu1du1, dde1dtdp1du1, dde1dtdr2du1, dde1dtdu2du1, dde1dtdp2du1,
    dde1dtdr1dp1, dde1dtdu1dp1, dde1dtdp1dp1, dde1dtdr2dp1, dde1dtdu2dp1, dde1dtdp2dp1,
    dde1dtdr1dr2, dde1dtdu1dr2, dde1dtdp1dr2, dde1dtdr2dr2, dde1dtdu2dr2, dde1dtdp2dr2,
    dde1dtdr1du2, dde1dtdu1du2, dde1dtdp1du2, dde1dtdr2du2, dde1dtdu2du2, dde1dtdp2du2,
    dde1dtdr1dp2, dde1dtdu1dp2, dde1dtdp1dp2, dde1dtdr2dp2, dde1dtdu2dp2, dde1dtdp2dp2;

    dde1dtdr1dr1 = (cv*Dt[dp1dt, {r1, 2}])/r + 
     (u1*(u1*Dt[dr1dt, {r1, 2}] + 4*Dt[du1dt, r1] + 2*r1*Dt[du1dt, {r1, 2}]))/2
    dde1dtdu1dr1 = du1dt + u1*Dt[dr1dt, r1] + r1*Dt[du1dt, r1] + u1*Dt[du1dt, u1] + 
     (cv*Dt[dp1dt, r1, u1])/r + (u1^2*Dt[dr1dt, r1, u1])/2 + r1*u1*Dt[du1dt, r1, u1]
    dde1dtdp1dr1 = u1*Dt[du1dt, p1] + (cv*Dt[dp1dt, p1, r1])/r + (u1^2*Dt[dr1dt, p1, r1])/2 + 
     r1*u1*Dt[du1dt, p1, r1]
    dde1dtdr2dr1 = u1*Dt[du1dt, r2] + (cv*Dt[dp1dt, r1, r2])/r + (u1^2*Dt[dr1dt, r1, r2])/2 + 
     r1*u1*Dt[du1dt, r1, r2]
    dde1dtdu2dr1 = u1*Dt[du1dt, u2] + (cv*Dt[dp1dt, r1, u2])/r + (u1^2*Dt[dr1dt, r1, u2])/2 + 
     r1*u1*Dt[du1dt, r1, u2]
    dde1dtdp2dr1 = u1*Dt[du1dt, p2] + (cv*Dt[dp1dt, p2, r1])/r + (u1^2*Dt[dr1dt, p2, r1])/2 + 
     r1*u1*Dt[du1dt, p2, r1]

    dde1dtdr1du1 = dde1dtdu1dr1;
    dde1dtdu1du1 = dr1dt + (cv*Dt[dp1dt, {u1, 2}])/r + 2*u1*Dt[dr1dt, u1] + (u1^2*Dt[dr1dt, {u1, 2}])/2 + 
     2*r1*Dt[du1dt, u1] + r1*u1*Dt[du1dt, {u1, 2}]
    dde1dtdp1du1 = u1*Dt[dr1dt, p1] + r1*Dt[du1dt, p1] + (cv*Dt[dp1dt, p1, u1])/r + 
     (u1^2*Dt[dr1dt, p1, u1])/2 + r1*u1*Dt[du1dt, p1, u1]
    dde1dtdr2du1 = u1*Dt[dr1dt, r2] + r1*Dt[du1dt, r2] + (cv*Dt[dp1dt, r2, u1])/r + 
     (u1^2*Dt[dr1dt, r2, u1])/2 + r1*u1*Dt[du1dt, r2, u1]
    dde1dtdu2du1 = u1*Dt[dr1dt, u2] + r1*Dt[du1dt, u2] + (cv*Dt[dp1dt, u1, u2])/r + 
     (u1^2*Dt[dr1dt, u1, u2])/2 + r1*u1*Dt[du1dt, u1, u2]
    dde1dtdp2du1 = u1*Dt[dr1dt, p2] + r1*Dt[du1dt, p2] + (cv*Dt[dp1dt, p2, u1])/r + 
     (u1^2*Dt[dr1dt, p2, u1])/2 + r1*u1*Dt[du1dt, p2, u1]


    dde1dtdr1dp1 = dde1dtdp1dr1;
    dde1dtdu1dp1 = dde1dtdp1du1;
    dde1dtdp1dp1 =(cv*Dt[dp1dt, {p1, 2}])/r + (u1^2*Dt[dr1dt, {p1, 2}])/2 + r1*u1*Dt[du1dt, {p1, 2}]
    dde1dtdr2dp1 =(cv*Dt[dp1dt, p1, r2])/r + (u1^2*Dt[dr1dt, p1, r2])/2 + r1*u1*Dt[du1dt, p1, r2]
    dde1dtdu2dp1 =(cv*Dt[dp1dt, p1, u2])/r + (u1^2*Dt[dr1dt, p1, u2])/2 + r1*u1*Dt[du1dt, p1, u2]
    dde1dtdp2dp1 =(cv*Dt[dp1dt, p1, p2])/r + (u1^2*Dt[dr1dt, p1, p2])/2 + r1*u1*Dt[du1dt, p1, p2]

    dde1dtdr1dr2 = dde1dtdr2dr1;
    dde1dtdu1dr2 = dde1dtdr2du1;
    dde1dtdp1dr2 = dde1dtdr2dp1; 
    dde1dtdr2dr2 =(cv*Dt[dp1dt, {r2, 2}])/r + (u1^2*Dt[dr1dt, {r2, 2}])/2 + r1*u1*Dt[du1dt, {r2, 2}]
    dde1dtdu2dr2 =(cv*Dt[dp1dt, r2, u2])/r + (u1^2*Dt[dr1dt, r2, u2])/2 + r1*u1*Dt[du1dt, r2, u2]
    dde1dtdp2dr2 =(cv*Dt[dp1dt, p2, r2])/r + (u1^2*Dt[dr1dt, p2, r2])/2 + r1*u1*Dt[du1dt, p2, r2]

    dde1dtdr1du2 = dde1dtdu2dr1;
    dde1dtdu1du2 = dde1dtdu2du1;
    dde1dtdp1du2 = dde1dtdu2dp1;
    dde1dtdr2du2 = dde1dtdu2dr2;
    dde1dtdu2du2 = (cv*Dt[dp1dt, {u2, 2}])/r + (u1^2*Dt[dr1dt, {u2, 2}])/2 + r1*u1*Dt[du1dt, {u2, 2}]
    dde1dtdp2du2 = (cv*Dt[dp1dt, p2, u2])/r + (u1^2*Dt[dr1dt, p2, u2])/2 + r1*u1*Dt[du1dt, p2, u2]

    dde1dtdr1dp2 = dde1dtdp2dr1;
    dde1dtdu1dp2 = dde1dtdp2du1;
    dde1dtdp1dp2 = dde1dtdp2dp1;
    dde1dtdr2dp2 = dde1dtdp2dr2;
    dde1dtdu2dp2 = dde1dtdp2du2;
    dde1dtdp2dp2 = (cv*Dt[dp1dt, {p2, 2}])/r + (u1^2*Dt[dr1dt, {p2, 2}])/2 + r1*u1*Dt[du1dt, {p2, 2}]

    // ***********************************************************************
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
    int Rik, Rikp;
    // Inlet
    int Ri = 0;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        Rik = Ri * 3 + Rk;
        Rikp = (Ri + 1) * 3 + Rk;

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
    Ri = nx - 1;
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
    if(Min > 1.0)
    {
        for(int Rk = 0; Rk < 3; Rk++)
        {
            // Supersonic Inlet
            for(int row = 0; row < 6; row++)
            for(int col = 0; col < 6; col++)
            {
                // ddRindWindWin
                ddRindWdW[Rk](row, col) = 0;
                if(row == col)
                    ddRindWdW[Rk](row, col) = 1;
            }
        }
    }
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
        dR3dp1 = eig3 * (1 + du * r1 * dc1dp1) - (dp - cr * du) * deig3dp1;
        dR3dr2 = -(dp - cr * du) * deig3dr2;
        dR3du2 = cr * eig3 - (dp - cr * du) * deig3du2;
        dR3dp2 = -eig3 - (dp - cr * du) * deig3dp2;
        // dp1
        double dp1du1_n, dp1du1du1;
        // Same Values
        dp1du1_n = -2.0 * gamr * ptin * u1 * pow(fu, 1.0 / (gam - 1.0)) * gam
                 / (a2 * (gam - 1.0));

        dp1du1du1 = 2.0 * gamr * ptin * pow(fu, gam/(gam - 1.0)) * gam
                    * (a2 - a2 * gam + gamr * uu * (gam + 1))
                    / pow((a2 - gamr * uu) * (gam - 1.0), 2);

        // du1
        du1dt = R3 / (dp1du1_n - cr);
        du1dtdr1 = dR3dr1 / (dp1du1_n - cr)
                   - R3 * (-c1 - r1 * dc1dr1) / pow((dp1du1_n - cr), 2);
        du1dtdu1 = dR3du1 / (dp1du1_n - cr)
                   - R3 * dp1du1du1 / pow((dp1du1_n - cr), 2);
        du1dtdp1 = dR3dp1 / (dp1du1_n - cr)
                   + (R3 * r1 * dc1dp1) / pow((dp1du1_n - cr), 2);
        du1dtdr2 = dR3dr2 / (dp1du1_n - cr);
        du1dtdu2 = dR3du2 / (dp1du1_n - cr);
        du1dtdp2 = dR3dp2 / (dp1du1_n - cr);

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
        dp1dtdu1 = dpnp1dunp1 * (du1dtdu1 + 1);
        dp1dtdp1 = dpnp1dunp1 * du1dtdp1 - 1;
        dp1dtdr2 = dpnp1dunp1 * du1dtdr2;
        dp1dtdu2 = dpnp1dunp1 * du1dtdu2;
        dp1dtdp2 = dpnp1dunp1 * du1dtdp2;

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
