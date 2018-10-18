#include "hessianOutlet.h"
#include<Eigen/Core>
#include<math.h>
#include<iostream>
#include"convert.h"
#include"structures.h"
//#include"residuald2.h"

using namespace Eigen;

void HessianOutlet(
    const struct Flow_options<double> &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddRoutdWdW)
{
    const int n_elem = flo_opts.n_elem;
    const double gam = flo_opts.gam;
    const double R = flo_opts.R;
    const double Cv = flo_opts.Cv;

    // First Derivatives Required for Second Derivatives
    Matrix3d dRodWd, dRodWo;
    for (int Rk = 0; Rk < 3; Rk++) {
        ddRoutdWdW[Rk].setZero();
    }

    // ************************
    // OUTLET HESSIAN
    // ************************

    const int i1 = n_elem - 1;
    const int i2 = n_elem - 2;
    const double r1 = W[i1*3+0];
    const double r2 = W[i2*3+0];
    const double p1 = get_p(gam, W[i1*3+0], W[i1*3+1], W[i1*3+2]);
    const double p2 = get_p(gam, W[i2*3+0], W[i2*3+1], W[i2*3+2]);
	const double u1 = W[i1*3+1] / W[i1*3+0];
	const double u2 = W[i2*3+1] / W[i2*3+0];
    const double c1 = get_c(gam, W[i1*3+0], W[i1*3+1], W[i1*3+2]);
    const double c2 = get_c(gam, W[i2*3+0], W[i2*3+1], W[i2*3+2]);

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

	UNUSED(ddR1dr2dr1);
	UNUSED(ddR1du2dr1);
	UNUSED(ddR1du2du1);
	UNUSED(ddR1dr1dp1);
	UNUSED(ddR1du1dp1);
	UNUSED(ddR1dr2dp1);
	UNUSED(ddR1du2dp1);
	UNUSED(ddR1dp2dp1);
	UNUSED(ddR1du1dr2);
	UNUSED(ddR1dr1dp2);
	UNUSED(ddR1du1dp2);
	UNUSED(ddR1dr2dp2);
	UNUSED(ddR1du2dp2);

	UNUSED(ddR2dr1dp1);
	UNUSED(ddR2du1dp1);
	UNUSED(ddR2du1dr2);
	UNUSED(ddR2dr1dp2);
	UNUSED(ddR2du1dp2);
	UNUSED(ddR2dr2dp2);
	UNUSED(ddR2du2dp2);

	UNUSED(ddR3dr1du1);
	UNUSED(ddR3dr1dp1);
	UNUSED(ddR3du1dp1);
	UNUSED(ddR3dr1dr2);
	UNUSED(ddR3du1dr2);
	UNUSED(ddR3dp1dr2);
	UNUSED(ddR3dr1du2);
	UNUSED(ddR3du1du2);
	UNUSED(ddR3dp1du2);
	UNUSED(ddR3dr2du2);
	UNUSED(ddR3dr1dp2);
	UNUSED(ddR3du1dp2);
	UNUSED(ddR3dp1dp2);
	UNUSED(ddR3dr2dp2);
	UNUSED(ddR3du2dp2);
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
    if (u1 < c1)
    {
        dp1dt = 0.0;

        dp1dtdr1 = 0.0;
        dp1dtdu1 = 0.0;
        dp1dtdp1 = 1.0;
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
	UNUSED(ddp1dtdr1dp1);
	UNUSED(ddp1dtdu1dp1);
	UNUSED(ddp1dtdu1dr2);
	UNUSED(ddp1dtdr1dp2);
	UNUSED(ddp1dtdu1dp2);
	UNUSED(ddp1dtdr2dp2);
	UNUSED(ddp1dtdu2dp2);
	UNUSED(ddp1dtdu1dr2);
	UNUSED(ddp1dtdr1dp2);
	UNUSED(ddp1dtdu1dp2);

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

	UNUSED(ddu1dtdu1dr2);
	UNUSED(ddu1dtdr1dp2);
	UNUSED(ddu1dtdu1dp2);

    // ***********************************************************************
    // d(ru)1/dt
//  double dru1dt;
//  dru1dt = r1 * du1dt + u1 * dr1dt + dr1dt * du1dt;
    double dru1dtdr1, dru1dtdu1, dru1dtdp1;
    double dru1dtdr2, dru1dtdu2, dru1dtdp2;
    dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1
                + dr1dtdr1 * du1dt + dr1dt * du1dtdr1;
    dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1
                + dr1dtdu1 * du1dt + dr1dt * du1dtdu1;
    dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1
                + dr1dtdp1 * du1dt + dr1dt * du1dtdp1;
    dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2
                + dr1dtdr2 * du1dt + dr1dt * du1dtdr2;
    dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2
                + dr1dtdu2 * du1dt + dr1dt * du1dtdu2;
    dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2
                + dr1dtdp2 * du1dt + dr1dt * du1dtdp2;

    // Second Derivative
    double
    ddru1dtdr1dr1, ddru1dtdu1dr1, ddru1dtdp1dr1, ddru1dtdr2dr1, ddru1dtdu2dr1, ddru1dtdp2dr1,
    ddru1dtdr1du1, ddru1dtdu1du1, ddru1dtdp1du1, ddru1dtdr2du1, ddru1dtdu2du1, ddru1dtdp2du1,
    ddru1dtdr1dp1, ddru1dtdu1dp1, ddru1dtdp1dp1, ddru1dtdr2dp1, ddru1dtdu2dp1, ddru1dtdp2dp1,
    ddru1dtdr1dr2, ddru1dtdu1dr2, ddru1dtdp1dr2, ddru1dtdr2dr2, ddru1dtdu2dr2, ddru1dtdp2dr2,
    ddru1dtdr1du2, ddru1dtdu1du2, ddru1dtdp1du2, ddru1dtdr2du2, ddru1dtdu2du2, ddru1dtdp2du2,
    ddru1dtdr1dp2, ddru1dtdu1dp2, ddru1dtdp1dp2, ddru1dtdr2dp2, ddru1dtdu2dp2, ddru1dtdp2dp2;

    // Top Left
    ddru1dtdr2dr2 = u1 * ddr1dtdr2dr2 + r1 * ddu1dtdr2dr2
                    + ddr1dtdr2dr2 * du1dt + dr1dtdr2 * du1dtdr2
                    + dr1dt * ddu1dtdr2dr2 + dr1dtdr2 * du1dtdr2;
    ddru1dtdu2du2 = u1 * ddr1dtdu2du2 + r1 * ddu1dtdu2du2
                    + ddr1dtdu2du2 * du1dt + dr1dtdu2 * du1dtdu2
                    + dr1dt * ddu1dtdu2du2 + dr1dtdu2 * du1dtdu2;
    ddru1dtdp2dp2 = u1 * ddr1dtdp2dp2 + r1 * ddu1dtdp2dp2
                    + ddr1dtdp2dp2 * du1dt + dr1dtdp2 * du1dtdp2
                    + dr1dt * ddu1dtdp2dp2 + dr1dtdp2 * du1dtdp2;

    ddru1dtdr2du2 = u1 * ddr1dtdr2du2 + r1 * ddu1dtdr2du2
                    + ddr1dtdr2du2 * du1dt + dr1dtdr2 * du1dtdu2
                    + dr1dt * ddu1dtdr2du2 + dr1dtdu2 * du1dtdr2;
    ddru1dtdu2dr2 = ddru1dtdr2du2;

    ddru1dtdr2dp2 = u1 * ddr1dtdr2dp2 + r1 * ddu1dtdr2dp2
                    + ddr1dtdr2dp2 * du1dt + dr1dtdr2 * du1dtdp2
                    + dr1dt * ddu1dtdr2dp2 + dr1dtdp2 * du1dtdr2;
    ddru1dtdp2dr2 = ddru1dtdr2dp2;

    ddru1dtdu2dp2 = u1 * ddr1dtdu2dp2 + r1 * ddu1dtdu2dp2
                    + ddr1dtdu2dp2 * du1dt + dr1dtdu2 * du1dtdp2
                    + dr1dt * ddu1dtdu2dp2 + dr1dtdp2 * du1dtdu2;
    ddru1dtdp2du2 = ddru1dtdu2dp2;

    // Bottom Right
    ddru1dtdr1dr1 = u1 * ddr1dtdr1dr1 + 2.0 * du1dtdr1 + r1 * ddu1dtdr1dr1
                    + ddr1dtdr1dr1 * du1dt + dr1dtdr1 * du1dtdr1
                    + dr1dt * ddu1dtdr1dr1 + dr1dtdr1 * du1dtdr1;
    ddru1dtdu1du1 = u1 * ddr1dtdu1du1 + 2.0 * dr1dtdu1 + r1 * ddu1dtdu1du1
                    + ddr1dtdu1du1 * du1dt + dr1dtdu1 * du1dtdu1
                    + dr1dt * ddu1dtdu1du1 + dr1dtdu1 * du1dtdu1;
    ddru1dtdp1dp1 = u1 * ddr1dtdp1dp1 + r1 * ddu1dtdp1dp1
                    + ddr1dtdp1dp1 * du1dt + dr1dtdp1 * du1dtdp1
                    + dr1dt * ddu1dtdp1dp1 + dr1dtdp1 * du1dtdp1;

    ddru1dtdr1du1 = u1 * ddr1dtdr1du1 + dr1dtdr1 + du1dtdu1 + r1 * ddu1dtdr1du1
                    + ddr1dtdr1du1 * du1dt + dr1dtdr1 * du1dtdu1
                    + dr1dt * ddu1dtdr1du1 + dr1dtdu1 * du1dtdr1;
    ddru1dtdu1dr1 = ddru1dtdr1du1;

    ddru1dtdr1dp1 = u1 * ddr1dtdp1dr1 + du1dtdp1 + r1 * ddu1dtdp1dr1
                    + ddr1dtdr1dp1 * du1dt + dr1dtdr1 * du1dtdp1
                    + dr1dt * ddu1dtdr1dp1 + dr1dtdp1 * du1dtdr1;
    ddru1dtdp1dr1 = ddru1dtdr1dp1;

    ddru1dtdu1dp1 = u1 * ddr1dtdp1du1 + dr1dtdp1 + r1 * ddu1dtdp1du1
                    + ddr1dtdu1dp1 * du1dt + dr1dtdu1 * du1dtdp1
                    + dr1dt * ddu1dtdu1dp1 + dr1dtdp1 * du1dtdu1;
    ddru1dtdp1du1 = ddru1dtdu1dp1;

    // Top Right

    ddru1dtdr2dr1 = u1 * ddr1dtdr1dr2 + du1dtdr2 + r1 * ddu1dtdr1dr2
                    + ddr1dtdr2dr1 * du1dt + dr1dtdr2 * du1dtdr1
                    + dr1dt * ddu1dtdr2dr1 + dr1dtdr1 * du1dtdr2;
    ddru1dtdu2du1 = u1 * ddr1dtdu2du1 + dr1dtdu2 + r1 * ddu1dtdu2du1
                    + ddr1dtdu2du1 * du1dt + dr1dtdu2 * du1dtdu1
                    + dr1dt * ddu1dtdu2du1 + dr1dtdu1 * du1dtdu2;
    ddru1dtdp2dp1 = u1 * ddr1dtdp2dp1 + r1 * ddu1dtdp2dp1
                    + ddr1dtdp2dp1 * du1dt + dr1dtdp2 * du1dtdp1
                    + dr1dt * ddu1dtdp2dp1 + dr1dtdp1 * du1dtdp2;

    ddru1dtdr2du1 = u1 * ddr1dtdr2du1 + dr1dtdr2 + r1 * ddu1dtdr2du1
                    + ddr1dtdr2du1 * du1dt + dr1dtdr2 * du1dtdu1
                    + dr1dt * ddu1dtdr2du1 + dr1dtdu1 * du1dtdr2;
    ddru1dtdu2dr1 = u1 * ddr1dtdu2dr1 + du1dtdu2 + r1 * ddu1dtdu2dr1
                    + ddr1dtdu2dr1 * du1dt + dr1dtdu2 * du1dtdr1
                    + dr1dt * ddu1dtdu2dr1 + dr1dtdr1 * du1dtdu2;

    ddru1dtdr2dp1 = u1 * ddr1dtdr2dp1 + r1 * ddu1dtdr2dp1
                    + ddr1dtdr2dp1 * du1dt + dr1dtdr2 * du1dtdp1
                    + dr1dt * ddu1dtdr2dp1 + dr1dtdp1 * du1dtdr2;
    ddru1dtdp2dr1 = u1 * ddr1dtdp2dr1 + du1dtdp2 + r1 * ddu1dtdp2dr1
                    + ddr1dtdp2dr1 * du1dt + dr1dtdp2 * du1dtdr1
                    + dr1dt * ddu1dtdp2dr1 + dr1dtdr1 * du1dtdp2;

    ddru1dtdu2dp1 = u1 * ddr1dtdu2dp1 + r1 * ddu1dtdu2dp1
                    + ddr1dtdu2dp1 * du1dt + dr1dtdu2 * du1dtdp1
                    + dr1dt * ddu1dtdu2dp1 + dr1dtdp1 * du1dtdu2;
    ddru1dtdp2du1 = u1 * ddr1dtdp2du1 + dr1dtdp2 + r1 * ddu1dtdp2du1
                    + ddr1dtdp2du1 * du1dt + dr1dtdp2 * du1dtdu1
                    + dr1dt * ddu1dtdp2du1 + dr1dtdu1 * du1dtdp2;


    // Bottom Left = Transpose of Top Right

    ddru1dtdr1dr2 = ddru1dtdr2dr1;
    ddru1dtdu1du2 = ddru1dtdu2du1;
    ddru1dtdp1dp2 = ddru1dtdp2dp1;

    ddru1dtdr1du2 = ddru1dtdu2dr1;
    ddru1dtdu1dr2 = ddru1dtdr2du1;

    ddru1dtdr1dp2 = ddru1dtdp2dr1;
    ddru1dtdp1dr2 = ddru1dtdr2dp1;

    ddru1dtdp1du2 = ddru1dtdu2dp1;
    ddru1dtdu1dp2 = ddru1dtdp2du1;


    // ***********************************************************************
    // de1/dt
//  double de1dt;
//  de1dt = dp1dt * Cv / R + u1 * r1 * du1dt + (u1 * u1) * dr1dt / 2.0;
    double de1dtdr1, de1dtdu1, de1dtdp1;
    double de1dtdr2, de1dtdu2, de1dtdp2;

    de1dtdr1 = (2.0 * Cv * dp1dtdr1
               + R * (du1dt * (du1dt + 2.0 * u1)
               + pow(du1dt + u1, 2.0) * dr1dtdr1
               + 2.0 * (dr1dt + r1) * (du1dt + u1) * du1dtdr1))
               / (2.0 * R);
    de1dtdu1 = (2.0 * Cv * dp1dtdu1
               + R * (pow(du1dt + u1, 2.0) * dr1dtdu1
               + 2.0 * (du1dt * r1 + dr1dt * (du1dt + u1)
               + (dr1dt + r1) * (du1dt + u1) * du1dtdu1)))
               / (2.0 * R);
    de1dtdp1 = (2.0 * Cv * dp1dtdp1
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp1
               + 2.0 * (dr1dt + r1) * du1dtdp1))
               / (2.0 * R);
    de1dtdr2 = (2.0 * Cv * dp1dtdr2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdr2
               + 2.0 * (dr1dt + r1) * du1dtdr2))
               / (2.0 * R);
    de1dtdu2 = (2.0 * Cv * dp1dtdu2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdu2
               + 2.0 * (dr1dt + r1) * du1dtdu2))
               / (2.0 * R);
    de1dtdp2 = (2.0 * Cv * dp1dtdp2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp2
               + 2.0 * (dr1dt + r1) * du1dtdp2))
               / (2.0 * R);

    // Second Derivative
    double
    dde1dtdr1dr1, dde1dtdu1dr1, dde1dtdp1dr1, dde1dtdr2dr1, dde1dtdu2dr1, dde1dtdp2dr1,
    dde1dtdr1du1, dde1dtdu1du1, dde1dtdp1du1, dde1dtdr2du1, dde1dtdu2du1, dde1dtdp2du1,
    dde1dtdr1dp1, dde1dtdu1dp1, dde1dtdp1dp1, dde1dtdr2dp1, dde1dtdu2dp1, dde1dtdp2dp1,
    dde1dtdr1dr2, dde1dtdu1dr2, dde1dtdp1dr2, dde1dtdr2dr2, dde1dtdu2dr2, dde1dtdp2dr2,
    dde1dtdr1du2, dde1dtdu1du2, dde1dtdp1du2, dde1dtdr2du2, dde1dtdu2du2, dde1dtdp2du2,
    dde1dtdr1dp2, dde1dtdu1dp2, dde1dtdp1dp2, dde1dtdr2dp2, dde1dtdu2dp2, dde1dtdp2dp2;

    dde1dtdr1dr1 = (2.0 * Cv * ddp1dtdr1dr1
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdr1dr1
                   + 2.0 * (2.0 * (du1dt + u1) * (1 + dr1dtdr1) * du1dtdr1
                   + (dr1dt + r1) * pow(du1dtdr1, 2.0)
                   + (dr1dt + r1) * (du1dt + u1) * ddu1dtdr1dr1)))
                   / (2.0 * R);
    dde1dtdu1dr1 = (2.0 * du1dt * R
                   + 2.0 * du1dt * R * du1dtdu1
                   + 2.0 * R * u1 * du1dtdu1
                   + 2.0 * R * (du1dt + u1) * dr1dtdr1 * (1 + du1dtdu1)
                   + 2.0 * R * du1dtdr1 * ((du1dt + u1) * dr1dtdu1
                   + (dr1dt + r1) * (1 + du1dtdu1))
                   + 2.0 * Cv * ddp1dtdr1du1
                   + pow(du1dt, 2.0) * R * ddr1dtdr1du1
                   + 2.0 * du1dt * R * u1 * ddr1dtdr1du1
                   + R * pow(u1, 2.0) * ddr1dtdr1du1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdr1du1
                   + 2.0 * du1dt * R * r1 * ddu1dtdr1du1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdr1du1
                   + 2.0 * R * r1 * u1 * ddu1dtdr1du1)
                   / (2.0 * R);
    dde1dtdp1dr1 = (2.0 * R * (du1dt + u1) * dr1dtdp1 * du1dtdr1
                   + 2.0 * R * du1dtdp1 * (du1dt + u1 + (du1dt + u1) * dr1dtdr1
                   + (dr1dt + r1) * du1dtdr1)
                   + 2.0 * Cv * ddp1dtdp1dr1
                   + pow(du1dt, 2.0) * R * ddr1dtdp1dr1
                   + 2.0 * du1dt * R * u1 * ddr1dtdp1dr1
                   + R * pow(u1, 2.0) * ddr1dtdp1dr1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdp1dr1
                   + 2.0 * du1dt * R * r1 * ddu1dtdp1dr1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdp1dr1
                   + 2.0 * R * r1 * u1 * ddu1dtdp1dr1)
                   / (2.0 * R);
    dde1dtdr2dr1 = (2.0 * R * (du1dt + u1) * dr1dtdr2 * du1dtdr1
                   + 2.0 * R * (du1dt + u1
                   + (du1dt + u1) * dr1dtdr1
                   + (dr1dt + r1) * du1dtdr1) * du1dtdr2
                   + 2.0 * Cv * ddp1dtdr1dr2
                   + pow(du1dt, 2.0) * R * ddr1dtdr1dr2
                   + 2.0 * du1dt * R * u1 * ddr1dtdr1dr2
                   + R * pow(u1, 2.0) * ddr1dtdr1dr2
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdr1dr2
                   + 2.0 * du1dt * R * r1 * ddu1dtdr1dr2
                   + 2.0 * dr1dt * R * u1 * ddu1dtdr1dr2
                   + 2.0 * R * r1 * u1 * ddu1dtdr1dr2)
                   / (2.0 * R);
    dde1dtdu2dr1 = (2.0 * R * (du1dt + u1) * dr1dtdu2 * du1dtdr1
                   + 2.0 * R * (du1dt + u1
                   + (du1dt + u1) * dr1dtdr1
                   + (dr1dt + r1) * du1dtdr1) * du1dtdu2
                   + 2.0 * Cv * ddp1dtdr1du2
                   + pow(du1dt, 2.0) * R * ddr1dtdr1du2
                   + 2.0 * du1dt * R * u1 * ddr1dtdr1du2
                   + R * pow(u1, 2.0) * ddr1dtdr1du2
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdr1du2
                   + 2.0 * du1dt * R * r1 * ddu1dtdr1du2
                   + 2.0 * dr1dt * R * u1 * ddu1dtdr1du2
                   + 2.0 * R * r1 * u1 * ddu1dtdr1du2)
                   / (2.0 * R);
    dde1dtdp2dr1 = (2.0 * R * (du1dt + u1) * dr1dtdp2 * du1dtdr1
                   + 2.0 * R * du1dtdp2 * (du1dt + u1 + (du1dt + u1) * dr1dtdr1
                   + (dr1dt + r1) * du1dtdr1)
                   + 2.0 * Cv * ddp1dtdp2dr1
                   + pow(du1dt, 2.0) * R * ddr1dtdp2dr1
                   + 2.0 * du1dt * R * u1 * ddr1dtdp2dr1
                   + R * pow(u1, 2.0) * ddr1dtdp2dr1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdp2dr1
                   + 2.0 * du1dt * R * r1 * ddu1dtdp2dr1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdp2dr1
                   + 2.0 * R * r1 * u1 * ddu1dtdp2dr1)
                   / (2.0 * R);

    dde1dtdr1du1 = dde1dtdu1dr1;
    dde1dtdu1du1 = (2.0 * Cv * ddp1dtdu1du1
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdu1du1
                   + 4.0 * (du1dt + u1) * dr1dtdu1 * (1 + du1dtdu1)
                   + 2.0 * (dr1dt + 2.0 * (dr1dt + r1) * du1dtdu1
                   + (dr1dt + r1) * pow(du1dtdu1, 2.0)
                   + (dr1dt + r1) * (du1dt + u1) * ddu1dtdu1du1)))
                   / (2.0 * R);
    dde1dtdp1du1 = (2.0 * R * (du1dt + u1) * dr1dtdp1 * (1 + du1dtdu1)
                   + 2.0 * R * du1dtdp1 * ((du1dt + u1) * dr1dtdu1
                   + (dr1dt + r1) * (1 + du1dtdu1))
                   + 2.0 * Cv * ddp1dtdp1du1
                   + pow(du1dt, 2.0) * R * ddr1dtdp1du1 + 2.0 * du1dt * R * u1 * ddr1dtdp1du1
                   + R * pow(u1, 2.0) * ddr1dtdp1du1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdp1du1
                   + 2.0 * du1dt * R * r1 * ddu1dtdp1du1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdp1du1
                   + 2.0 * R * r1 * u1 * ddu1dtdp1du1)
                   / (2.0 * R);
    dde1dtdr2du1 = (2.0 * R * (du1dt + u1) * dr1dtdr2 * (1 + du1dtdu1)
                   + 2.0 * R * du1dtdr2 * ((du1dt + u1) * dr1dtdu1
                   + (dr1dt + r1) * (1 + du1dtdu1))
                   + 2.0 * Cv * ddp1dtdr2du1
                   + pow(du1dt, 2.0) * R * ddr1dtdr2du1
                   + 2.0 * du1dt * R * u1 * ddr1dtdr2du1
                   + R * pow(u1, 2.0) * ddr1dtdr2du1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdr2du1
                   + 2.0 * du1dt * R * r1 * ddu1dtdr2du1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdr2du1
                   + 2.0 * R * r1 * u1 * ddu1dtdr2du1)
                   / (2.0 * R);
    dde1dtdu2du1 = (2.0 * R * (du1dt + u1) * dr1dtdu2 * (1 + du1dtdu1)
                   + 2.0 * R * ((du1dt + u1) * dr1dtdu1
                   + (dr1dt + r1) * (1 + du1dtdu1)) * du1dtdu2
                   + 2.0 * Cv * ddp1dtdu1du2
                   + pow(du1dt, 2.0) * R * ddr1dtdu1du2
                   + 2.0 * du1dt * R * u1 * ddr1dtdu1du2
                   + R * pow(u1, 2.0) * ddr1dtdu1du2
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdu1du2
                   + 2.0 * du1dt * R * r1 * ddu1dtdu1du2
                   + 2.0 * dr1dt * R * u1 * ddu1dtdu1du2
                   + 2.0 * R * r1 * u1 * ddu1dtdu1du2)
                   / (2.0 * R);
    dde1dtdp2du1 = (2.0 * R * (du1dt + u1) * dr1dtdp2 * (1 + du1dtdu1)
                   + 2.0 * R * du1dtdp2 * ((du1dt + u1) * dr1dtdu1
                   + (dr1dt + r1) * (1 + du1dtdu1))
                   + 2.0 * Cv * ddp1dtdp2du1
                   + pow(du1dt, 2.0) * R * ddr1dtdp2du1
                   + 2.0 * du1dt * R * u1 * ddr1dtdp2du1
                   + R * pow(u1, 2.0) * ddr1dtdp2du1
                   + 2.0 * dr1dt * du1dt * R * ddu1dtdp2du1
                   + 2.0 * du1dt * R * r1 * ddu1dtdp2du1
                   + 2.0 * dr1dt * R * u1 * ddu1dtdp2du1
                   + 2.0 * R * r1 * u1 * ddu1dtdp2du1)
                   / (2.0 * R);

    dde1dtdr1dp1 = dde1dtdp1dr1;
    dde1dtdu1dp1 = dde1dtdp1du1;
    dde1dtdp1dp1 = (2.0 * Cv * ddp1dtdp1dp1
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdp1dp1
                   + 2.0 * (2.0 * (du1dt + u1) * dr1dtdp1 * du1dtdp1
                   + (dr1dt + r1) * (pow(du1dtdp1, 2.0)
                   + (du1dt + u1) * ddu1dtdp1dp1))))/(2.0 * R);
    dde1dtdr2dp1 = (du1dt + u1) * dr1dtdr2 * du1dtdp1
                   + (du1dt + u1) * dr1dtdp1 * du1dtdr2
                   + dr1dt * du1dtdp1 * du1dtdr2
                   + r1 * du1dtdp1 * du1dtdr2
                   + (Cv * ddp1dtdp1dr2) / R
                   + (pow(du1dt, 2.0) * ddr1dtdp1dr2) / 2.0
                   + du1dt * u1 * ddr1dtdp1dr2
                   + (pow(u1, 2.0) * ddr1dtdp1dr2) / 2.0
                   + dr1dt * du1dt * ddu1dtdp1dr2
                   + du1dt * r1 * ddu1dtdp1dr2
                   + dr1dt * u1 * ddu1dtdp1dr2
                   + r1 * u1 * ddu1dtdp1dr2;
    dde1dtdu2dp1 = (du1dt + u1) * dr1dtdu2 * du1dtdp1
                   + (du1dt + u1) * dr1dtdp1 * du1dtdu2
                   + dr1dt * du1dtdp1 * du1dtdu2
                   + r1 * du1dtdp1 * du1dtdu2
                   + (Cv * ddp1dtdp1du2)/R
                   + (pow(du1dt, 2.0) * ddr1dtdp1du2) / 2.0
                   + du1dt * u1 * ddr1dtdp1du2
                   + (pow(u1, 2.0) * ddr1dtdp1du2) / 2.0
                   + dr1dt * du1dt * ddu1dtdp1du2
                   + du1dt * r1 * ddu1dtdp1du2
                   + dr1dt * u1 * ddu1dtdp1du2
                   + r1 * u1 * ddu1dtdp1du2;
    dde1dtdp2dp1 = (du1dt + u1) * dr1dtdp2 * du1dtdp1
                   + (du1dt + u1) * dr1dtdp1 * du1dtdp2
                   + dr1dt * du1dtdp1 * du1dtdp2
                   + r1 * du1dtdp1 * du1dtdp2
                   + (Cv * ddp1dtdp1dp2) / R
                   + (pow(du1dt, 2.0) * ddr1dtdp1dp2) / 2.0
                   + du1dt * u1 * ddr1dtdp1dp2
                   + (pow(u1, 2.0) * ddr1dtdp1dp2) / 2.0
                   + dr1dt * du1dt * ddu1dtdp1dp2
                   + du1dt * r1 * ddu1dtdp1dp2
                   + dr1dt * u1 * ddu1dtdp1dp2
                   + r1 * u1 * ddu1dtdp1dp2;

    dde1dtdr1dr2 = dde1dtdr2dr1;
    dde1dtdu1dr2 = dde1dtdr2du1;
    dde1dtdp1dr2 = dde1dtdr2dp1;
    dde1dtdr2dr2 = (2.0 * Cv * ddp1dtdr2dr2
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdr2dr2
                   + 2.0 * (2.0 * (du1dt + u1) * dr1dtdr2 * du1dtdr2
                   + (dr1dt + r1) * (pow(du1dtdr2, 2.0)
                   + (du1dt + u1) * ddu1dtdr2dr2))))/(2.0 * R);
    dde1dtdu2dr2 = (du1dt + u1) * dr1dtdu2 * du1dtdr2
                   + (du1dt + u1) * dr1dtdr2 * du1dtdu2
                   + dr1dt * du1dtdr2 * du1dtdu2
                   + r1 * du1dtdr2 * du1dtdu2
                   + (Cv * ddp1dtdr2du2)/R
                   + (pow(du1dt, 2.0) * ddr1dtdr2du2) / 2.0
                   + du1dt * u1 * ddr1dtdr2du2
                   + (pow(u1, 2.0) * ddr1dtdr2du2) / 2.0
                   + dr1dt * du1dt * ddu1dtdr2du2
                   + du1dt * r1 * ddu1dtdr2du2
                   + dr1dt * u1 * ddu1dtdr2du2
                   + r1 * u1 * ddu1dtdr2du2;
    dde1dtdp2dr2 = (du1dt + u1) * dr1dtdr2 * du1dtdp2
                   + (du1dt + u1) * dr1dtdp2 * du1dtdr2
                   + dr1dt * du1dtdp2 * du1dtdr2
                   + r1 * du1dtdp2 * du1dtdr2
                   + (Cv * ddp1dtdp2dr2)/R
                   + (pow(du1dt, 2.0) * ddr1dtdp2dr2) / 2.0
                   + du1dt * u1 * ddr1dtdp2dr2
                   + (pow(u1, 2.0) * ddr1dtdp2dr2) / 2.0
                   + dr1dt * du1dt * ddu1dtdp2dr2
                   + du1dt * r1 * ddu1dtdp2dr2
                   + dr1dt * u1 * ddu1dtdp2dr2
                   + r1 * u1 * ddu1dtdp2dr2;

    dde1dtdr1du2 = dde1dtdu2dr1;
    dde1dtdu1du2 = dde1dtdu2du1;
    dde1dtdp1du2 = dde1dtdu2dp1;
    dde1dtdr2du2 = dde1dtdu2dr2;
    dde1dtdu2du2 = (2.0 * Cv * ddp1dtdu2du2
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdu2du2
                   + 2.0 * (2.0 * (du1dt + u1) * dr1dtdu2 * du1dtdu2
                   + (dr1dt + r1) * (pow(du1dtdu2, 2.0)
                   + (du1dt + u1) * ddu1dtdu2du2))))/(2.0 * R);
    dde1dtdp2du2 = (du1dt + u1) * dr1dtdu2 * du1dtdp2
                   + (du1dt + u1) * dr1dtdp2 * du1dtdu2
                   + dr1dt * du1dtdp2 * du1dtdu2
                   + r1 * du1dtdp2 * du1dtdu2
                   + (Cv * ddp1dtdp2du2) / R
                   + (pow(du1dt, 2.0) * ddr1dtdp2du2) / 2.0
                   + du1dt * u1 * ddr1dtdp2du2
                   + (pow(u1, 2.0) * ddr1dtdp2du2) / 2.0
                   + dr1dt * du1dt * ddu1dtdp2du2
                   + du1dt * r1 * ddu1dtdp2du2
                   + dr1dt * u1 * ddu1dtdp2du2
                   + r1 * u1 * ddu1dtdp2du2;

    dde1dtdr1dp2 = dde1dtdp2dr1;
    dde1dtdu1dp2 = dde1dtdp2du1;
    dde1dtdp1dp2 = dde1dtdp2dp1;
    dde1dtdr2dp2 = dde1dtdp2dr2;
    dde1dtdu2dp2 = dde1dtdp2du2;
    dde1dtdp2dp2 = (2.0 * Cv * ddp1dtdp2dp2
                   + R * (pow(du1dt + u1, 2.0) * ddr1dtdp2dp2
                   + 2.0 * (2.0 * (du1dt + u1) * dr1dtdp2 * du1dtdp2
                   + (dr1dt + r1) * (pow(du1dtdp2, 2.0)
                   + (du1dt + u1) * ddu1dtdp2dp2))))
                   / (2.0 * R);


    //***********************************************************************
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

    // ****************
    // 1st Residual
    // ****************
    // Top Left
    ddRoutdWdW[0](0, 0) = ddr1dtdr2dr2;
    ddRoutdWdW[0](1, 1) = ddr1dtdu2du2;
    ddRoutdWdW[0](2, 2) = ddr1dtdp2dp2;

    ddRoutdWdW[0](0, 1) = ddr1dtdr2du2;
    ddRoutdWdW[0](1, 0) = ddr1dtdu2dr2;

    ddRoutdWdW[0](0, 2) = ddr1dtdr2dp2;
    ddRoutdWdW[0](2, 0) = ddr1dtdp2dr2;

    ddRoutdWdW[0](1, 2) = ddr1dtdu2dp2;
    ddRoutdWdW[0](2, 1) = ddr1dtdp2du2;

    // Bottom Right
    ddRoutdWdW[0](3, 3) = ddr1dtdr1dr1;
    ddRoutdWdW[0](4, 4) = ddr1dtdu1du1;
    ddRoutdWdW[0](5, 5) = ddr1dtdp1dp1;

    ddRoutdWdW[0](3, 4) = ddr1dtdr1du1;
    ddRoutdWdW[0](4, 3) = ddr1dtdu1dr1;

    ddRoutdWdW[0](3, 5) = ddr1dtdr1dp1;
    ddRoutdWdW[0](5, 3) = ddr1dtdp1dr1;

    ddRoutdWdW[0](4, 5) = ddr1dtdu1dp1;
    ddRoutdWdW[0](5, 4) = ddr1dtdp1du1;

    // Top Right
    ddRoutdWdW[0](0, 3) = ddr1dtdr2dr1;
    ddRoutdWdW[0](1, 4) = ddr1dtdu2du1;
    ddRoutdWdW[0](2, 5) = ddr1dtdp2dp1;

    ddRoutdWdW[0](0, 4) = ddr1dtdr2du1;
    ddRoutdWdW[0](1, 3) = ddr1dtdu2dr1;

    ddRoutdWdW[0](0, 5) = ddr1dtdr2dp1;
    ddRoutdWdW[0](2, 3) = ddr1dtdp2dr1;

    ddRoutdWdW[0](1, 5) = ddr1dtdu2dp1;
    ddRoutdWdW[0](2, 4) = ddr1dtdp2du1;


    // Bottom Left
    ddRoutdWdW[0](3, 0) = ddr1dtdr1dr2;
    ddRoutdWdW[0](4, 1) = ddr1dtdu1du2;
    ddRoutdWdW[0](5, 2) = ddr1dtdp1dp2;

    ddRoutdWdW[0](4, 0) = ddr1dtdu1dr2;
    ddRoutdWdW[0](3, 1) = ddr1dtdr1du2;

    ddRoutdWdW[0](5, 0) = ddr1dtdp1dr2;
    ddRoutdWdW[0](3, 2) = ddr1dtdr1dp2;

    ddRoutdWdW[0](5, 1) = ddr1dtdp1du2;
    ddRoutdWdW[0](4, 2) = ddr1dtdu1dp2;

    // ****************
    // 2nd Residual
    // ****************
    // Top Left
    ddRoutdWdW[1](0, 0) = ddru1dtdr2dr2;
    ddRoutdWdW[1](1, 1) = ddru1dtdu2du2;
    ddRoutdWdW[1](2, 2) = ddru1dtdp2dp2;

    ddRoutdWdW[1](0, 1) = ddru1dtdr2du2;
    ddRoutdWdW[1](1, 0) = ddru1dtdu2dr2;

    ddRoutdWdW[1](0, 2) = ddru1dtdr2dp2;
    ddRoutdWdW[1](2, 0) = ddru1dtdp2dr2;

    ddRoutdWdW[1](1, 2) = ddru1dtdu2dp2;
    ddRoutdWdW[1](2, 1) = ddru1dtdp2du2;

    // Bottom Right
    ddRoutdWdW[1](3, 3) = ddru1dtdr1dr1;
    ddRoutdWdW[1](4, 4) = ddru1dtdu1du1;
    ddRoutdWdW[1](5, 5) = ddru1dtdp1dp1;

    ddRoutdWdW[1](3, 4) = ddru1dtdr1du1;
    ddRoutdWdW[1](4, 3) = ddru1dtdu1dr1;

    ddRoutdWdW[1](3, 5) = ddru1dtdr1dp1;
    ddRoutdWdW[1](5, 3) = ddru1dtdp1dr1;

    ddRoutdWdW[1](4, 5) = ddru1dtdu1dp1;
    ddRoutdWdW[1](5, 4) = ddru1dtdp1du1;

    // Top Right
    ddRoutdWdW[1](0, 3) = ddru1dtdr2dr1;
    ddRoutdWdW[1](1, 4) = ddru1dtdu2du1;
    ddRoutdWdW[1](2, 5) = ddru1dtdp2dp1;

    ddRoutdWdW[1](0, 4) = ddru1dtdr2du1;
    ddRoutdWdW[1](1, 3) = ddru1dtdu2dr1;

    ddRoutdWdW[1](0, 5) = ddru1dtdr2dp1;
    ddRoutdWdW[1](2, 3) = ddru1dtdp2dr1;

    ddRoutdWdW[1](1, 5) = ddru1dtdu2dp1;
    ddRoutdWdW[1](2, 4) = ddru1dtdp2du1;


    // Bottom Left
    ddRoutdWdW[1](3, 0) = ddru1dtdr1dr2;
    ddRoutdWdW[1](4, 1) = ddru1dtdu1du2;
    ddRoutdWdW[1](5, 2) = ddru1dtdp1dp2;

    ddRoutdWdW[1](4, 0) = ddru1dtdu1dr2;
    ddRoutdWdW[1](3, 1) = ddru1dtdr1du2;

    ddRoutdWdW[1](5, 0) = ddru1dtdp1dr2;
    ddRoutdWdW[1](3, 2) = ddru1dtdr1dp2;

    ddRoutdWdW[1](5, 1) = ddru1dtdp1du2;
    ddRoutdWdW[1](4, 2) = ddru1dtdu1dp2;

    // ****************
    // 3rd Residual
    // ****************
    // Top Left
    ddRoutdWdW[2](0, 0) = dde1dtdr2dr2;
    ddRoutdWdW[2](1, 1) = dde1dtdu2du2;
    ddRoutdWdW[2](2, 2) = dde1dtdp2dp2;

    ddRoutdWdW[2](0, 1) = dde1dtdr2du2;
    ddRoutdWdW[2](1, 0) = dde1dtdu2dr2;

    ddRoutdWdW[2](0, 2) = dde1dtdr2dp2;
    ddRoutdWdW[2](2, 0) = dde1dtdp2dr2;

    ddRoutdWdW[2](1, 2) = dde1dtdu2dp2;
    ddRoutdWdW[2](2, 1) = dde1dtdp2du2;

    // Bottom Right
    ddRoutdWdW[2](3, 3) = dde1dtdr1dr1;
    ddRoutdWdW[2](4, 4) = dde1dtdu1du1;
    ddRoutdWdW[2](5, 5) = dde1dtdp1dp1;

    ddRoutdWdW[2](3, 4) = dde1dtdr1du1;
    ddRoutdWdW[2](4, 3) = dde1dtdu1dr1;

    ddRoutdWdW[2](3, 5) = dde1dtdr1dp1;
    ddRoutdWdW[2](5, 3) = dde1dtdp1dr1;

    ddRoutdWdW[2](4, 5) = dde1dtdu1dp1;
    ddRoutdWdW[2](5, 4) = dde1dtdp1du1;

    // Top Right
    ddRoutdWdW[2](0, 3) = dde1dtdr2dr1;
    ddRoutdWdW[2](1, 4) = dde1dtdu2du1;
    ddRoutdWdW[2](2, 5) = dde1dtdp2dp1;

    ddRoutdWdW[2](0, 4) = dde1dtdr2du1;
    ddRoutdWdW[2](1, 3) = dde1dtdu2dr1;

    ddRoutdWdW[2](0, 5) = dde1dtdr2dp1;
    ddRoutdWdW[2](2, 3) = dde1dtdp2dr1;

    ddRoutdWdW[2](1, 5) = dde1dtdu2dp1;
    ddRoutdWdW[2](2, 4) = dde1dtdp2du1;


    // Bottom Left
    ddRoutdWdW[2](3, 0) = dde1dtdr1dr2;
    ddRoutdWdW[2](4, 1) = dde1dtdu1du2;
    ddRoutdWdW[2](5, 2) = dde1dtdp1dp2;

    ddRoutdWdW[2](4, 0) = dde1dtdu1dr2;
    ddRoutdWdW[2](3, 1) = dde1dtdr1du2;

    ddRoutdWdW[2](5, 0) = dde1dtdp1dr2;
    ddRoutdWdW[2](3, 2) = dde1dtdr1dp2;

    ddRoutdWdW[2](5, 1) = dde1dtdp1du2;
    ddRoutdWdW[2](4, 2) = dde1dtdu1dp2;


    // Get Transformation Matrices
    MatrixXd dwpdw2 = eval_dWpdW(gam, W[i2*3+0], W[i2*3+1]);
    std::vector <Matrix3d> ddwpdwdwp = ddWpdWdWp(gam, W, n_elem - 2);

    MatrixXd temp(3, 3);
    for (int Ri = 0; Ri < 3; Ri++)
    {
        temp.setZero();
        temp += dwpdw2.transpose() * ddRoutdWdW[Ri].topLeftCorner(3, 3);
        for (int Wpi = 0; Wpi < 3; Wpi++)
        {
            temp += dRodWd(Ri, Wpi) * ddwpdwdwp[Wpi];
        }
        ddRoutdWdW[Ri].topLeftCorner(3, 3) = temp * dwpdw2;
    }

    // Get Transformation Matrices
    MatrixXd dwpdw1 = eval_dWpdW(gam, W[i1*3+0], W[i1*3+1]);
    ddwpdwdwp = ddWpdWdWp(gam, W, n_elem - 1);

    for (int Ri = 0; Ri < 3; Ri++)
    {
        temp.setZero();
        temp += dwpdw1.transpose() * ddRoutdWdW[Ri].bottomRightCorner(3, 3);
        for (int Wpi = 0; Wpi < 3; Wpi++)
        {
            temp += dRodWo(Ri, Wpi) * ddwpdwdwp[Wpi];
        }
        ddRoutdWdW[Ri].bottomRightCorner(3, 3) = temp * dwpdw1;
    }


    for (int Ri = 0; Ri < 3; Ri++)
    {
        ddRoutdWdW[Ri].topRightCorner(3, 3) =
            dwpdw1.transpose() * ddRoutdWdW[Ri].topRightCorner(3, 3) * dwpdw2;
    }
    for (int Ri = 0; Ri < 3; Ri++)
    {
        ddRoutdWdW[Ri].bottomLeftCorner(3, 3) =
            dwpdw2.transpose() * ddRoutdWdW[Ri].bottomLeftCorner(3, 3) * dwpdw1;
    }

}
