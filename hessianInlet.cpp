#include<Eigen/Core>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"residuald2.h"

using namespace Eigen;

void HessianInlet(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW)
{
    // First Derivatives Required for Second Derivatives
    Matrix3d dRidWi, dRidWd;
    for(int Rk = 0; Rk < 3; Rk++)
    {
        ddRindWdW[Rk].setZero();
    }

    std::vector <double> rho(nx), u(nx), e(nx), p(nx), c(nx), T(nx);
    WtoP(W, rho, u, e, p, c, T);

    // ************************
    // SUBSONIC INLET HESSIAN
    // ************************
    if(u[0] < c[0])
    {
        // Values at time-step N
        double i1, i2;
        double r1, r2, p1, p2, u1, u2, c1, c2;
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
        double gamr = (gam - 1.0) / (gam + 1.0);

        // Speed of Sound
        double dc1dr1, dc2dr2, dc1dp1, dc2dp2;
        dc1dr1 = - p1 * gam / (2.0 * r1 * c1 * r1);
        dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
        dc1dp1 = gam / (2.0 * r1 * c1);
        dc2dp2 = gam / (2.0 * c2 * r2);

        double ddc1dr1dr1, ddc1dp1dp1, ddc1dr1dp1,
               ddc2dr2dr2, ddc2dp2dp2, ddc2dr2dp2;
        ddc1dr1dr1 = 3.0 * c1 / (4.0 * r1 * r1);
        ddc1dp1dp1 = -c1 / (4.0 * p1 * p1);
        ddc1dr1dp1 = -gam / (4.0 * r1 * r1 * c1);

        ddc2dr2dr2 = 3.0 * c2 / (4.0 * r2 * r2);
        ddc2dp2dp2 = -c2 * (4.0 * p2 * p2);
        ddc2dr2dp2 = -gam / (4.0 * r2 * r2 * c2);

        // Eigenvalue
        double eig1, eig3;
        eig1 = (u1 + u2) / 2.0;
        eig3 = eig1 - (c1 + c2) / 2.0;

        double deig1du1, deig1du2;
        double deig3dr1, deig3du1, deig3dp1, deig3dr2, deig3du2, deig3dp2;
        deig1du1 = 0.5;
        deig1du2 = 0.5;

        deig3dr1 = - dc1dr1 / 2.0;
        deig3du1 = deig1du1;
        deig3dp1 = - dc1dp1 / 2.0;
        deig3dr2 = - dc2dr2 / 2.0;
        deig3du2 = deig1du2;
        deig3dp2 = - dc2dp2 / 2.0;

        double ddeig3dr1dr1, ddeig3dp1dp1, ddeig3dr1dp1;
        double ddeig3dr2dr2, ddeig3dp2dp2, ddeig3dr2dp2;
        ddeig3dr1dr1 = -ddc1dr1dr1 / 2.0;
        ddeig3dp1dp1 = -ddc1dp1dp1 / 2.0;
        ddeig3dr1dp1 = -ddc1dr1dp1 / 2.0;

        ddeig3dr2dr2 = -ddc2dr2dr2 / 2.0;
        ddeig3dp2dp2 = -ddc2dp2dp2 / 2.0;
        ddeig3dr2dp2 = -ddc2dr2dp2 / 2.0;

        // Riemann Invariants
        double R3;
        R3 = - eig3 * ((p2 - p1) - r1 * c1 * (u2 - u1));

        double dR3dr1, dR3du1, dR3dp1, dR3dr2, dR3du2, dR3dp2;
        dR3dr1 = -eig3 * (-c1 * (u2 - u1) - (u2 - u1) * r1 * dc1dr1)
                 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dr1;
        dR3du1 = -r1 * c1 * eig3 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3du1;
        dR3dp1 = eig3 * (1.0 + (u2 - u1) * r1 * dc1dp1)
                 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dp1;
        dR3dr2 = -((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dr2;
        dR3du2 = r1 * c1 * eig3 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3du2;
        dR3dp2 = -eig3 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dp2;

        double ddR3dr1dr1, ddR3du1dr1, ddR3dp1dr1, ddR3dr2dr1, ddR3du2dr1, ddR3dp2dr1,
               ddR3dr1du1, ddR3du1du1, ddR3dp1du1, ddR3dr2du1, ddR3du2du1, ddR3dp2du1,
               ddR3dr1dp1, ddR3du1dp1, ddR3dp1dp1, ddR3dr2dp1, ddR3du2dp1, ddR3dp2dp1,
               ddR3dr1dr2, ddR3du1dr2, ddR3dp1dr2, ddR3dr2dr2, ddR3du2dr2, ddR3dp2dr2,
               ddR3dr1du2, ddR3du1du2, ddR3dp1du2, ddR3dr2du2, ddR3du2du2, ddR3dp2du2,
               ddR3dr1dp2, ddR3du1dp2, ddR3dp1dp2, ddR3dr2dp2, ddR3du2dp2, ddR3dp2dp2;
        ddR3dr1dr1 = (p1 - p2) * ddeig3dr1dr1
                     - (u1 - u2) * (eig3 * r1 * ddc1dr1dr1
                     + 2.0 * c1 * deig3dr1
                     + 2.0 * dc1dr1 * (eig3 + r1 * deig3dr1)
                     + c1 * r1 * ddeig3dr1dr1);

        ddR3du1dr1 = - (r1 * dc1dr1 * (eig3 + (u1 - u2) * deig3du1))
                     - c1 * (eig3 + r1 * deig3dr1 + (u1 - u2) * deig3du1);

        ddR3dp1dr1 = deig3dr1 + (p1 - p2) * ddeig3dr1dp1
                     - (u1 - u2) * ((c1 + r1 * dc1dr1) * deig3dp1
                     + dc1dp1 * (eig3 + r1 * deig3dr1)
                     + r1 * (eig3 * ddc1dr1dp1 + c1 * ddeig3dr1dp1));
        ddR3dr2dr1 = -((u1 - u2) * (c1 + r1 * dc1dr1) * deig3dr2);
        ddR3du2dr1 = c1 * eig3 + eig3 * r1 * dc1dr1
                     + c1 * r1 * deig3dr1
                     - (u1 - u2) * (c1 + r1 * dc1dr1) * deig3du2;
        ddR3dp2dr1 = -((u1 - u2) * (c1 + r1 * dc1dr1) * deig3dp2)
                     - deig3dr1;

        ddR3dr1du1 = ddR3du1dr1;
        ddR3du1du1 = - 2.0 * c1 * r1 * deig3du1;
        ddR3dp1du1 = deig3du1
                     - r1 * (c1 * deig3dp1 + dc1dp1 * (eig3 + (u1 - u2) * deig3du1));
        ddR3dr2du1 = - (c1 * r1 * deig3dr2);
        ddR3du2du1 = c1 * r1 * deig3du1 - c1 * r1 * deig3du2;
        ddR3dp2du1 = - (c1 * r1 * deig3dp2) - deig3du1;

        ddR3dr1dp1 = ddR3dp1dr1;
        ddR3du1dp1 = ddR3dp1du1;
        ddR3dp1dp1 = 2.0 * deig3dp1 + (p1 - p2) * ddeig3dp1dp1
                     - r1 * (u1 - u2) * (eig3 * ddc1dp1dp1
                     + 2.0 * dc1dp1 * deig3dp1 + c1 * ddeig3dp1dp1);
        ddR3dr2dp1 = (1.0 + r1 * (-u1 + u2) * dc1dp1) * deig3dr2;
        ddR3du2dp1 = c1 * r1 * deig3dp1 + deig3du2
                     + r1 * dc1dp1 * (eig3 + (-u1 + u2) * deig3du2);
        ddR3dp2dp1 = -deig3dp1 + (1.0 + r1 * (-u1 + u2) * dc1dp1) * deig3dp2;

        ddR3dr1dr2 = ddR3dr2dr1;
        ddR3du1dr2 = ddR3dr2du1;
        ddR3dp1dr2 = ddR3dr2dp1;
        ddR3dr2dr2 = (p1 - p2 + c1 * r1 * (-u1 + u2)) * ddeig3dr2dr2;
        ddR3du2dr2 = c1 * r1 * deig3dr2;
        ddR3dp2dr2 = -deig3dr2
                     + (p1 - p2 + c1 * r1 * (-u1 + u2)) * ddeig3dr2dp2;

        ddR3dr1du2 = ddR3du2dr1;
        ddR3du1du2 = ddR3du2du1;
        ddR3dp1du2 = ddR3du2dp1;
        ddR3dr2du2 = ddR3du2dr2;
        ddR3du2du2 = 2.0 * c1 * r1 * deig3du2;
        ddR3dp2du2 = c1 * r1 * deig3dp2 - deig3du2;

        ddR3dr1dp2 = ddR3dp2dr1;
        ddR3du1dp2 = ddR3dp2du1;
        ddR3dp1dp2 = ddR3dp2dp1;
        ddR3dr2dp2 = ddR3dp2dr2;
        ddR3du2dp2 = ddR3dp2du2;
        ddR3dp2dp2 = -2.0 * deig3dp2
                     + (p1 - p2 + c1 * r1 * (-u1 + u2)) * ddeig3dp2dp2;

		UNUSED(ddR3dr1dp1);
		UNUSED(ddR3dr1dp2);
		UNUSED(ddR3du1dp1);
		UNUSED(ddR3du1dr2);
		UNUSED(ddR3du1dp2);
		UNUSED(ddR3dr2dp2);
		UNUSED(ddR3du2dp2);

        // dp1
        double dp1du1, ddp1du1du1, dddp1du1du1du1;
        // Same Values
        dp1du1 = (inlet_total_p * u1 * pow(1.0 - (gamr * u1 * u1) / a2, gam/(-1.0 + gam))
                 * gam * (2.0 * gamr)) / ((-a2 + gamr * u1 * u1) * (-1.0 + gam));

        ddp1du1du1 = (2.0 * gamr * inlet_total_p *
                      pow(1.0 - (gamr * pow(u1, 2)) / a2
                      , gam/(-1 + gam)) * gam *
                      (a2 - a2 * gam + gamr * pow(u1, 2) * (1 + gam)))
                      /(pow(a2 - gamr * pow(u1,2), 2) * pow(-1 + gam, 2));

        dddp1du1du1du1 = (4.0 * pow(gamr, 2) * inlet_total_p * u1 *
                         pow(1 - (gamr * pow(u1, 2)) / a2
                         , gam / (-1 + gam)) * gam *
                         (-3.0 * a2 * (-1.0 + gam)
                         + gamr * pow(u1, 2) * (1.0 + gam)))
                         /(pow(-a2 + gamr * pow(u1, 2), 3) *
                         pow(-1.0 + gam, 3));

        // du1
        double du1dt;
        du1dt = R3 / (dp1du1 - r1 * c1);
        double du1dtdr1, du1dtdu1, du1dtdp1;
        double du1dtdr2, du1dtdu2, du1dtdp2;

        du1dtdr1 = dR3dr1 / (dp1du1 - r1 * c1)
                   - R3 * (-c1 - r1 * dc1dr1) / pow((dp1du1 - r1 * c1), 2);
        du1dtdu1 = dR3du1 / (dp1du1 - r1 * c1)
                   - R3 * ddp1du1du1 / pow((dp1du1 - r1 * c1), 2);
        du1dtdp1 = dR3dp1 / (dp1du1 - r1 * c1)
                   + (R3 * r1 * dc1dp1) / pow((dp1du1 - r1 * c1), 2);
        du1dtdr2 = dR3dr2 / (dp1du1 - r1 * c1);
        du1dtdu2 = dR3du2 / (dp1du1 - r1 * c1);
        du1dtdp2 = dR3dp2 / (dp1du1 - r1 * c1);

        double
        ddu1dtdr1dr1, ddu1dtdu1dr1, ddu1dtdp1dr1, ddu1dtdr2dr1, ddu1dtdu2dr1, ddu1dtdp2dr1,
        ddu1dtdr1du1, ddu1dtdu1du1, ddu1dtdp1du1, ddu1dtdr2du1, ddu1dtdu2du1, ddu1dtdp2du1,
        ddu1dtdr1dp1, ddu1dtdu1dp1, ddu1dtdp1dp1, ddu1dtdr2dp1, ddu1dtdu2dp1, ddu1dtdp2dp1,
        ddu1dtdr1dr2, ddu1dtdu1dr2, ddu1dtdp1dr2, ddu1dtdr2dr2, ddu1dtdu2dr2, ddu1dtdp2dr2,
        ddu1dtdr1du2, ddu1dtdu1du2, ddu1dtdp1du2, ddu1dtdr2du2, ddu1dtdu2du2, ddu1dtdp2du2,
        ddu1dtdr1dp2, ddu1dtdu1dp2, ddu1dtdp1dp2, ddu1dtdr2dp2, ddu1dtdu2dp2, ddu1dtdp2dp2;
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
                       / pow(dp1du1 - c1 * r1, 2.0);

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

		UNUSED(ddu1dtdu1dr2);
		UNUSED(ddu1dtdu1dp2);
		UNUSED(ddu1dtdr1dp2);

        // Primitive values at time-step n+1
        double unp1, pnp1, rnp1, tnp1;
        unp1 = u1 + du1dt;
        pnp1 = inlet_total_p * pow(1.0 - gamr * pow(unp1, 2) / a2, gam / (gam - 1.0));
        tnp1 = inlet_total_T * (1.0 - gamr * unp1 * unp1 / a2);
        rnp1 = pnp1 / (R * tnp1);
        double dpnp1dunp1, ddpnp1dunp1dunp1;
//      dpnp1dunp1 = -2.0 * gamr * inlet_total_p * unp1
//                   * pow((1.0 - gamr * unp1 * unp1 / a2), 1.0 / (gam - 1.0))
//                   * gam / (a2 * (gam - 1.0));

        dpnp1dunp1 = (inlet_total_p * unp1
                     * pow(1.0 - (gamr * unp1 * unp1) / a2,
                           gam/(-1.0 + gam))
                     * gam * (2.0 * gamr))
                     / ((-a2 + gamr * unp1 * unp1)
                     * (-1.0 + gam));
        ddpnp1dunp1dunp1 = (2.0 * gamr * inlet_total_p *
                           pow(1.0 - (gamr * pow(unp1, 2)) / a2,
                               gam / (-1.0 + gam)) * gam
                           * (a2 - a2 * gam + gamr * pow(unp1, 2)
                           * (1.0 + gam))) / (pow(a2 - gamr * pow(unp1, 2), 2)
                           * pow(-1.0 + gam, 2));


        double dunp1dr1, dunp1du1, dunp1dp1, dunp1dr2, dunp1du2, dunp1dp2;
        dunp1dr1 = du1dtdr1;
        dunp1du1 = 1.0 + du1dtdu1;
        dunp1dp1 = du1dtdp1;
        dunp1dr2 = du1dtdr2;
        dunp1du2 = du1dtdu2;
        dunp1dp2 = du1dtdp2;

        // dp1
        // First Derivative
        double dp1dt;
        dp1dt  = pnp1 - p1;
        double dp1dtdr1, dp1dtdu1, dp1dtdp1;
        double dp1dtdr2, dp1dtdu2, dp1dtdp2;
        dp1dtdr1 = dpnp1dunp1 * dunp1dr1;
        dp1dtdu1 = dpnp1dunp1 * dunp1du1;
        dp1dtdp1 = dpnp1dunp1 * dunp1dp1 - 1.0; // -1.0 due to dp1dp1
        dp1dtdr2 = dpnp1dunp1 * dunp1dr2;
        dp1dtdu2 = dpnp1dunp1 * dunp1du2;
        dp1dtdp2 = dpnp1dunp1 * dunp1dp2;

		UNUSED(dp1dt);
        
        // Second Derivative
        // NOTE: second derivative of (unp1) = second derivative of (du1dt)
        double
        ddp1dtdr1dr1, ddp1dtdu1dr1, ddp1dtdp1dr1, ddp1dtdr2dr1, ddp1dtdu2dr1, ddp1dtdp2dr1,
        ddp1dtdr1du1, ddp1dtdu1du1, ddp1dtdp1du1, ddp1dtdr2du1, ddp1dtdu2du1, ddp1dtdp2du1,
        ddp1dtdr1dp1, ddp1dtdu1dp1, ddp1dtdp1dp1, ddp1dtdr2dp1, ddp1dtdu2dp1, ddp1dtdp2dp1,
        ddp1dtdr1dr2, ddp1dtdu1dr2, ddp1dtdp1dr2, ddp1dtdr2dr2, ddp1dtdu2dr2, ddp1dtdp2dr2,
        ddp1dtdr1du2, ddp1dtdu1du2, ddp1dtdp1du2, ddp1dtdr2du2, ddp1dtdu2du2, ddp1dtdp2du2,
        ddp1dtdr1dp2, ddp1dtdu1dp2, ddp1dtdp1dp2, ddp1dtdr2dp2, ddp1dtdu2dp2, ddp1dtdp2dp2;

        ddp1dtdr1dr1 = dpnp1dunp1 * ddu1dtdr1dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1dr1;
        ddp1dtdu1dr1 = dpnp1dunp1 * ddu1dtdu1dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1du1;
        ddp1dtdp1dr1 = dpnp1dunp1 * ddu1dtdp1dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1dp1;
        ddp1dtdr2dr1 = dpnp1dunp1 * ddu1dtdr2dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1dr2;
        ddp1dtdu2dr1 = dpnp1dunp1 * ddu1dtdu2dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1du2;
        ddp1dtdp2dr1 = dpnp1dunp1 * ddu1dtdp2dr1
                       + ddpnp1dunp1dunp1 * dunp1dr1 * dunp1dp2;

        ddp1dtdr1du1 = ddp1dtdu1dr1;
        ddp1dtdu1du1 = dpnp1dunp1 * ddu1dtdu1du1
                       + ddpnp1dunp1dunp1 * dunp1du1 * dunp1du1;
        ddp1dtdp1du1 = dpnp1dunp1 * ddu1dtdp1du1
                       + ddpnp1dunp1dunp1 * dunp1du1 * dunp1dp1;
        ddp1dtdr2du1 = dpnp1dunp1 * ddu1dtdr2du1
                       + ddpnp1dunp1dunp1 * dunp1du1 * dunp1dr2;
        ddp1dtdu2du1 = dpnp1dunp1 * ddu1dtdu2du1
                       + ddpnp1dunp1dunp1 * dunp1du1 * dunp1du2;
        ddp1dtdp2du1 = dpnp1dunp1 * ddu1dtdp2du1
                       + ddpnp1dunp1dunp1 * dunp1du1 * dunp1dp2;

        ddp1dtdr1dp1 = ddp1dtdp1dr1;
        ddp1dtdu1dp1 = ddp1dtdp1du1;
        ddp1dtdp1dp1 = dpnp1dunp1 * ddu1dtdp1dp1
                       + ddpnp1dunp1dunp1 * dunp1dp1 * dunp1dp1;
        ddp1dtdr2dp1 = dpnp1dunp1 * ddu1dtdr2dp1
                       + ddpnp1dunp1dunp1 * dunp1dp1 * dunp1dr2;
        ddp1dtdu2dp1 = dpnp1dunp1 * ddu1dtdu2dp1
                       + ddpnp1dunp1dunp1 * dunp1dp1 * dunp1du2;
        ddp1dtdp2dp1 = dpnp1dunp1 * ddu1dtdp2dp1
                       + ddpnp1dunp1dunp1 * dunp1dp1 * dunp1dp2;

        ddp1dtdr1dr2 = ddp1dtdr2dr1;
        ddp1dtdu1dr2 = ddp1dtdr2du1;
        ddp1dtdp1dr2 = ddp1dtdr2dp1;
        ddp1dtdr2dr2 = dpnp1dunp1 * ddu1dtdr2dr2
                       + ddpnp1dunp1dunp1 * dunp1dr2 * dunp1dr2;
        ddp1dtdu2dr2 = dpnp1dunp1 * ddu1dtdu2dr2
                       + ddpnp1dunp1dunp1 * dunp1dr2 * dunp1du2;
        ddp1dtdp2dr2 = dpnp1dunp1 * ddu1dtdp2dr2
                       + ddpnp1dunp1dunp1 * dunp1dr2 * dunp1dp2;

        ddp1dtdr1du2 = ddp1dtdu2dr1;
        ddp1dtdu1du2 = ddp1dtdu2du1;
        ddp1dtdp1du2 = ddp1dtdu2dp1;
        ddp1dtdr2du2 = ddp1dtdu2dr2;
        ddp1dtdu2du2 = dpnp1dunp1 * ddu1dtdu2du2
                       + ddpnp1dunp1dunp1 * dunp1du2 * dunp1du2;
        ddp1dtdp2du2 = dpnp1dunp1 * ddu1dtdp2du2
                       + ddpnp1dunp1dunp1 * dunp1du2 * dunp1dp2;

        ddp1dtdr1dp2 = ddp1dtdp2dr1;
        ddp1dtdu1dp2 = ddp1dtdp2du1;
        ddp1dtdp1dp2 = ddp1dtdp2dp1;
        ddp1dtdr2dp2 = ddp1dtdp2dr2;
        ddp1dtdu2dp2 = ddp1dtdp2du2;
        ddp1dtdp2dp2 = dpnp1dunp1 * ddu1dtdp2dp2
                       + ddpnp1dunp1dunp1 * dunp1dp2 * dunp1dp2;

		UNUSED(ddp1dtdr1dp1);
		UNUSED(ddp1dtdu1dp1);
		UNUSED(ddp1dtdu1dr2);
		UNUSED(ddp1dtdr1dp2);
		UNUSED(ddp1dtdu1dp2);
		UNUSED(ddp1dtdr2dp2);
		UNUSED(ddp1dtdu2dp2);

        // dr1
        // Total derivative from rho_n+1 to p_n+1 and u_n+1
        double drnp1dpnp1, drnp1dtnp1, dtnp1dpnp1;
        drnp1dpnp1 = 1.0 / (R * tnp1);
        drnp1dtnp1 = -pnp1 / (R * tnp1 * tnp1);
        dtnp1dpnp1 = inlet_total_T / inlet_total_p * (gam - 1.0) / gam * pow(pnp1 / inlet_total_p, - 1.0 / gam);
        
        double Drnp1Dpnp1 = drnp1dpnp1 + drnp1dtnp1 * dtnp1dpnp1;
        double drnp1dunp1 = Drnp1Dpnp1 * dpnp1dunp1;

        drnp1dunp1 = 
            (-2.0 * a2 * (1.0 + gam) * inlet_total_p *unp1 
            * pow(1.0 + (pow(unp1, 2) - gam * pow(unp1, 2)) / (a2 + a2 * gam),
            gam/(-1 + gam))) / (R * inlet_total_T * pow(a2 + a2 * gam + pow(unp1, 2) 
            - gam * pow(unp1, 2), 2));
        double ddrnp1dunp1dunp1 = 
            (-2.0 * a2 * (1.0 + gam) * inlet_total_p * (a2 * (1.0 + gam) 
            + (-3.0 + gam) * pow(unp1, 2)) * pow(1.0 + (pow(unp1, 2) - gam 
            * pow(unp1, 2)) / (a2 + a2 * gam), gam/(-1.0 + gam))) 
            / (R * inlet_total_T * pow(a2 + a2 * gam 
            + pow(unp1, 2) - gam * pow(unp1, 2), 3));
    
    
        double dr1dt;
        dr1dt = rnp1 - r1;

        
        double dr1dtdr1, dr1dtdu1, dr1dtdp1;
        double dr1dtdr2, dr1dtdu2, dr1dtdp2;
        dr1dtdr1 = drnp1dunp1 * dunp1dr1 - 1.0;
        dr1dtdu1 = drnp1dunp1 * dunp1du1;
        dr1dtdp1 = drnp1dunp1 * dunp1dp1;
        dr1dtdr2 = drnp1dunp1 * dunp1dr2;
        dr1dtdu2 = drnp1dunp1 * dunp1du2;
        dr1dtdp2 = drnp1dunp1 * dunp1dp2;
    
        double
        ddr1dtdr1dr1, ddr1dtdu1dr1, ddr1dtdp1dr1, ddr1dtdr2dr1, ddr1dtdu2dr1, ddr1dtdp2dr1,
        ddr1dtdr1du1, ddr1dtdu1du1, ddr1dtdp1du1, ddr1dtdr2du1, ddr1dtdu2du1, ddr1dtdp2du1,
        ddr1dtdr1dp1, ddr1dtdu1dp1, ddr1dtdp1dp1, ddr1dtdr2dp1, ddr1dtdu2dp1, ddr1dtdp2dp1,
        ddr1dtdr1dr2, ddr1dtdu1dr2, ddr1dtdp1dr2, ddr1dtdr2dr2, ddr1dtdu2dr2, ddr1dtdp2dr2,
        ddr1dtdr1du2, ddr1dtdu1du2, ddr1dtdp1du2, ddr1dtdr2du2, ddr1dtdu2du2, ddr1dtdp2du2,
        ddr1dtdr1dp2, ddr1dtdu1dp2, ddr1dtdp1dp2, ddr1dtdr2dp2, ddr1dtdu2dp2, ddr1dtdp2dp2;

        ddr1dtdr1dr1 = drnp1dunp1 * ddu1dtdr1dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1dr1;
        ddr1dtdu1dr1 = drnp1dunp1 * ddu1dtdu1dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1du1;
        ddr1dtdp1dr1 = drnp1dunp1 * ddu1dtdp1dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1dp1;
        ddr1dtdr2dr1 = drnp1dunp1 * ddu1dtdr2dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1dr2;
        ddr1dtdu2dr1 = drnp1dunp1 * ddu1dtdu2dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1du2;
        ddr1dtdp2dr1 = drnp1dunp1 * ddu1dtdp2dr1
                       + ddrnp1dunp1dunp1 * dunp1dr1 * dunp1dp2;

        ddr1dtdr1du1 = ddr1dtdu1dr1;
        ddr1dtdu1du1 = drnp1dunp1 * ddu1dtdu1du1
                       + ddrnp1dunp1dunp1 * dunp1du1 * dunp1du1;
        ddr1dtdp1du1 = drnp1dunp1 * ddu1dtdp1du1
                       + ddrnp1dunp1dunp1 * dunp1du1 * dunp1dp1;
        ddr1dtdr2du1 = drnp1dunp1 * ddu1dtdr2du1
                       + ddrnp1dunp1dunp1 * dunp1du1 * dunp1dr2;
        ddr1dtdu2du1 = drnp1dunp1 * ddu1dtdu2du1
                       + ddrnp1dunp1dunp1 * dunp1du1 * dunp1du2;
        ddr1dtdp2du1 = drnp1dunp1 * ddu1dtdp2du1
                       + ddrnp1dunp1dunp1 * dunp1du1 * dunp1dp2;

        ddr1dtdr1dp1 = ddr1dtdp1dr1;
        ddr1dtdu1dp1 = ddr1dtdp1du1;
        ddr1dtdp1dp1 = drnp1dunp1 * ddu1dtdp1dp1
                       + ddrnp1dunp1dunp1 * dunp1dp1 * dunp1dp1;
        ddr1dtdr2dp1 = drnp1dunp1 * ddu1dtdr2dp1
                       + ddrnp1dunp1dunp1 * dunp1dp1 * dunp1dr2;
        ddr1dtdu2dp1 = drnp1dunp1 * ddu1dtdu2dp1
                       + ddrnp1dunp1dunp1 * dunp1dp1 * dunp1du2;
        ddr1dtdp2dp1 = drnp1dunp1 * ddu1dtdp2dp1
                       + ddrnp1dunp1dunp1 * dunp1dp1 * dunp1dp2;

        ddr1dtdr1dr2 = ddr1dtdr2dr1;
        ddr1dtdu1dr2 = ddr1dtdr2du1;
        ddr1dtdp1dr2 = ddr1dtdr2dp1;
        ddr1dtdr2dr2 = drnp1dunp1 * ddu1dtdr2dr2
                       + ddrnp1dunp1dunp1 * dunp1dr2 * dunp1dr2;
        ddr1dtdu2dr2 = drnp1dunp1 * ddu1dtdu2dr2
                       + ddrnp1dunp1dunp1 * dunp1dr2 * dunp1du2;
        ddr1dtdp2dr2 = drnp1dunp1 * ddu1dtdp2dr2
                       + ddrnp1dunp1dunp1 * dunp1dr2 * dunp1dp2;

        ddr1dtdr1du2 = ddr1dtdu2dr1;
        ddr1dtdu1du2 = ddr1dtdu2du1;
        ddr1dtdp1du2 = ddr1dtdu2dp1;
        ddr1dtdr2du2 = ddr1dtdu2dr2;
        ddr1dtdu2du2 = drnp1dunp1 * ddu1dtdu2du2
                       + ddrnp1dunp1dunp1 * dunp1du2 * dunp1du2;
        ddr1dtdp2du2 = drnp1dunp1 * ddu1dtdp2du2
                       + ddrnp1dunp1dunp1 * dunp1du2 * dunp1dp2;

        ddr1dtdr1dp2 = ddr1dtdp2dr1;
        ddr1dtdu1dp2 = ddr1dtdp2du1;
        ddr1dtdp1dp2 = ddr1dtdp2dp1;
        ddr1dtdr2dp2 = ddr1dtdp2dr2;
        ddr1dtdu2dp2 = ddr1dtdp2du2;
        ddr1dtdp2dp2 = drnp1dunp1 * ddu1dtdp2dp2
                       + ddrnp1dunp1dunp1 * dunp1dp2 * dunp1dp2;

        // dru1/dt
//      dru1dt = r1 * du1dt + u1 * dr1dt;

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

        double
        ddru1dtdr1dr1, ddru1dtdu1dr1, ddru1dtdp1dr1, ddru1dtdr2dr1, ddru1dtdu2dr1, ddru1dtdp2dr1,
        ddru1dtdr1du1, ddru1dtdu1du1, ddru1dtdp1du1, ddru1dtdr2du1, ddru1dtdu2du1, ddru1dtdp2du1,
        ddru1dtdr1dp1, ddru1dtdu1dp1, ddru1dtdp1dp1, ddru1dtdr2dp1, ddru1dtdu2dp1, ddru1dtdp2dp1,
        ddru1dtdr1dr2, ddru1dtdu1dr2, ddru1dtdp1dr2, ddru1dtdr2dr2, ddru1dtdu2dr2, ddru1dtdp2dr2,
        ddru1dtdr1du2, ddru1dtdu1du2, ddru1dtdp1du2, ddru1dtdr2du2, ddru1dtdu2du2, ddru1dtdp2du2,
        ddru1dtdr1dp2, ddru1dtdu1dp2, ddru1dtdp1dp2, ddru1dtdr2dp2, ddru1dtdu2dp2, ddru1dtdp2dp2;

        ddru1dtdr1dr1 = u1 * ddr1dtdr1dr1 + 2.0 * du1dtdr1 + r1 * ddu1dtdr1dr1
                        + ddr1dtdr1dr1 * du1dt + dr1dtdr1 * du1dtdr1
                        + dr1dt * ddu1dtdr1dr1 + dr1dtdr1 * du1dtdr1;
        ddru1dtdu1dr1 = u1 * ddr1dtdr1du1 + dr1dtdr1 + du1dtdu1 + r1 * ddu1dtdr1du1
                        + ddr1dtdr1du1 * du1dt + dr1dtdr1 * du1dtdu1
                        + dr1dt * ddu1dtdr1du1 + dr1dtdu1 * du1dtdr1;
        ddru1dtdp1dr1 = u1 * ddr1dtdp1dr1 + du1dtdp1 + r1 * ddu1dtdp1dr1
                        + ddr1dtdr1dp1 * du1dt + dr1dtdr1 * du1dtdp1
                        + dr1dt * ddu1dtdr1dp1 + dr1dtdp1 * du1dtdr1;
        ddru1dtdr2dr1 = u1 * ddr1dtdr1dr2 + du1dtdr2 + r1 * ddu1dtdr1dr2
                        + ddr1dtdr2dr1 * du1dt + dr1dtdr2 * du1dtdr1
                        + dr1dt * ddu1dtdr2dr1 + dr1dtdr1 * du1dtdr2;
        ddru1dtdu2dr1 = u1 * ddr1dtdu2dr1 + du1dtdu2 + r1 * ddu1dtdu2dr1
                        + ddr1dtdu2dr1 * du1dt + dr1dtdu2 * du1dtdr1
                        + dr1dt * ddu1dtdu2dr1 + dr1dtdr1 * du1dtdu2;
        ddru1dtdp2dr1 = u1 * ddr1dtdp2dr1 + du1dtdp2 + r1 * ddu1dtdp2dr1
                        + ddr1dtdp2dr1 * du1dt + dr1dtdp2 * du1dtdr1
                        + dr1dt * ddu1dtdp2dr1 + dr1dtdr1 * du1dtdp2;

        ddru1dtdr1du1 = ddru1dtdu1dr1;
        ddru1dtdu1du1 = u1 * ddr1dtdu1du1 + 2.0 * dr1dtdu1 + r1 * ddu1dtdu1du1
                        + ddr1dtdu1du1 * du1dt + dr1dtdu1 * du1dtdu1
                        + dr1dt * ddu1dtdu1du1 + dr1dtdu1 * du1dtdu1;
        ddru1dtdp1du1 = u1 * ddr1dtdp1du1 + dr1dtdp1 + r1 * ddu1dtdp1du1
                        + ddr1dtdu1dp1 * du1dt + dr1dtdu1 * du1dtdp1
                        + dr1dt * ddu1dtdu1dp1 + dr1dtdp1 * du1dtdu1;
        ddru1dtdr2du1 = u1 * ddr1dtdr2du1 + dr1dtdr2 + r1 * ddu1dtdr2du1
                        + ddr1dtdr2du1 * du1dt + dr1dtdr2 * du1dtdu1
                        + dr1dt * ddu1dtdr2du1 + dr1dtdu1 * du1dtdr2;
        ddru1dtdu2du1 = u1 * ddr1dtdu2du1 + dr1dtdu2 + r1 * ddu1dtdu2du1
                        + ddr1dtdu2du1 * du1dt + dr1dtdu2 * du1dtdu1
                        + dr1dt * ddu1dtdu2du1 + dr1dtdu1 * du1dtdu2;
        ddru1dtdp2du1 = u1 * ddr1dtdp2du1 + dr1dtdp2 + r1 * ddu1dtdp2du1
                        + ddr1dtdp2du1 * du1dt + dr1dtdp2 * du1dtdu1
                        + dr1dt * ddu1dtdp2du1 + dr1dtdu1 * du1dtdp2;


        ddru1dtdr1dp1 = ddru1dtdp1dr1;
        ddru1dtdu1dp1 = ddru1dtdp1du1;
        ddru1dtdp1dp1 = u1 * ddr1dtdp1dp1 + r1 * ddu1dtdp1dp1
                        + ddr1dtdp1dp1 * du1dt + dr1dtdp1 * du1dtdp1
                        + dr1dt * ddu1dtdp1dp1 + dr1dtdp1 * du1dtdp1;
        ddru1dtdr2dp1 = u1 * ddr1dtdr2dp1 + r1 * ddu1dtdr2dp1
                        + ddr1dtdr2dp1 * du1dt + dr1dtdr2 * du1dtdp1
                        + dr1dt * ddu1dtdr2dp1 + dr1dtdp1 * du1dtdr2;
        ddru1dtdu2dp1 = u1 * ddr1dtdu2dp1 + r1 * ddu1dtdu2dp1
                        + ddr1dtdu2dp1 * du1dt + dr1dtdu2 * du1dtdp1
                        + dr1dt * ddu1dtdu2dp1 + dr1dtdp1 * du1dtdu2;
        ddru1dtdp2dp1 = u1 * ddr1dtdp2dp1 + r1 * ddu1dtdp2dp1
                        + ddr1dtdp2dp1 * du1dt + dr1dtdp2 * du1dtdp1
                        + dr1dt * ddu1dtdp2dp1 + dr1dtdp1 * du1dtdp2;
        
        ddru1dtdr1dr2 = ddru1dtdr2dr1;
        ddru1dtdu1dr2 = ddru1dtdr2du1;
        ddru1dtdp1dr2 = ddru1dtdr2dp1;
        ddru1dtdr2dr2 = u1 * ddr1dtdr2dr2 + r1 * ddu1dtdr2dr2
                        + ddr1dtdr2dr2 * du1dt + dr1dtdr2 * du1dtdr2
                        + dr1dt * ddu1dtdr2dr2 + dr1dtdr2 * du1dtdr2;
        ddru1dtdu2dr2 = u1 * ddr1dtdr2du2 + r1 * ddu1dtdr2du2
                        + ddr1dtdr2du2 * du1dt + dr1dtdr2 * du1dtdu2
                        + dr1dt * ddu1dtdr2du2 + dr1dtdu2 * du1dtdr2;
        ddru1dtdp2dr2 = u1 * ddr1dtdr2dp2 + r1 * ddu1dtdr2dp2
                        + ddr1dtdr2dp2 * du1dt + dr1dtdr2 * du1dtdp2
                        + dr1dt * ddu1dtdr2dp2 + dr1dtdp2 * du1dtdr2;


        ddru1dtdr1du2 = ddru1dtdu2dr1;
        ddru1dtdu1du2 = ddru1dtdu2du1;
        ddru1dtdp1du2 = ddru1dtdu2dp1;
        ddru1dtdr2du2 = ddru1dtdu2dr2;
        ddru1dtdu2du2 = u1 * ddr1dtdu2du2 + r1 * ddu1dtdu2du2
                        + ddr1dtdu2du2 * du1dt + dr1dtdu2 * du1dtdu2
                        + dr1dt * ddu1dtdu2du2 + dr1dtdu2 * du1dtdu2;
        ddru1dtdp2du2 = u1 * ddr1dtdu2dp2 + r1 * ddu1dtdu2dp2
                        + ddr1dtdu2dp2 * du1dt + dr1dtdu2 * du1dtdp2
                        + dr1dt * ddu1dtdu2dp2 + dr1dtdp2 * du1dtdu2;

        ddru1dtdr1dp2 = ddru1dtdp2dr1;
        ddru1dtdu1dp2 = ddru1dtdp2du1;
        ddru1dtdp1dp2 = ddru1dtdp2dp1;
        ddru1dtdr2dp2 = ddru1dtdp2dr2;        
        ddru1dtdu2dp2 = ddru1dtdp2du2;
        ddru1dtdp2dp2 = u1 * ddr1dtdp2dp2 + r1 * ddu1dtdp2dp2
                        + ddr1dtdp2dp2 * du1dt + dr1dtdp2 * du1dtdp2
                        + dr1dt * ddu1dtdp2dp2 + dr1dtdp2 * du1dtdp2;

        // de1/dt
//      de1dt = dp1dt / (gam - 1.0) + r1 * u1 * du1dt + u1 * u1 * dr1dt / 2.0;

        // First Derivative
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
//
// ***********************************************************************
//
        dRidWi(0,0) = dr1dtdr1;
        dRidWi(0,1) = dr1dtdu1;
        dRidWi(0,2) = dr1dtdp1;
        dRidWi(1,0) = dru1dtdr1;
        dRidWi(1,1) = dru1dtdu1;
        dRidWi(1,2) = dru1dtdp1;
        dRidWi(2,0) = de1dtdr1;
        dRidWi(2,1) = de1dtdu1;
        dRidWi(2,2) = de1dtdp1;

        dRidWd(0,0) = dr1dtdr2;
        dRidWd(0,1) = dr1dtdu2;
        dRidWd(0,2) = dr1dtdp2;
        dRidWd(1,0) = dru1dtdr2;
        dRidWd(1,1) = dru1dtdu2;
        dRidWd(1,2) = dru1dtdp2;
        dRidWd(2,0) = de1dtdr2;
        dRidWd(2,1) = de1dtdu2;
        dRidWd(2,2) = de1dtdp2;

        // ****************
        // 1st Residual
        // ****************
        // Top Left
        ddRindWdW[0](0, 0) = ddr1dtdr1dr1;
        ddRindWdW[0](1, 1) = ddr1dtdu1du1;
        ddRindWdW[0](2, 2) = ddr1dtdp1dp1;

        ddRindWdW[0](0, 1) = ddr1dtdr1du1;
        ddRindWdW[0](1, 0) = ddr1dtdu1dr1;

        ddRindWdW[0](0, 2) = ddr1dtdr1dp1;
        ddRindWdW[0](2, 0) = ddr1dtdp1dr1;

        ddRindWdW[0](1, 2) = ddr1dtdu1dp1;
        ddRindWdW[0](2, 1) = ddr1dtdp1du1;

        // Bottom Right
        ddRindWdW[0](3, 3) = ddr1dtdr2dr2;
        ddRindWdW[0](4, 4) = ddr1dtdu2du2;
        ddRindWdW[0](5, 5) = ddr1dtdp2dp2;

        ddRindWdW[0](3, 4) = ddr1dtdr2du2;
        ddRindWdW[0](4, 3) = ddr1dtdu2dr2;

        ddRindWdW[0](3, 5) = ddr1dtdr2dp2;
        ddRindWdW[0](5, 3) = ddr1dtdp2dr2;

        ddRindWdW[0](4, 5) = ddr1dtdu2dp2;
        ddRindWdW[0](5, 4) = ddr1dtdp2du2;

        // Top Right
        ddRindWdW[0](0, 3) = ddr1dtdr1dr2;
        ddRindWdW[0](1, 4) = ddr1dtdu1du2;
        ddRindWdW[0](2, 5) = ddr1dtdp1dp2;

        ddRindWdW[0](0, 4) = ddr1dtdr1du2;
        ddRindWdW[0](1, 3) = ddr1dtdu1dr2;

        ddRindWdW[0](0, 5) = ddr1dtdr1dp2;
        ddRindWdW[0](2, 3) = ddr1dtdp1dr2;

        ddRindWdW[0](1, 5) = ddr1dtdu1dp2;
        ddRindWdW[0](2, 4) = ddr1dtdp1du2;


        // Bottom Left
        ddRindWdW[0](3, 0) = ddr1dtdr2dr1;
        ddRindWdW[0](4, 1) = ddr1dtdu2du1;
        ddRindWdW[0](5, 2) = ddr1dtdp2dp1;

        ddRindWdW[0](4, 0) = ddr1dtdu2dr1;
        ddRindWdW[0](3, 1) = ddr1dtdr2du1;

        ddRindWdW[0](5, 0) = ddr1dtdp2dr1;
        ddRindWdW[0](3, 2) = ddr1dtdr2dp1;

        ddRindWdW[0](5, 1) = ddr1dtdp2du1;
        ddRindWdW[0](4, 2) = ddr1dtdu2dp1;

        // ****************
        // 2nd Residual
        // ****************
        // Top Left
        ddRindWdW[1](0, 0) = ddru1dtdr1dr1;
        ddRindWdW[1](1, 1) = ddru1dtdu1du1;
        ddRindWdW[1](2, 2) = ddru1dtdp1dp1;

        ddRindWdW[1](0, 1) = ddru1dtdr1du1;
        ddRindWdW[1](1, 0) = ddru1dtdu1dr1;

        ddRindWdW[1](0, 2) = ddru1dtdr1dp1;
        ddRindWdW[1](2, 0) = ddru1dtdp1dr1;

        ddRindWdW[1](1, 2) = ddru1dtdu1dp1;
        ddRindWdW[1](2, 1) = ddru1dtdp1du1;

        // Bottom Right
        ddRindWdW[1](3, 3) = ddru1dtdr2dr2;
        ddRindWdW[1](4, 4) = ddru1dtdu2du2;
        ddRindWdW[1](5, 5) = ddru1dtdp2dp2;

        ddRindWdW[1](3, 4) = ddru1dtdr2du2;
        ddRindWdW[1](4, 3) = ddru1dtdu2dr2;

        ddRindWdW[1](3, 5) = ddru1dtdr2dp2;
        ddRindWdW[1](5, 3) = ddru1dtdp2dr2;

        ddRindWdW[1](4, 5) = ddru1dtdu2dp2;
        ddRindWdW[1](5, 4) = ddru1dtdp2du2;

        // Top Right
        ddRindWdW[1](0, 3) = ddru1dtdr1dr2;
        ddRindWdW[1](1, 4) = ddru1dtdu1du2;
        ddRindWdW[1](2, 5) = ddru1dtdp1dp2;

        ddRindWdW[1](0, 4) = ddru1dtdr1du2;
        ddRindWdW[1](1, 3) = ddru1dtdu1dr2;

        ddRindWdW[1](0, 5) = ddru1dtdr1dp2;
        ddRindWdW[1](2, 3) = ddru1dtdp1dr2;

        ddRindWdW[1](1, 5) = ddru1dtdu1dp2;
        ddRindWdW[1](2, 4) = ddru1dtdp1du2;


        // Bottom Left
        ddRindWdW[1](3, 0) = ddru1dtdr2dr1;
        ddRindWdW[1](4, 1) = ddru1dtdu2du1;
        ddRindWdW[1](5, 2) = ddru1dtdp2dp1;

        ddRindWdW[1](4, 0) = ddru1dtdu2dr1;
        ddRindWdW[1](3, 1) = ddru1dtdr2du1;

        ddRindWdW[1](5, 0) = ddru1dtdp2dr1;
        ddRindWdW[1](3, 2) = ddru1dtdr2dp1;

        ddRindWdW[1](5, 1) = ddru1dtdp2du1;
        ddRindWdW[1](4, 2) = ddru1dtdu2dp1;

        // ****************
        // 3rd Residual
        // ****************
        // Top Left
        ddRindWdW[2](0, 0) = dde1dtdr1dr1;
        ddRindWdW[2](1, 1) = dde1dtdu1du1;
        ddRindWdW[2](2, 2) = dde1dtdp1dp1;

        ddRindWdW[2](0, 1) = dde1dtdr1du1;
        ddRindWdW[2](1, 0) = dde1dtdu1dr1;

        ddRindWdW[2](0, 2) = dde1dtdr1dp1;
        ddRindWdW[2](2, 0) = dde1dtdp1dr1;

        ddRindWdW[2](1, 2) = dde1dtdu1dp1;
        ddRindWdW[2](2, 1) = dde1dtdp1du1;

        // Bottom Right
        ddRindWdW[2](3, 3) = dde1dtdr2dr2;
        ddRindWdW[2](4, 4) = dde1dtdu2du2;
        ddRindWdW[2](5, 5) = dde1dtdp2dp2;

        ddRindWdW[2](3, 4) = dde1dtdr2du2;
        ddRindWdW[2](4, 3) = dde1dtdu2dr2;

        ddRindWdW[2](3, 5) = dde1dtdr2dp2;
        ddRindWdW[2](5, 3) = dde1dtdp2dr2;

        ddRindWdW[2](4, 5) = dde1dtdu2dp2;
        ddRindWdW[2](5, 4) = dde1dtdp2du2;

        // Top Right
        ddRindWdW[2](0, 3) = dde1dtdr1dr2;
        ddRindWdW[2](1, 4) = dde1dtdu1du2;
        ddRindWdW[2](2, 5) = dde1dtdp1dp2;

        ddRindWdW[2](0, 4) = dde1dtdr1du2;
        ddRindWdW[2](1, 3) = dde1dtdu1dr2;

        ddRindWdW[2](0, 5) = dde1dtdr1dp2;
        ddRindWdW[2](2, 3) = dde1dtdp1dr2;

        ddRindWdW[2](1, 5) = dde1dtdu1dp2;
        ddRindWdW[2](2, 4) = dde1dtdp1du2;


        // Bottom Left
        ddRindWdW[2](3, 0) = dde1dtdr2dr1;
        ddRindWdW[2](4, 1) = dde1dtdu2du1;
        ddRindWdW[2](5, 2) = dde1dtdp2dp1;

        ddRindWdW[2](4, 0) = dde1dtdu2dr1;
        ddRindWdW[2](3, 1) = dde1dtdr2du1;

        ddRindWdW[2](5, 0) = dde1dtdp2dr1;
        ddRindWdW[2](3, 2) = dde1dtdr2dp1;

        ddRindWdW[2](5, 1) = dde1dtdp2du1;
        ddRindWdW[2](4, 2) = dde1dtdu2dp1;

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
                temp += dRidWi(Ri, Wpi) * ddwpdwdwp[Wpi];
            }
            ddRindWdW[Ri].topLeftCorner(3, 3) = temp * dwpdw;
        }

        // Get Transformation Matrices
        dwpdw = dWpdW(W, 1);
        ddwpdwdwp = ddWpdWdWp(W, 1);

        for(int Ri = 0; Ri < 3; Ri++)
        {
            temp.setZero();
            temp += dwpdw.transpose() * ddRindWdW[Ri].bottomRightCorner(3, 3);
            for(int Wpi = 0; Wpi < 3; Wpi++)
            {
                temp += dRidWd(Ri, Wpi) * ddwpdwdwp[Wpi];
            }
            ddRindWdW[Ri].bottomRightCorner(3, 3) = temp * dwpdw;
        }

        // Get Transformation Matrices
        dwpdw = dWpdW(W, 0);
        MatrixXd dwpdw2 = dWpdW(W, 1);

        for(int Ri = 0; Ri < 3; Ri++)
        {
            ddRindWdW[Ri].topRightCorner(3, 3) =
                dwpdw.transpose() * ddRindWdW[Ri].topRightCorner(3, 3) * dwpdw2;
        }
        for(int Ri = 0; Ri < 3; Ri++)
        {
            ddRindWdW[Ri].bottomLeftCorner(3, 3) =
                dwpdw2.transpose() * ddRindWdW[Ri].bottomLeftCorner(3, 3) * dwpdw;
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
