#include <math.h>
#include <iostream>
#include <vector>
#include "structures.h"
#include "convert.h"
// Jacobian at the Boundaries
void dRdW_BC_inlet(
	const struct Flow_options &flow_options,
	const std::vector<double> &W,
    std::vector<double> &dBidWi,
    std::vector<double> &dBidWd)
{
	const double gam = flow_options.gam;
	const double Cv  = flow_options.Cv;
	const double R   = flow_options.R;
	const double a2  = flow_options.a2;
	const double inlet_total_p  = flow_options.inlet_total_p;
	const double inlet_total_T  = flow_options.inlet_total_T;
    std::vector<double> dbdwp(9, 0), dwpdw(9);

    for (int i = 0; i < 9; i++) {
        dBidWi[i] = 0;
        dBidWd[i] = 0;
    }

    // *********************
    // INLET JACOBIANS
    // *********************
    // Subsonic Inlet
	const int i1 = 0;
	const int i2 = 1;
	const double r1 = W[i1*3+0];
	const double r2 = W[i2*3+0];
	const double p1 = get_p(gam, W[i1*3+0], W[i1*3+1], W[i1*3+2]);
	const double p2 = get_p(gam, W[i2*3+0], W[i2*3+1], W[i2*3+2]);
	const double u1 = W[i1*3+1] / r1;
	const double u2 = W[i2*3+1] / r2;
	const double c1 = get_c(gam, W[i1*3+0], W[i1*3+1], W[i1*3+2]);
	const double c2 = get_c(gam, W[i2*3+0], W[i2*3+1], W[i2*3+2]);
    if (u1 < c1) {
        // Shorthand
        const double gamr = (gam - 1.0) / (gam + 1.0);
        const double fu = 1.0 - gamr * u1 * u1 / a2;

        // Speed of Sound
        const double dc1dr1 = - p1 * gam / (2.0 * r1 * c1 * r1);
        const double dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
        const double dc1dp1 = gam / (2.0 * r1 * c1);
        const double dc2dp2 = gam / (2.0 * c2 * r2);

        // Eigenvalue
        const double eig1 = (u1 + u2) / 2.0;
        const double eig3 = eig1 - (c1 + c2) / 2.0;

        const double deig1du1 = 0.5;
        const double deig1du2 = 0.5;

        const double deig3dr1 = - dc1dr1 / 2.0;
        const double deig3du1 = deig1du1;
        const double deig3dp1 = - dc1dp1 / 2.0;
        const double deig3dr2 = - dc2dr2 / 2.0;
        const double deig3du2 = deig1du2;
        const double deig3dp2 = - dc2dp2 / 2.0;

        // Riemann Invariants
        const double R3 = - eig3 * ((p2 - p1) - r1 * c1 * (u2 - u1));

        const double dR3dr1 = -eig3 * (-c1 * (u2 - u1) - (u2 - u1) * r1 * dc1dr1)
                              - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dr1;
        const double dR3du1 = -r1 * c1 * eig3
                              - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3du1;
        const double dR3dp1 = eig3 * (1.0 + (u2 - u1) * r1 * dc1dp1)
                              - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dp1;
        const double dR3dr2 = - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dr2;
        const double dR3du2 = r1 * c1 * eig3
                              - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3du2;
        const double dR3dp2 = - eig3
                 - ((p2 - p1) - r1 * c1 * (u2 - u1)) * deig3dp2;
        // dp1
        // Same Values
        const double dp1du1 = -2.0 * gamr * inlet_total_p * u1 * pow(fu, 1.0 / (gam - 1.0)) * gam
                 / (a2 * (gam - 1.0));

        const double dp1du1du1 = 2.0 * gamr * inlet_total_p * pow(fu, gam / (gam - 1.0)) * gam
                    * (a2 - a2 * gam + gamr * u1 * u1 * (gam + 1))
                    / pow((a2 - gamr * u1 * u1) * (gam - 1.0), 2);

        // du1
        const double du1dt = R3 / (dp1du1 - r1 * c1);
        const double du1dtdr1 = dR3dr1 / (dp1du1 - r1 * c1)
                   - R3 * (-c1 - r1 * dc1dr1) / pow((dp1du1 - r1 * c1), 2);
        const double du1dtdu1 = dR3du1 / (dp1du1 - r1 * c1)
                   - R3 * dp1du1du1 / pow((dp1du1 - r1 * c1), 2);
        const double du1dtdp1 = dR3dp1 / (dp1du1 - r1 * c1)
                   + (R3 * r1 * dc1dp1) / pow((dp1du1 - r1 * c1), 2);
        const double du1dtdr2 = dR3dr2 / (dp1du1 - r1 * c1);
        const double du1dtdu2 = dR3du2 / (dp1du1 - r1 * c1);
        const double du1dtdp2 = dR3dp2 / (dp1du1 - r1 * c1);

        // Primitive values at time-step n+1
        const double unp1 = u1 + du1dt;
        const double pnp1 = inlet_total_p * pow(1 - gamr * pow(unp1, 2) / a2, gam / (gam - 1.0));
        const double tnp1 = inlet_total_T * ( 1 - gamr * unp1 * unp1 / a2 );
        const double rnp1 = pnp1 / (R * tnp1);
        const double dpnp1dunp1 = -2.0 * gamr * inlet_total_p * unp1
                     * pow((1.0 - gamr * unp1 * unp1 / a2), 1.0 / (gam - 1.0))
                     * gam / (a2 * (gam - 1.0));

        const double dunp1dr1 = du1dtdr1;
        const double dunp1du1 = 1.0 + du1dtdu1;
        const double dunp1dp1 = du1dtdp1;
        const double dunp1dr2 = du1dtdr2;
        const double dunp1du2 = du1dtdu2;
        const double dunp1dp2 = du1dtdp2;

        // dp1
        //const double dp1dt  = pnp1 - p1;
        const double dp1dtdr1 = dpnp1dunp1 * dunp1dr1;
        const double dp1dtdu1 = dpnp1dunp1 * dunp1du1;
        const double dp1dtdp1 = dpnp1dunp1 * dunp1dp1 - 1.0; // -1.0 due to dp1dp1
        const double dp1dtdr2 = dpnp1dunp1 * dunp1dr2;
        const double dp1dtdu2 = dpnp1dunp1 * dunp1du2;
        const double dp1dtdp2 = dpnp1dunp1 * dunp1dp2;

        // dr1
        // Total derivative from rho_n+1 to p_n+1 and u_n+1
        //const double drnp1dpnp1 = 1.0 / (R * tnp1);
        //const double drnp1dtnp1 = -pnp1 / (R * tnp1 * tnp1);
        //const double dtnp1dpnp1 = inlet_total_T / inlet_total_p * (gam - 1.0) / gam * pow(pnp1 / inlet_total_p, - 1.0 / gam);
        //double Drnp1Dpnp1 = drnp1dpnp1 + drnp1dtnp1 * dtnp1dpnp1;
        //double drnp1dunp1 = Drnp1Dpnp1 * dpnp1dunp1;
        const double drnp1dunp1 = (-2.0 * a2 * (1.0 + gam) * inlet_total_p *unp1 
                     * pow(1.0 + (pow(unp1, 2) 
                     - gam * pow(unp1, 2)) / (a2 + a2*gam),
                     gam/(-1 + gam)))
                     / (R * inlet_total_T * pow(a2 + a2 * gam + pow(unp1, 2) 
                     - gam * pow(unp1, 2), 2));

        const double dr1dt = rnp1 - r1;

        const double dr1dtdr1 = drnp1dunp1 * du1dtdr1 - 1.0;
        const double dr1dtdu1 = drnp1dunp1 * (1.0 + du1dtdu1);
        const double dr1dtdp1 = drnp1dunp1 * du1dtdp1;

        const double dr1dtdr2 = drnp1dunp1 * du1dtdr2;
        const double dr1dtdu2 = drnp1dunp1 * du1dtdu2;
        const double dr1dtdp2 = drnp1dunp1 * du1dtdp2;

        // dru1/dt
//      dru1dt = r1 * du1dt + u1 * dr1dt + dr1dt * du1dt;

        const double dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1
                    + dr1dtdr1 * du1dt + dr1dt * du1dtdr1;
        const double dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1
                    + dr1dtdu1 * du1dt + dr1dt * du1dtdu1;
        const double dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1
                    + dr1dtdp1 * du1dt + dr1dt * du1dtdp1;
        const double dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2
                    + dr1dtdr2 * du1dt + dr1dt * du1dtdr2;
        const double dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2
                    + dr1dtdu2 * du1dt + dr1dt * du1dtdu2;
        const double dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2
                    + dr1dtdp2 * du1dt + dr1dt * du1dtdp2;
        // de1/dt

        const double de1dtdr1 = (2.0 * Cv * dp1dtdr1
                   + R * (du1dt * (du1dt + 2.0 * u1)
                   + pow(du1dt + u1, 2.0) * dr1dtdr1
                   + 2.0 * (dr1dt + r1) * (du1dt + u1) * du1dtdr1))
                   / (2.0 * R);
        const double de1dtdu1 = (2.0 * Cv * dp1dtdu1
                   + R * (pow(du1dt + u1, 2.0) * dr1dtdu1
                   + 2.0 * (du1dt * r1 + dr1dt * (du1dt + u1)
                   + (dr1dt + r1) * (du1dt + u1) * du1dtdu1)))
                   / (2.0 * R);
        const double de1dtdp1 = (2.0 * Cv * dp1dtdp1
                   + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp1
                   + 2.0 * (dr1dt + r1) * du1dtdp1))
                   / (2.0 * R);
        const double de1dtdr2 = (2.0 * Cv * dp1dtdr2
                   + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdr2
                   + 2.0 * (dr1dt + r1) * du1dtdr2))
                   / (2.0 * R);
        const double de1dtdu2 = (2.0 * Cv * dp1dtdu2
                   + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdu2
                   + 2.0 * (dr1dt + r1) * du1dtdu2))
                   / (2.0 * R);
        const double de1dtdp2 = (2.0 * Cv * dp1dtdp2
                   + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp2
                   + 2.0 * (dr1dt + r1) * du1dtdp2))
                   / (2.0 * R);
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
        eval_dWpdW(gam, W[i1*3+0], W[i1*3+1], &dwpdw);

        for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
        for (int k = 0; k < 3; k++)
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
        eval_dWpdW(gam, W[i2*3+0], W[i2*3+1], &dwpdw);

        for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
        for (int k = 0; k < 3; k++)
            dBidWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
    } else { // Supersonic Inlet
		for (int i = 0; i < 9; i++) {
			dBidWi[i] = 0;
			dBidWd[i] = 0;
			if (i % 4 == 0)
				dBidWi[i] = 1;
		}
    }
}

void dRdW_BC_outlet(
	const struct Flow_options &flow_options,
	const std::vector<double> &W,
    std::vector<double> &dBodWd,
    std::vector<double> &dBodWo)
{
	int n_elem = flow_options.n_elem;
	double gam = flow_options.gam;
	double Cv  = flow_options.Cv;
	double R   = flow_options.R;
    std::vector<double> dbdwp(9, 0), dwpdw(9);

    for (int i = 0; i < 9; i++) {
        dBodWd[i] = 0;
        dBodWo[i] = 0;
    }

    // ************************
    // OUTLET JACOBIANS
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

    // Speed of Sound
    const double dc1dr1 = - p1 * gam / (2.0 * r1 * c1 * r1);
    const double dc2dr2 = - p2 * gam / (2.0 * c2 * r2 * r2);
    const double dc1dp1 = gam / (2.0 * r1 * c1);
    const double dc2dp2 = gam / (2.0 * c2 * r2);

    // Eigenvalue
    const double eig1 = (u1 + u2) / 2.0;
    const double eig2 = eig1 + (c1 + c2) / 2.0;
    const double eig3 = eig1 - (c1 + c2) / 2.0;

    const double deig1du1 = 0.5;
    const double deig1du2 = 0.5;

    const double deig2dr1 = dc1dr1 / 2.0;
    const double deig2du1 = deig1du1;
    const double deig2dp1 = dc1dp1 / 2.0;
    const double deig2dr2 = dc2dr2 / 2.0;
    const double deig2du2 = deig1du2;
    const double deig2dp2 = dc2dp2 / 2.0;

    const double deig3dr1 = - dc1dr1 / 2.0;
    const double deig3du1 = deig1du1;
    const double deig3dp1 = - dc1dp1 / 2.0;
    const double deig3dr2 = - dc2dr2 / 2.0;
    const double deig3du2 = deig1du2;
    const double deig3dp2 = - dc2dp2 / 2.0;

    // Riemann invariants
    const double R1 = - eig1 * ((r1 - r2) - (p1 - p2) / (c1 * c1));
    const double R2 = - eig2 * ((p1 - p2) + r1 * c1 * (u1 - u2));
    const double R3 = - eig3 * ((p1 - p2) - r1 * c1 * (u1 - u2));

    const double dR1dr1 = - eig1 * (1.0 + 2.0 * (p1 - p2) * dc1dr1 / pow(c1, 3) );
    const double dR1du1 = deig1du1 * ((p1 - p2) - c1 * c1 * (r1 - r2)) / (c1 * c1);
    const double dR1dp1 = eig1 * (c1 - 2.0 * (p1 - p2) * dc1dp1) / pow(c1, 3);
    const double dR1dr2 = eig1;
    const double dR1du2 = deig1du2 * ((p1 - p2) - c1 * c1 * (r1 - r2)) / (c1 * c1);
    const double dR1dp2 = - eig1 / (c1 * c1);

    const double dR2dr1 = - (u1 - u2) * eig2 * (c1 + r1 * dc1dr1)
                          - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dr1;
    const double dR2du1 = - r1 * c1 * eig2
                          - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2du1;
    const double dR2dp1 = - eig2 * (1.0 + (u1 - u2) * r1 * dc1dp1)
                          - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dp1;
    const double dR2dr2 = - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dr2;
    const double dR2du2 = r1 * c1 * eig2
                          - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2du2;
    const double dR2dp2 = eig2
                          - ((p1 - p2) + r1 * c1 * (u1 - u2)) * deig2dp2;

    const double dR3dr1 = eig3 * (u1 - u2) * (c1 + r1 * dc1dr1)
                          - ((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dr1;
    const double dR3du1 = r1 * c1 * eig3
                          - ((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3du1;
    const double dR3dp1 = - eig3 - (p1 - p2) * deig3dp1
                          + (u1 - u2) * r1 * eig3 * dc1dp1;
    const double dR3dr2 = - ((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dr2;
    const double dR3du2 = - r1 * c1 * eig3
                          - ((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3du2;
    const double dR3dp2 = eig3
                          - ((p1 - p2) - r1 * c1 * (u1 - u2)) * deig3dp2;

    // dp1/dt
    double dp1dt;
    double dp1dtdr1, dp1dtdu1, dp1dtdp1;
    double dp1dtdr2, dp1dtdu2, dp1dtdp2;
    if (u1 < c1) {
        dp1dt = 0;
        dp1dtdr1 = 0;
        dp1dtdu1 = 0;
        dp1dtdp1 = 1;
        dp1dtdr2 = 0;
        dp1dtdu2 = 0;
        dp1dtdp2 = 0;
    } else {
        dp1dt = (R2 + R3) / 2.0;
        dp1dtdr1 = (dR2dr1 + dR3dr1) / 2.0;
        dp1dtdu1 = (dR2du1 + dR3du1) / 2.0;
        dp1dtdp1 = (dR2dp1 + dR3dp1) / 2.0;
        dp1dtdr2 = (dR2dr2 + dR3dr2) / 2.0;
        dp1dtdu2 = (dR2du2 + dR3du2) / 2.0;
        dp1dtdp2 = (dR2dp2 + dR3dp2) / 2.0;
    }

    // drho1/dt
    const double dr1dt = R1 + dp1dt / (c1 * c1);

    const double dr1dtdr1 = dR1dr1 + dp1dtdr1 / (c1 * c1) - 2.0 * dp1dt * dc1dr1 / pow(c1, 3);
    const double dr1dtdu1 = dR1du1 + dp1dtdu1 / (c1 * c1);
    const double dr1dtdp1 = dR1dp1 + dp1dtdp1 / (c1 * c1) - 2.0 * dp1dt * dc1dp1 / pow(c1, 3);
    const double dr1dtdr2 = dR1dr2 + dp1dtdr2 / (c1 * c1);
    const double dr1dtdu2 = dR1du2 + dp1dtdu2 / (c1 * c1);
    const double dr1dtdp2 = dR1dp2 + dp1dtdp2 / (c1 * c1);

    // du1/dt
    const double du1dt = (R2 - dp1dt) / (r1 * c1);

    const double du1dtdr1 = ( (dp1dt - R2) * r1 * dc1dr1
                            + c1 * (dp1dt - R2 - r1 * dp1dtdr1 + r1 * dR2dr1) )
                            / (r1 * c1 * r1 * c1);
    const double du1dtdu1 = (dR2du1 - dp1dtdu1) / (r1 * c1);
    const double du1dtdp1 = ( (dp1dt - R2) * dc1dp1 + c1 * (dR2dp1 - dp1dtdp1) ) / (r1 * c1 * c1);
    const double du1dtdr2 = (dR2dr2 - dp1dtdr2) / (r1 * c1);
    const double du1dtdu2 = (dR2du2 - dp1dtdu2) / (r1 * c1);
    const double du1dtdp2 = (dR2dp2 - dp1dtdp2) / (r1 * c1);

    // d(ru)1/dt
//  double dru1dt;
//  dru1dt = r1 * du1dt + u1 * dr1dt + dr1dt * du1dt;
    const double dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1
                + dr1dtdr1 * du1dt + dr1dt * du1dtdr1;
    const double dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1
                + dr1dtdu1 * du1dt + dr1dt * du1dtdu1;
    const double dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1
                + dr1dtdp1 * du1dt + dr1dt * du1dtdp1;
    const double dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2
                + dr1dtdr2 * du1dt + dr1dt * du1dtdr2;
    const double dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2
                + dr1dtdu2 * du1dt + dr1dt * du1dtdu2;
    const double dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2
                + dr1dtdp2 * du1dt + dr1dt * du1dtdp2;

    // de1/dt
//  double de1dt;
//  de1dt = dp1dt * Cv / R + u1 * r1 * du1dt + u1 * u1 * dr1dt / 2.0;

    //const double de1dtdr1 = dp1dtdr1 * Cv / R + u1 * u1 * dr1dtdr1 / 2.0 + r1 * u1 * du1dtdr1
    //           + du1dt * u1;
    //const double de1dtdu1 = dp1dtdu1 * Cv / R + u1 * u1 * dr1dtdu1 / 2.0 + r1 * u1 * du1dtdu1
    //           + du1dt * r1 + dr1dt * u1;
    //const double de1dtdp1 = dp1dtdp1 / (gam - 1) + u1 * u1 * dr1dtdp1 / 2.0 + r1 * u1 * du1dtdp1;
    //const double de1dtdr2 = dp1dtdr2 / (gam - 1) + u1 * u1 * dr1dtdr2 / 2.0 + r1 * u1 * du1dtdr2;
    //const double de1dtdu2 = dp1dtdu2 / (gam - 1) + u1 * u1 * dr1dtdu2 / 2.0 + r1 * u1 * du1dtdu2;
    //const double de1dtdp2 = dp1dtdp2 / (gam - 1) + u1 * u1 * dr1dtdp2 / 2.0 + r1 * u1 * du1dtdp2;

    const double de1dtdr1 = (2.0 * Cv * dp1dtdr1
               + R * (du1dt * (du1dt + 2.0 * u1)
               + pow(du1dt + u1, 2.0) * dr1dtdr1
               + 2.0 * (dr1dt + r1) * (du1dt + u1) * du1dtdr1))
               / (2.0 * R);
    const double de1dtdu1 = (2.0 * Cv * dp1dtdu1
               + R * (pow(du1dt + u1, 2.0) * dr1dtdu1
               + 2.0 * (du1dt * r1 + dr1dt * (du1dt + u1)
               + (dr1dt + r1) * (du1dt + u1) * du1dtdu1)))
               / (2.0 * R);
    const double de1dtdp1 = (2.0 * Cv * dp1dtdp1
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp1
               + 2.0 * (dr1dt + r1) * du1dtdp1))
               / (2.0 * R);
    const double de1dtdr2 = (2.0 * Cv * dp1dtdr2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdr2
               + 2.0 * (dr1dt + r1) * du1dtdr2))
               / (2.0 * R);
    const double de1dtdu2 = (2.0 * Cv * dp1dtdu2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdu2
               + 2.0 * (dr1dt + r1) * du1dtdu2))
               / (2.0 * R);
    const double de1dtdp2 = (2.0 * Cv * dp1dtdp2
               + R * (du1dt + u1) * ((du1dt + u1) * dr1dtdp2
               + 2.0 * (dr1dt + r1) * du1dtdp2))
               / (2.0 * R);

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
	eval_dWpdW(gam, W[i1*3+0], W[i1*3+1], &dwpdw);

    for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			for (int k = 0; k < 3; k++) {
				dBodWo[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
			}
		}
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

    // Get Transformation Matrix
	eval_dWpdW(gam, W[i2*3+0], W[i2*3+1], &dwpdw);

    for (int row = 0; row < 3; row++)
    for (int col = 0; col < 3; col++)
    for (int k = 0; k < 3; k++)
        dBodWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];

}

