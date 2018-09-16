#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <vector>
#include "flux.h"
#include "timestep.h"
#include "globals.h"
#include "output.h"
#include "fitness.h"

double isenP(double pt, double M);

double isenT(double Tt, double M);

void inletBC(std::vector<double> &W, std::vector<double> &Resi, double dt0, double dx0);
void outletBC(std::vector<double> &W, std::vector<double> &Resi, double dt0, double dx0);

double quasiOneD(std::vector<double> x, std::vector<double> S, std::vector<double> &W) {
    std::vector<double> dx(nx);
    for (int i = 1; i < nx - 1; i++) {
        dx[i] =  (x[i] - x[i - 1])/2  +  (x[i + 1] - x[i])/2 ;
    }
    dx[0] = x[1] - x[0];
    dx[nx-1] = x[x.size() - 1] - x[x.size() - 2];

    std::vector<double> rho(nx), u(nx), e(nx);
    std::vector<double> p(nx), c(nx);
    std::vector<double> Resi(3 * nx, 0);

    std::vector<double> dt(nx);

    std::vector <int> itV(maxIt/printIt);
    std::vector<double> normV(maxIt/printIt);
    std::vector<double> timeVec(maxIt/printIt);


    double normR = 1.0;
    int iterations = 0;


    // Inlet flow properties
    double inlet_T = isenT(inlet_total_T, inlet_mach);
    p[0] = isenP(inlet_total_p, inlet_mach);
    rho[0] = p[0] / (R * inlet_T);
    c[0] = sqrt(gam * p[0] / rho[0]);
    u[0] = inlet_mach * c[0];
    e[0] = rho[0] * (Cv * inlet_T + 0.5 * pow(u[0], 2));

    // Flow Properties Initialization
    for (int i = 1; i < nx; i++) {
        p[i] = outlet_p;
        rho[i] = p[i] / (R * inlet_T);
        c[i] = sqrt(gam * p[i] / rho[i]);
        u[i] = c[i] * inlet_mach;
        e[i] = rho[i] * (Cv * inlet_T + 0.5 * pow(u[i], 2));
    }

    // State Vectors Initialization
    for (int i = 0; i < nx; i++) {
        W[i * 3 + 0] = rho[i];
        W[i * 3 + 1] = rho[i] * u[i];
        W[i * 3 + 2] = e[i];
    }

    initializeTimeStep(nx); // Initiliaze Vector Sizes in timestep.cpp
    initializeFlux(nx); // Initiliaze Vector Sizes in flux.cpp

    clock_t tic = clock();
    clock_t toc;
    double elapsed;

    while(normR > flowConv && iterations < maxIt) {
        iterations++;

        // Calculate Time Step
        for (int i = 0; i < nx; i++) {
            dt[i] = (CFL * dx[i]) / fabs(u[i] + c[i]);
		}

        // Step in Time
        stepInTime(S, dx, dt, Resi, W);

        inletBC(W, Resi, dt[0], dx[0]);

        outletBC(W, Resi, dt[nx - 1], dx[nx - 1]);

        // Calculating the norm of the density residual
        normR = 0;
        for (int i = 0; i < nx; i++)
            normR = normR + pow(Resi[i * 3 + 0], 2);
        normR = sqrt(normR);

        // Monitor Convergence
        if (iterations%printIt == 0) {
            if (printConv == 1) {
                std::cout<<"Iteration "<<iterations <<"   NormR "<<std::setprecision(15)<<normR<<std::endl;
            }
            itV[iterations / printIt - 1] = iterations;
            normV[iterations / printIt - 1] = normR;

            toc = clock();
            elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
            timeVec[iterations / printIt - 1] = elapsed;
        }
    }

    if (printW == 1) {
        for (int k = 0; k < 3; k++) {
            std::cout<<"W"<<k + 1<<std::endl;
            for (int i = 0; i < nx; i++) {
                std::cout<<W[i * 3 + k]<<std::endl;
            }
        }
    }
    // Update flow properties
    for (int i = 0; i < nx ; i++) {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        e[i] = W[i * 3 + 2];
        p[i] = (gam - 1) * (e[i] - rho[i] * pow(u[i], 2) / 2);
        c[i] = sqrt((gam * p[i]) / rho[i]);
    }
    std::cout<<"Flow iterations = "<<iterations<<"   Density Residual = "<<normR<<std::endl;

    outVec("Geom.dat", "w", x);
    outVec("Geom.dat", "a", S);
    std::vector<double> pn(nx);
    for (int i = 0; i < nx; i++) pn[i] = p[i]/inlet_total_p;
    outVec("Flow.dat", "w", pn);
    outVec("Flow.dat", "a", rho);
	int nItPrint = iterations/printIt;
    std::vector<double> conv(nItPrint);
    for (int i = 0; i < nItPrint; i++) conv[i] = normV[i];
    outVec("FlowConv.dat", "w", conv);
    for (int i = 0; i < nItPrint; i++) conv[i] = timeVec[i];
    outVec("FlowTime.dat", "w", conv);

    return 1;
}

double isenP(double pt, double M) {
    return pt * pow((1 + (gam - 1) / 2 * pow(M, 2)), ( - gam / (gam - 1)));
}

double isenT(double Tt, double M) {
    return Tt * pow((1 + (gam - 1) / 2 * pow(M, 2)), - 1);
}


void inletBC(
    std::vector<double> &W,
    std::vector<double> &Resi,
    double dt0, double dx0)
{
    double dpdu, dpdx, dudx, dtdx, du, eigenvalue, T0;
    double rho[2], u[2], e[2], p[2], c[2];
    for (int i = 0; i < 2; i++) {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        e[i] = W[i * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
    if (u[0] < c[0]) {
        dpdu = inlet_total_p * (gam / (gam - 1.0))
             * pow(1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2,
                   1.0 / (gam - 1.0))
             * ( - 2.0 * ((gam - 1.0) / (gam + 1.0)) * u[0] / a2);
        dtdx = dt0 / dx0;
        eigenvalue = ((u[1] + u[0] - c[1] - c[0]) / 2.0) * dtdx;

        dpdx = p[1] - p[0];
        dudx = u[1] - u[0];
        du = -eigenvalue * (dpdx - rho[0] * c[0] * dudx)
                / (dpdu - rho[0] * c[0]);

//      Resi[0 * 3 + 1] = -((u[0] + du) - u[0]) / dtdx;
        u[0] = u[0] + du;

        T0 = inlet_total_T * (1.0 - ((gam - 1.0) / (gam + 1.0)) * u[0] * u[0] / a2);
//      Resi[0 * 3 + 2] = -(inlet_total_p * pow(T0 / inlet_total_T, gam / (gam - 1.0)) - p[0]) / dtdx;
        p[0] = inlet_total_p * pow(T0 / inlet_total_T, gam / (gam - 1.0));
        rho[0] = p[0] / (R * T0);
        e[0] = rho[0] * (Cv * T0 + 0.5 * u[0] * u[0]);

        Resi[0 * 3 + 0] = -(rho[0] - W[0 * 3 + 0]) / dtdx;
        Resi[0 * 3 + 1] = -(rho[0] * u[0] - W[0 * 3 + 1]) / dtdx;
        Resi[0 * 3 + 2] = -(e[0] - W[0 * 3 + 2]) / dtdx;

        W[0 * 3 + 0] = rho[0];
        W[0 * 3 + 1] = rho[0] * u[0];
        W[0 * 3 + 2] = e[0];
    } else {
        Resi[0 * 3 + 0] = 0;
        Resi[0 * 3 + 1] = 0;
        Resi[0 * 3 + 2] = 0;
    }
}

void outletBC(
    std::vector<double> &W,
    std::vector<double> &Resi,
    double dt0, double dx0)
{
    double avgc, avgu, dtdx, MachOut;
    double eigenvalues[3], Ri[3];
    double dpdx, dudx, du, drho, dp, T;
    double rho[2], u[2], e[2], p[2], c[2];

    for (int i = 0; i < 2; i++) {
        rho[i] = W[(i + (nx - 2)) * 3 + 0];
        u[i] = W[(i + (nx - 2)) * 3 + 1] / rho[i];
        e[i] = W[(i + (nx - 2)) * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }

    // Exit boundary condition
    avgu = (u[1] + u[0]) / 2;
    avgc = (c[1] + c[0]) / 2;
    dtdx = dt0 / dx0;
    eigenvalues[0] = avgu * dtdx;
    eigenvalues[1] = (avgu + avgc) * dtdx;
    eigenvalues[2] = (avgu - avgc) * dtdx;

    dpdx = p[1] - p[0];
    dudx = u[1] - u[0];

    Ri[0] = -eigenvalues[0] * ( rho[1] - rho[0] - dpdx / pow(c[1], 2) );
    Ri[1] = -eigenvalues[1] * ( dpdx + rho[1] * c[1] * dudx );
    Ri[2] = -eigenvalues[2] * ( dpdx - rho[1] * c[1] * dudx );

    MachOut = avgu / avgc;
    if (MachOut > 1) {
        dp = 0.5 * (Ri[1] + Ri[2]);
    } else {
        dp = 0;
    }

    drho = Ri[0] + dp / (pow(c[1], 2));
    du = (Ri[1] - dp) / (rho[1] * c[1]);

    u[1] = u[1] + du;
    rho[1] = rho[1] + drho;
    p[1] = p[1] + dp;
    T = p[1] / (rho[1] * R);
    e[1] = rho[1] * (Cv * T + 0.5 * pow(u[1], 2));

    Resi[(nx - 1) * 3 + 0] = (W[(nx - 1) * 3 + 0] - rho[1]) / dtdx;
    Resi[(nx - 1) * 3 + 1] = (W[(nx - 1) * 3 + 1] - rho[1] * u[1]) / dtdx;
    Resi[(nx - 1) * 3 + 2] = (W[(nx - 1) * 3 + 2] - e[1]) / dtdx;

    W[(nx - 1) * 3 + 0] = rho[1];
    W[(nx - 1) * 3 + 1] = rho[1] * u[1];
    W[(nx - 1) * 3 + 2] = e[1];
}
