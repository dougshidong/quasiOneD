#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <vector>
#include "flux.h"

const double PI = atan(1.0) * 4.0;

// Constants
double gam = 1.4;
double R = 1716;
double Cv = R / (gam - 1);

// Problem parameters
// Inlet
double Min = 0.12;
double Ttin = 531.2;
double ptin = 2117;
// Outlet
double pexit = 0.8 * ptin;
// Constant
double a2 = 2 * gam * Cv * Ttin * ((gam - 1) / (gam + 1)); // used in isentropic nozzle


// Convergence Settings
double CFL = 2;
double conv = 1e-13;
int maxIt = 60000;
int printIt = 1000; 
int printConv = 1; // 0 to hide real - time convergence
int printW = 0;
// Stepping Scheme
// 0   -   Euler Explicit
// 1   -   Runge - Kutta 4th order
int StepScheme = 1;
// Create Target Pressure
// 0   -   Do NOT Create Target Pressure
// 1   -   Create Target Pressure
int createTarget = 1;

double isenP(double pt, double M);

double isenT(double Tt, double M);

std::vector <double> calcVolume(std::vector <double> S, std::vector <double> dx);

void EulerExplicitStep(int nx, std::vector <double> S,
			std::vector <double> V,
			std::vector <double> dt,
			std::vector <double> Flux,
			std::vector <double> Q,
			std::vector <double> &Resi,
			std::vector <double> &W);

void rk4(int nx, std::vector <double> dx, std::vector <double> S, 
		std::vector <double> dt, 
	       	std::vector <double> &W,
	        std::vector <double> Q,	
		std::vector <double> &Resi,
		std::vector <double> &Flux);


double TotalPressureLoss(int nx, std::vector <double> W);

void ioTargetPressure(int io, int nx, std::vector <double> &p);

double inverseFitness(int nx, std::vector <double> pcurrent, std::vector <double> ptarget,
		std::vector <double> dx);

double quasiOneD(int nx, 
		std::vector <double> x, 
		std::vector <double> dx, 
		std::vector <double> S,
		int fitnessFun,
		std::vector <double> designVar,
		std::vector <double> &W)
{
	std::vector <double> rho(nx), u(nx), e(nx);
	std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
	std::vector <double> F(3 * nx, 0),Q(3 * nx, 0), Resi(3 * nx, 0);
	std::vector <double> Flux(3 * (nx + 1), 0);
//	std::vector <std::vector <double> > W(3,std::vector <double> (nx, 0)),
//					    F(3, std::vector <double> (nx, 0)),
//					    Q(3, std::vector <double> (nx, 0));

//	std::vector <std::vector <double> > Flux(3, std::vector <double> (nx + 1));

//	std::vector <std::vector <double> > Resi(3, std::vector <double> (nx, 0)),

	std::vector <double> dt(nx), V(nx);

	std::vector <int> itV(maxIt/printIt);
	std::vector <double> normV(maxIt/printIt);

	double dpdu, avgu, avgc, dtdx, dpdx, dudx;
	std::vector <double> eigenvalues(3);
	std::vector <double> charRel(3);
	double dp, drho, du;
	double MachBound;

	int iterlength;
	double normR = 1.0;
	int iterations = 0;


	V = calcVolume(S, dx);


	// Inlet flow properties
	Mach[0]=Min;
	T[0] = isenT(Ttin, Mach[0]);
	p[0] = isenP(ptin, Mach[0]);
	rho[0] = p[0] / (R * T[0]);
	c[0] = sqrt(gam * p[0] / rho[0]);
	u[0] = Mach[0] * c[0];
	e[0] = rho[0] * (Cv * T[0] + 0.5 * pow(u[0], 2));

	// Flow Properties Initialization
	for(int i = 1; i < nx; i++)
	{
		Mach[i] = Min;
		p[i] = pexit;
		T[i] = T[0];
		rho[i] = p[i] / (R * T[i]);
		c[i] = sqrt(gam * p[i] / rho[i]);
		u[i] = c[i] * Mach[i];
		e[i] = rho[i] * (Cv * T[i] + 0.5 * pow(u[i], 2));
	}

	// State and Flux Vectors Initialization
	for(int i = 0; i < nx; i++)
	{
		W[0 * nx + i] = rho[i];
		W[1 * nx + i] = rho[i] * u[i];
		W[2 * nx + i] = e[i];

		F[0 * nx + i] = rho[i] * u[i];
		F[1 * nx + i] = rho[i] * u[i] * u[i] + p[i];
		F[2 * nx + i] = ( e[i] + p[i] ) * u[i];
	}

	for(int i = 1; i < nx - 1; i++)
	{
		Q[1 * nx + i] = p[i] * (S[i] - S[i - 1]);
	}

	while(normR > conv && iterations < maxIt)
	{
		iterations++;

		if(iterations%printIt == 0) 
		{
			if(printConv == 1)
			{
				std::cout<<"Iteration "<<iterations
				         <<"   NormR "<<std::setprecision(15)<<normR<<std::endl;
			}
			itV[iterations / printIt - 1] = iterations;
			normV[iterations / printIt - 1] = normR;
		}

		for(int i = 0; i < nx; i++)
			dt[i] = (CFL * dx[i]) / fabs(u[i] + c[i]);

		Flux_StegerWarmingV(nx, Flux, W, u, c, rho);

		
		if(StepScheme == 0)
		{
			EulerExplicitStep(nx, S, V, dt, Flux, Q, Resi, W);
		}
		else if(StepScheme == 1)
		{
			rk4(nx, dx, S, dt, W, Q, Resi, Flux);
		}


		// Inlet Boundary Condition
		if(Mach[0] < 1)
		{
			dpdu = ptin * (gam / (gam - 1))
				 * pow(1 - ((gam - 1) / (gam + 1)) * u[0] * u[0] / a2,
                       1 / (gam - 1))
				 * ( - 2 * ((gam - 1) / (gam + 1)) * u[0] / a2);
                 
			eigenvalues[0] = ((u[1] + u[0] - c[1] - c[0]) / 2) * (dt[0] / dx[0]);
            
            dpdx = p[1] - p[0];
            dudx = u[1] - u[0];
			du = -eigenvalues[0] * (dpdx - rho[0] * c[0] * dudx)
				 / (dpdu - rho[0] * c[0]);

			u[0] = u[0] + du;
			T[0] = Ttin * (1 - ((gam - 1) / (gam + 1)) * u[0] * u[0] / a2);
			p[0] = ptin * pow(T[0] / Ttin, gam / (gam - 1));
			rho[0] = p[0] / (R * T[0]);
			e[0] = rho[0] * (Cv * T[0] + 0.5 * u[0] * u[0]);
			c[0] = sqrt(gam * p[0] / rho[0]);
			Mach[0] = u[0] / c[0];
		}

		// Exit boundary condition
        avgu = (u[nx - 1] + u[nx - 2]) / 2;
        avgc = (c[nx - 1] + u[nx - 2]) / 2;
        dtdx = dt[nx - 1] / dx[nx - 1];
		eigenvalues[0] = avgu * dtdx;
		eigenvalues[1] = (avgu + avgc) * dtdx;
		eigenvalues[2] = (avgu - avgc) * dtdx;

        dpdx = p[nx - 1] - p[nx - 2];
        dudx = u[nx - 1] - u[nx - 2];
		charRel[0] = -eigenvalues[0] * (rho[nx - 1] - rho[nx - 2]
                                        - (1 / (c[nx - 1] * c[nx - 1])) * dpdx);
		charRel[1] = -eigenvalues[1] * (dpdx + rho[nx - 1] * c[nx - 1] * dudx);
		charRel[2] = -eigenvalues[2] * (dpdx - rho[nx - 1] * c[nx - 1] * dudx);

		MachBound = avgu / avgc;
		//MachBound = (u[nx - 1] / (c[nx - 1]));
 		if(MachBound > 1)
        {
		    dp = 0.5 * (charRel[1] + charRel[2]);
        }
		else
        {
		    dp = 0;
        }
	
		drho = charRel[0] + dp / (pow(c[nx - 1], 2));
		du = (charRel[1] - dp) / (rho[nx - 1] * c[nx - 1]);
	
		u[nx - 1] = u[nx - 1] + du;
		rho[nx - 1] = rho[nx - 1] + drho;
		p[nx - 1] = p[nx - 1] + dp;
		T[nx - 1] = p[nx - 1] / (rho[nx - 1] * R);
		e[nx - 1] = rho[nx - 1] * (Cv * T[nx - 1] + 0.5 * pow(u[nx - 1], 2));
		c[nx - 1] = sqrt((gam * p[nx - 1]) / rho[nx - 1]);
		Mach[nx - 1] = u[nx - 1] / c[nx - 1];

		// Update flow properties
		for(int i = 1; i < nx - 1; i++)
		{
		    rho[i] = W[0 * nx + i];		// rho
		    u[i] = W[1 * nx + i] / rho[i];	// U
		    p[i] = (gam - 1) * (W[2 * nx + i] - rho[i] * pow(u[i], 2) / 2);  // Pressure
		    T[i] = p[i] / (rho[i] * R);	// Temperature
		    e[i] = W[2 * nx + i];		// Energy
		    c[i] = sqrt((gam * p[i]) / rho[i]);// Speed of sound
		    Mach[i] = u[i] / c[i];		// Mach number
		}

		// Update vectors at boundaries
		for(int i = 0; i < nx; i += nx - 1)
		{
		    W[0 * nx + i] = rho[i];
		    W[1 * nx + i] = rho[i] * u[i];
		    W[2 * nx + i] = e[i];
		    F[0 * nx + i] = rho[i] * u[i];
		    F[1 * nx + i] = rho[i] * u[i] * u[i] + p[i];
		    F[2 * nx + i] = (e[i] + p[i]) * u[i];
		    Q[0 * nx + i] = 0;
		    Q[1 * nx + i] = p[i] * (S[i+1] - S[i]);
		    Q[2 * nx + i] = 0;
		}
        
		// Calculating the norm of the density residual
		normR = 0;
		for(int i = 0; i < nx; i++)
		    normR = normR + Resi[0 * nx + i] * Resi[0 * nx + i];
		normR = sqrt(normR);
	}

	if(printW == 1)
	{
		for(int k = 0; k < 3; k++)
		{
			std::cout<<"W"<<k + 1<<std::endl;
			for(int i = 0; i < nx; i++)
				std::cout<<W[k * nx + i]<<std::endl;
		}
	}
	std::cout<<"Flow iterations = "<<iterations<<"   Density Residual = "<<normR<<std::endl;
	

	FILE  * Results;
	Results = fopen("Results.dat", "w");
	fprintf(Results, "%d\n", nx);
	for(int i = 0; i < nx; i++)
		fprintf(Results, "%.15f\n", x[i]);
	for(int i = 0; i < nx; i++)
		fprintf(Results, "%.15f\n", p[i] / ptin);
	for(int i = 0; i < nx; i++)
		fprintf(Results, "%.15f\n", rho[i]);
	for(int i = 0; i < nx; i++)
		fprintf(Results, "%.15f\n", Mach[i]);
	for(int i = 0; i < nx; i++)
		fprintf(Results, "%.15f\n", x[i] - dx[i] / 2);
	fprintf(Results, "%f\n", x.back() + dx.back() / 2);

	iterlength = itV.size();
	for(int i = 0; i < nx + 1; i++)
		fprintf(Results, "%.15f\n", S[i]);
	for(int i = 0; i < iterlength; i++)
		fprintf(Results, "%.15d\n", itV[i]);
	for(int i = 0; i < iterlength; i++)
		fprintf(Results, "%.15f\n", normV[i]);



	fclose(Results);

	// Create Target Pressure
	if(createTarget == 1) ioTargetPressure(1, nx, p);


	if(fitnessFun == 0)
		return TotalPressureLoss(nx, W);
	else if(fitnessFun == 1)
	{
		std::vector <double> ptarget(nx, 0);
		ioTargetPressure( - 1, nx, ptarget);
		return inverseFitness(nx, p, ptarget, dx);
	}


	return  - 9999.99;

}

double isenP(double pt, double M)
{
	return pt * pow((1 + (gam - 1) / 2 * pow(M, 2)), ( - gam / (gam - 1)));
}

double isenT(double Tt, double M)
{
	return Tt * pow((1 + (gam - 1) / 2 * pow(M, 2)), - 1);
}

double TotalPressureLoss(int nx, std::vector <double> W)
{
	double rhoout = W[0 * nx + (nx - 1)];
	double uout = W[1 * nx + (nx - 1)] / rhoout;
	double pout = (gam - 1) * (W[2 * nx + (nx - 1)] - rhoout * pow(uout, 2) / 2);
	//double Tout = pout/(rhoout * R);

	double ptout_normalized;

	double ToverTt = 1 - pow(uout, 2) / a2 * (gam - 1) / (gam + 1);

	double poverpt = pow(ToverTt, (gam / (gam - 1)));

	ptout_normalized = 1 - (pout / poverpt) / ptin;

	return ptout_normalized;
}

// Define Volume
std::vector <double> calcVolume(std::vector <double> S, std::vector <double> dx)
{
	std::vector <double> V;
	int ndx = dx.size();
	for(int i = 0; i < ndx; i++)
		V.push_back((S[i] + S[i + 1]) / 2 * dx[i]);

	return V;
}

// Input/Output Target Pressure Distribution
void ioTargetPressure(int io, int nx, std::vector <double> &p)
{

	FILE  * TargetP;
	// Output
	if(io > 0)
	{
		TargetP = fopen("targetP.dat", "w");
		fprintf(TargetP, "%d\n", nx);
//		for(int i = 0; i < nx; i++)
//			fprintf(TargetP, "%.15f\n", x[i]);
		for(int i = 0; i < nx; i++)
			fprintf(TargetP, "%.15f\n", p[i] / ptin);
	}
	// Input
	else
	{
		int nxT;

		TargetP = fopen("targetP.dat", "r");
		rewind(TargetP);
		fscanf(TargetP, "%d", &nxT);
        if(nxT!=nx) std::cout<< "nx and nxT are different for targetP";
		for(int iT = 0; iT < nxT; iT++)
		{
			fscanf(TargetP, "%lf", &p[iT]);
		}
	}	


	fclose(TargetP);

}

// Return Inverse Design Fitness

double inverseFitness(int nx, std::vector <double> pcurrent, std::vector <double> ptarget,
		std::vector <double> dx)
{
	double fit = 0;
	for(int i = 0; i < nx; i++)
	{
		fit += pow(pcurrent[i] / ptin - ptarget[i], 2) * dx[i];
	}
	std::cout<<"InverseFitness =  "<<fit / 2<<std::endl;
	return fit / 2;
}

// Jameson's 4th order Runge - Kutta Stepping Scheme
void rk4(int nx, std::vector <double> dx, std::vector <double> S, 
		std::vector <double> dt, 
	       	std::vector <double> &W,
	        std::vector <double> Q,	
		std::vector <double> &Resi,
		std::vector <double> &Flux)
{
	double ki;
	std::vector <double> Resi0(3 * nx, 0), Resi1(3 * nx, 0), Resi2(3 * nx, 0);
	std::vector <double> W1(3 * nx, 0), W2(3 * nx, 0), W3(3 * nx, 0);
	std::vector <double> Q1(3 * nx, 0), Q2(3 * nx, 0);
	std::vector <double> utemp(nx), rhotemp(nx), ptemp(nx), ctemp(nx);

	// Residual 0
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			Resi0[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q[ki + 1];
		}
	}
	// RK1
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			W1[ki] = W[ki] - (dt[i] / 2) * Resi0[ki] / dx[i];
		}
		W1[k * nx] = W[k * nx + 0];
		W1[k * nx + (nx - 1)] = W[k * nx + (nx - 1)];
	}
	for(int i = 0; i < nx; i++)
	{
		rhotemp[i] = W1[i];
		utemp[i] = W1[1 * nx + i] / rhotemp[i];
		ptemp[i] = (gam - 1) * (W1[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
		ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

		Q1[nx + i] = ptemp[i] * (S[i] - S[i - 1]);
	}

	Flux_StegerWarmingV(nx, Flux, W1, utemp, ctemp, rhotemp);

	// Residual 1
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			Resi1[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q1[ki];
		}
	}

	// RK2
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			W2[ki] = W[ki] - (dt[i] / 2) * Resi1[ki] / dx[i];
		}
		W2[k * nx] = W[k * nx + 0];
		W2[k * nx + (nx - 1)] = W[k * nx + (nx - 1)];
	}
	for(int i = 0; i < nx; i++)
	{
		rhotemp[i] = W1[i];
		utemp[i] = W1[1 * nx + i] / rhotemp[i];
		ptemp[i] = (gam - 1) * (W1[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
		ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

		Q2[nx + i] = ptemp[i] * (S[i] - S[i - 1]);
	}

	Flux_StegerWarmingV(nx, Flux, W2, utemp, ctemp, rhotemp);

	// Residual 2
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			Resi2[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q2[ki];
		}
	}

	// RK3
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			W3[ki] = W[ki] - (dt[i] / 2) * Resi2[ki] / dx[i];
		}
	}

	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			W[ki] = ((double)1.0 / 6.0) * (W[ki] + 2 * W1[ki] + 2 * W2[ki] + W3[ki]);
			Resi[ki] = (2 * Resi0[ki] + 2 * Resi1[ki] + Resi2[ki]) / 6.0;
		}
	}

	return;
}

// Euler Explicit
void EulerExplicitStep(int nx, std::vector <double> S,
			std::vector <double> V,
			std::vector <double> dt,
			std::vector <double> Flux,
			std::vector <double> Q,
			std::vector <double> &Resi,
			std::vector <double> &W)
{
	int ki;
	for(int k = 0; k < 3; k++)
	{
		for(int i = 1; i < nx - 1; i++)
		{
			ki = k * nx + i;
			Resi[ki] = Flux[ki + 1] * S[i + 1]
				 - Flux[ki] * S[i] - Q[ki];
		}
		Resi[k * nx + 0] = 0;
		Resi[k * nx + (nx - 1)] = 0;
	}
	for(int k = 0; k < 3; k++)
	for(int i = 0; i < nx - 1; i++)
	{
		ki = k * nx + i;
		W[ki] = W[ki] - (dt[i] / V[i - 1]) * Resi[ki];
	}

	return;
}
