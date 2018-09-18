#include <vector>
#include <iostream>
#include <stdio.h>
#include "globals.h"
#include "convert.h"

double TotalPressureLoss(std::vector<double> W);

void ioTargetPressure(int io, std::vector<double> &p);

double inverseFitness(
    std::vector<double> pcurrent,
    std::vector<double> ptarget,
    std::vector<double> dx);

double evalFitness(
    std::vector<double> dx,
    std::vector<double> W)
{
    // Compute Fitness
    if (fitnessFun == 0)
        return TotalPressureLoss(W);
    else if (fitnessFun == 1)
    {
        std::vector<double> pcurrent(n_elem, 0);
        getp(W, pcurrent);
        std::vector<double> ptarget(n_elem, 0);
        ioTargetPressure(-1, ptarget);
        return inverseFitness(pcurrent, ptarget, dx);
    }
	abort();
	return -999.999;
}

double TotalPressureLoss(std::vector<double> W)
{
    double rhoout = W[(n_elem - 1) * 3 + 0];
    double uout = W[(n_elem - 1) * 3 + 1] / rhoout;
    double pout = (gam - 1) * (W[(n_elem - 1) * 3 + 2] - rhoout * pow(uout, 2) / 2);
    //double Tout = pout/(rhoout * R);

    double ptout_normalized;

    double ToverTt = 1 - pow(uout, 2) / a2 * (gam - 1) / (gam + 1);

    double poverpt = pow(ToverTt, (gam / (gam - 1)));

    ptout_normalized = 1 - (pout / poverpt) / inlet_total_p;

    return ptout_normalized;
}

// Input/Output Target Pressure Distribution
void ioTargetPressure(int io, std::vector<double> &p)
{
    FILE *TargetP;
    int err;
    // Output
    if (io > 0)
    {
        TargetP = fopen("targetP.dat", "w");
        fprintf(TargetP, "%d\n", n_elem);
        for (int i = 0; i < n_elem; i++)
            fprintf(TargetP, "%.15f\n", p[i] / inlet_total_p);
    }
    // Input
    else
    {
        int nxT;

        TargetP = fopen("targetP.dat", "r");
        rewind(TargetP);
        err = fscanf(TargetP, "%d", &nxT);
        if (nxT!=n_elem) std::cout<< "n_elem and nxT are different for targetP";
        for (int iT = 0; iT < nxT; iT++)
        {
            err = fscanf(TargetP, "%lf", &p[iT]);
        }
        if (err != 1) std::cout<< "Err";
    }

    fclose(TargetP);
}

// Return Inverse Design Fitness

double inverseFitness(
    std::vector<double> pcurrent,
    std::vector<double> ptarget,
    std::vector<double> dx)
{
    double fit = 0;
    for (int i = 0; i < n_elem; i++)
    {
        fit += pow(pcurrent[i] / inlet_total_p - ptarget[i], 2) * dx[i];
    }
//  std::cout<<"InverseFitness =  "<<fit / 2<<std::endl;
    return fit / 2;
}
