#include"fitness.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include "structures.h"
#include "convert.h"

double TotalPressureLoss(std::vector<double> W);

void ioTargetPressure(int io, std::vector<double> &p);

double inverse_pressure_design(
	const struct Flow_options<double> &flow_options,
    const std::vector<double> &W,
    const std::vector<double> &p_target,
    const std::vector<double> &dx);

double evalFitness(
    const std::vector<double> &dx,
	const struct Flow_options<double> &flow_options,
    const std::vector<double> &W,
	const struct Optimization_options<double> &opt_options)
{
    // Compute Fitness
    if (opt_options.cost_function == 0) {
        abort();//return TotalPressureLoss(W);
	} else if (opt_options.cost_function == 1) {
		double fitness = inverse_pressure_design(flow_options, W, opt_options.target_pressure, dx);
        return fitness;
    }
	abort();
}

//double TotalPressureLoss(std::vector<double> W)
//{
//    double rhoout = W[(n_elem - 1) * 3 + 0];
//    double uout = W[(n_elem - 1) * 3 + 1] / rhoout;
//    double pout = (gam - 1) * (W[(n_elem - 1) * 3 + 2] - rhoout * pow(uout, 2) / 2);
//    //double Tout = pout/(rhoout * R);
//
//    double ptout_normalized;
//
//    double ToverTt = 1 - pow(uout, 2) / a2 * (gam - 1) / (gam + 1);
//
//    double poverpt = pow(ToverTt, (gam / (gam - 1)));
//
//    ptout_normalized = 1 - (pout / poverpt) / inlet_total_p;
//
//    return ptout_normalized;
//}
//
//// Input/Output Target Pressure Distribution
//void ioTargetPressure(int io, std::vector<double> &p)
//{
//    FILE *TargetP;
//    int err;
//    // Output
//    if (io > 0)
//    {
//        TargetP = fopen("targetP.dat", "w");
//        fprintf(TargetP, "%d\n", n_elem);
//        for (int i = 0; i < n_elem; i++)
//            fprintf(TargetP, "%.15f\n", p[i] / inlet_total_p);
//    }
//    // Input
//    else
//    {
//        int nxT;
//
//        TargetP = fopen("targetP.dat", "r");
//        rewind(TargetP);
//        err = fscanf(TargetP, "%d", &nxT);
//        if (nxT!=n_elem) std::cout<< "n_elem and nxT are different for targetP";
//        for (int iT = 0; iT < nxT; iT++)
//        {
//            err = fscanf(TargetP, "%lf", &p[iT]);
//        }
//        if (err != 1) std::cout<< "Err";
//    }
//
//    fclose(TargetP);
//}

// Return Inverse Design Fitness

double inverse_pressure_design(
	const struct Flow_options<double> &flow_options,
    const std::vector<double> &W,
    const std::vector<double> &p_target,
    const std::vector<double> &dx)
{
    double fit = 0;
    const int n_elem = flow_options.n_elem;
    assert(p_target.size() == n_elem+2);
    for (int i = 1; i < n_elem+1; i++) {
		double p_current = get_p(flow_options.gam, W[i*3+0], W[i*3+1], W[i*3+2]);
        fit += pow(p_current / flow_options.inlet_total_p - p_target[i], 2) * dx[i];
    }
    return fit / 2;
}
