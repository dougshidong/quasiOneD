#include"finitedifferences.hpp"
#include <vector>
#include <Eigen/Core>
#include "structures.hpp"
#include "fitness.hpp"
#include "grid.hpp"
#include "quasiOneD.hpp"
#include "gradient.hpp"

using namespace Eigen;

MatrixXd hessian_central_gradient(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design,
    double pert)
{
	int n_design_variables = design.n_design_variables;
    MatrixXd Hessian(n_design_variables, n_design_variables);

	// Copy design, area, and flow to perturbed
    struct Design<double> pert_design = design;
	pert_design.design_variables = design.design_variables; 
    std::vector<double> pert_area = area;
	class Flow_data<double> pert_flow = flow_data;
	pert_flow.W          = flow_data.W;
	pert_flow.W_stage    = flow_data.W_stage;
	pert_flow.fluxes     = flow_data.fluxes;
	pert_flow.residual   = flow_data.residual;

    VectorXd gradp(n_design_variables);
    VectorXd gradn(n_design_variables);

    Hessian.setZero();
    int k = 0;
    for (int i = 0; i < n_design_variables; i++)
    {
        printf("Solving gradient %d out of %d \n", k++, n_design_variables+1);
        pert_design.design_variables = design.design_variables;
		double dh = pert_design.design_variables[i] * pert;
		if(dh == 0) dh = pert;
        pert_design.design_variables[i] += dh;
        pert_area = evalS(pert_design, x, dx);
        quasiOneD(x, pert_area, flo_opts, &pert_flow);
		const int gradient_type = 1;
        gradp = getGradient(gradient_type, opt_opts.cost_function, x, dx, pert_area, flo_opts, pert_flow, opt_opts, pert_design);

        pert_design.design_variables = design.design_variables;
        pert_design.design_variables[i] -= dh;
        pert_area = evalS(pert_design, x, dx);
        quasiOneD(x, pert_area, flo_opts, &pert_flow);
        gradn = getGradient(gradient_type, opt_opts.cost_function, x, dx, pert_area, flo_opts, pert_flow, opt_opts, pert_design);
        for (int j = 0; j < n_design_variables; j++)
        {
            Hessian(i, j) += (gradp(j) - gradn(j)) / (2*dh);
            Hessian(j, i) += (gradp(j) - gradn(j)) / (2*dh);
        }
    }
    Hessian = Hessian / 2.0;
    return Hessian;
}
MatrixXd hessian_central(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design,
    double pert)
{
	int n_design_variables = design.n_design_variables;
    MatrixXd Hessian(n_design_variables, n_design_variables);

	// Copy design and area to perturbed
    struct Design<double> pert_design = design;
	pert_design.design_variables = design.design_variables; 
    std::vector<double> pert_area = area;
	class Flow_data<double> pert_flow = flow_data;
	//pert_flow.W          = flow_data.W;
	//pert_flow.W_stage    = flow_data.W_stage;
	//pert_flow.fluxes     = flow_data.fluxes;
	//pert_flow.residual   = flow_data.residual;

	quasiOneD(x, area, flo_opts, &pert_flow);
	double I = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

    int k = 0;
    for (int i = 0; i < n_design_variables; i++)
    for (int j = i; j < n_design_variables; j++)
    {
        printf("Solving flow %d out of %d \n", k++, (n_design_variables+1)*(n_design_variables)/2);
        double dhi = design.design_variables[i] * pert;
        double dhj = design.design_variables[j] * pert;
		if(dhi == 0) dhi = pert;
		if(dhj == 0) dhj = pert;
        if (i == j) {
            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] += dhi; pert_design.design_variables[j] += dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I1 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] += dhi;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I2 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] -= dhi;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I3 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] -= dhi; pert_design.design_variables[j] -= dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I4 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);
            Hessian(i, j) = (-I1 + 16*I2 - 30*I + 16*I3 - I4) / (12 * dhi * dhj);
        } else {
            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] += dhi; pert_design.design_variables[j] += dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I1 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] += dhi; pert_design.design_variables[j] -= dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I2 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] -= dhi; pert_design.design_variables[j] += dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I3 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            pert_design.design_variables = design.design_variables; pert_design.design_variables[i] -= dhi; pert_design.design_variables[j] -= dhj;
			pert_area = evalS(pert_design, x, dx);
            quasiOneD(x, pert_area, flo_opts, &pert_flow);
            double I4 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

            Hessian(i, j) = (I1 - I2 - I3 + I4) / (4 * dhi * dhj);
            Hessian(j, i) = Hessian(i, j);
        }
    }
    //Hessian = Hessian + Hessian.transposeInPlace();
    //Hessian = Hessian / 2.0;
    return Hessian;
}

