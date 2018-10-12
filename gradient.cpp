#include "structures.h"
#include "quasiOneD.h"
#include "fitness.h"
#include "grid.h"
#include <iomanip>
#include <iostream>
#include "flux.h"
#include "residuald1.h"
#include "parametrization.h"
#include "cost_derivative.h"

#include "solve_linear.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>

using namespace Eigen;

VectorXd gradient_adjoint(
	const int cost_function,
    const std::vector<double> &x,
	const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design);
VectorXd gradient_directDifferentiation(
	const int cost_function,
    const std::vector<double> &x,
	const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design);
VectorXd gradient_FD(
	const int FD_type,
	const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design,
	const double pert);

VectorXd getGradient(
	const int gradient_type,
	const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design)
{
    VectorXd grad(opt_opts.n_design_variables);
    if (gradient_type < 0) {;
        double pert = 1e-6;
		// -1 FFD
		// -2 BFD
		// -3 CFD
        grad = gradient_FD(gradient_type, cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design, pert);
    } else if (gradient_type == 1) {
        grad = gradient_adjoint(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);
    } else if (gradient_type == 2) {
        grad = gradient_directDifferentiation(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);
    }
    return grad;
}

VectorXd gradient_FD(
	const int FD_type,
	const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design,
	const double pert)
{
	int n_design_variables = design.n_design_variables;
    VectorXd grad(n_design_variables); // --> output

	// Copy design, area, and flow to perturbed
    struct Design pert_design = design;
	pert_design.design_variables = design.design_variables; 
    std::vector<double> pert_area = area;
	struct Flow_data pert_flow = flow_data;
	//pert_flow.W          = flow_data.W;
	//pert_flow.W_stage    = flow_data.W_stage;
	//pert_flow.fluxes     = flow_data.fluxes;
	//pert_flow.residual   = flow_data.residual;

    double I0;
    if (FD_type != -3) {
        quasiOneD(x, pert_area, flo_opts, &pert_flow);
		I0 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);
	}
    for (int i = 0; i < n_design_variables; i++) {
        double dh = design.design_variables[i] * pert;
		if(dh == 0) dh = pert;

		pert_design.design_variables = design.design_variables;
        if (FD_type == -1) // FFD
        {
			pert_design.design_variables[i] += dh;
            pert_area = evalS(pert_design, x, dx);
			quasiOneD(x, pert_area, flo_opts, &pert_flow);

			double I1 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);
            grad[i] = (I1 - I0) / dh;
        }
        else if (FD_type == -2) // BFD
        {
			pert_design.design_variables[i] -= dh;
            pert_area = evalS(pert_design, x, dx);
			quasiOneD(x, pert_area, flo_opts, &pert_flow);
			double I2 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);
            grad[i] = (I0 - I2) / dh;
        }
        else if (FD_type == -3) // CFD
        {
			pert_design.design_variables[i] += dh;
            pert_area = evalS(pert_design, x, dx);
			quasiOneD(x, pert_area, flo_opts, &pert_flow);
			double I1 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);

			pert_design.design_variables = design.design_variables;
			pert_design.design_variables[i] -= dh;
            pert_area = evalS(pert_design, x, dx);
			quasiOneD(x, pert_area, flo_opts, &pert_flow);
			double I2 = evalFitness(dx, flo_opts, pert_flow.W, opt_opts);
            grad[i] = (I1 - I2) / (2 * dh);
        }
    }

    return grad;
}

MatrixXd evaldWdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
    const struct Design &design,
	const int lsolver,
	const double tolerance)
{
    // DR   dR   dR DW           DW   -( dR ) ^ (-1) ( dR )
    // -- = -- + -- -- = 0  -->  -- =  ( -- )        ( -- )
    // DS   dS   dW DS           DS    ( dW )        ( dS )

	int n_elem = x.size();
    // Evaluate dRdArea
    MatrixXd dRdArea(3 * n_elem, n_elem + 1);
    dRdArea = evaldRdArea(flo_opts, flow_data);

    // Evaluate dAreadDes
    MatrixXd dAreadDes(n_elem + 1, design.n_design_variables);
    dAreadDes = evaldAreadDes(x, dx, design);

    //Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, design.n_design_variables);
    dRdDes = dRdArea * dAreadDes;

    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);

    // Solve DWDS
    MatrixXd dWdDes(3 * n_elem, design.n_design_variables);
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
	// 3 = GMRES
    dWdDes = solve_linear(-dRdW, dRdDes, lsolver, tolerance);

    return dWdDes;
}

VectorXd gradient_directDifferentiation(
	const int cost_function,
    const std::vector<double> &x,
	const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design)
{
	int n_elem = flo_opts.n_elem;
    // Direct Differentiation
    // I = Ic(W, area)
    // R = R(W, area) = 0 @ SS
    // W = W(W0, area)
    // DI   dCost   dCost DW
    // -- = --- + --- --
    // DS   dArea    dW  DS
    //
    // Evaluate dCostdW
    VectorXd dCostdW(3 * n_elem);
    dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate dAreadDes
    MatrixXd dAreadDes(n_elem + 1, design.n_design_variables);
    dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate dCostdArea
    VectorXd dCostdArea(n_elem + 1);
    dCostdArea = evaldCostdArea(n_elem);
    VectorXd dCostdDes(design.n_design_variables);
    dCostdDes = dCostdArea.transpose() * dAreadDes;

    // Evaluate dWdDes
    MatrixXd dWdDes(3 * n_elem, design.n_design_variables);
    dWdDes = evaldWdDes(x, dx, area, flo_opts, flow_data, design, 0, 1e-16); // Used for validate, so fast = better
    // Evaluate dCostdDes
    VectorXd dIdDes(design.n_design_variables);
    dIdDes = (dCostdDes.transpose() + dCostdW.transpose() * dWdDes);

    VectorXd grad(design.n_design_variables);
    grad = dIdDes;
    return grad;
}

VectorXd gradient_adjoint(
	const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
    const struct Design &design)
{
	int n_elem = flo_opts.n_elem;
	int n_design_variables = design.n_design_variables;
    VectorXd psi(n_elem*3);
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes
    MatrixXd dAreadDes(n_elem + 1, n_design_variables);
    dAreadDes = evaldAreadDes(x, dx, design);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dCostdW
    VectorXd dCostdW(3 * n_elem);
    dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate dCostdArea
    VectorXd dCostdArea(n_elem + 1);
    dCostdArea = evaldCostdArea(n_elem);
    // Evaluate dCostdDes
    VectorXd dCostdDes(n_design_variables);
    dCostdDes = dCostdArea.transpose() * dAreadDes;

    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
    MatrixXd dRdArea(3 * n_elem, n_elem + 1);
    dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, n_design_variables);
    dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);

    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    //VectorXd psi(3 * n_elem);
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(dRdW.transpose());
    if (slusolver1.info() != 0)
		std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    psi = slusolver1.solve(dCostdW);

    VectorXd dIdDes = dCostdDes.transpose() - psi.transpose()*dRdDes;

	return dIdDes;
}

