#include"analyticHessian.hpp"
#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include"structures.hpp"
#include"convert.hpp"
#include"parametrization.hpp"
#include"residuald1.hpp"
#include"residuald2.hpp"
#include"cost_derivative.hpp"
#include"boundary_gradient.hpp"
#include"quasiOneD.hpp"
#include"flux.hpp"
#include"grid.hpp"
#include"solve_linear.hpp"
#include"gradient.hpp"
#include"finitedifferences.hpp"

using namespace Eigen;

std::vector <MatrixXd> evalddWdDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);

std::vector <MatrixXd> evalddWdDesdDes_FD(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Design<double> &design);

MatrixXd directAdjointHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);
MatrixXd adjointAdjointHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);
MatrixXd adjointDirectHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);
MatrixXd directDirectHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);

MatrixXd getAnalyticHessian(
    const int hessian_type,
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_dvar = design.n_design_variables;
    MatrixXd Hessian(n_dvar, n_dvar);

    switch(hessian_type) {
        case 0: Hessian = directDirectHessian(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design); break;
        case 1: Hessian = adjointDirectHessian(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design); break;
        case 2: Hessian = adjointAdjointHessian(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design); break;
        case 3: Hessian = directAdjointHessian(cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design); break;
        case -1: Hessian = hessian_central_gradient(x, dx, area, flo_opts, flow_data, opt_opts, design, opt_opts.default_pert); break;
        case -2: Hessian = hessian_central(x, dx, area, flo_opts, flow_data, opt_opts, design, opt_opts.default_pert); break;
        default: abort();
    }

    return Hessian;
}

MatrixXd directAdjointHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes [ n_elem + 1, n_dvar ]
    MatrixXd dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate ddAreadDesdDes [ n_elem+1, n_dvar, n_dvar ]
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1); // by (nDes, nDes)
    ddAreadDesdDes = evalddAreadDesdDes(x, dx, design);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dCostdW [3*n_elem]
    VectorXd dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdW [ 3*n_elem, 3*n_elem ]
    SparseMatrix<double> ddCostdWdW = evaldCostdWdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdDes [ 3*n_elem, n_dvar ]
    MatrixXd ddCostdWdDes = evalddCostdWdArea(n_elem) * dAreadDes;
    // Evaluate pCostpArea [ n_elem+1 ]
    VectorXd pCostpArea = evaldCostdArea(n_elem);
    // Evaluate dCostdDes [ n_dvar ]
    VectorXd dCostdDes = pCostpArea.transpose() * dAreadDes;
    // Evaluate ddCostdDesdDes [ n_dvar, n_dvar ]
    MatrixXd ddCostdDesdDes = dAreadDes.transpose() * evalddCostdAreadArea(n_elem) * dAreadDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
	MatrixXd dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdArea
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdArea = evalddRdWdArea(area, flo_opts, flow_data);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, n_dvar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdArea[Ri] * dAreadDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);

    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
      SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
      slusolver1.compute(-dRdW.transpose());
      if (slusolver1.info() != 0)
          std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    psi = slusolver1.solve(dCostdW);
    //psi = solveGMRES(-dRdW.transpose(),dCostdW);
    std::cout<<"adjoint"<<std::endl;
    std::cout<<psi<<std::endl;

    // *************************************
    // Evaluate dWdDes (n_dvar Flow Eval)
    // *************************************
    MatrixXd dWdDes(3 * n_elem, n_dvar);
    if (opt_opts.exact_hessian == 1) {
        dWdDes = solve_linear(-dRdW, dRdDes, 0, opt_opts.dwdx_tol);
        std::cout<<"dwddes"<<std::endl;
        std::cout<<dWdDes<<std::endl;
    }
    else if (opt_opts.exact_hessian < 0)
    {
        // Iterative Solution of dWdDes
        dWdDes = solve_linear(-dRdW, dRdDes, 3, opt_opts.dwdx_tol);

        std::cout<<"dWdDes ||r||/||b|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()/dRdDes.norm()<<std::endl;

        // Direct Solution of dWdDes
        MatrixXd realdWdDes(3*n_elem,n_dvar);
        realdWdDes = solve_linear(-dRdW, dRdDes, 0, opt_opts.dwdx_tol);

        std::cout<<"Relative error of approximate dWdDes vs exact dWdDes:"<<std::endl;
        std::cout<<(realdWdDes - dWdDes).norm()/realdWdDes.norm()<<std::endl;
        for (int icol = 0; icol < n_dvar; icol++) {
            std::cout << icol << "\t" << (realdWdDes.col(icol) - dWdDes.col(icol)).norm()/realdWdDes.col(icol).norm()<<std::endl;
        }
    }


    // *************************************
    // Evaluate total derivative DDCostDDesDDes
    // *************************************
    MatrixXd DDCostDDesDDes(n_dvar, n_dvar);
    DDCostDDesDDes.setZero();
    DDCostDDesDDes = ddCostdDesdDes;
    DDCostDDesDDes += dWdDes.transpose() * ddCostdWdDes;
    DDCostDDesDDes += (dWdDes.transpose() * ddCostdWdDes).transpose();
    DDCostDDesDDes += dWdDes.transpose() * ddCostdWdW * dWdDes;
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDCostDDesDDes += pCostpArea(Si) * ddAreadDesdDes[Si];
        DDCostDDesDDes += psi.dot(dRdArea.col(Si)) * ddAreadDesdDes[Si];
    }
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
//        DDCostDDesDDes += psi(Ri) * ddRdDesdDes; //ddRdDesdDes is 0
        DDCostDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri]);
        DDCostDDesDDes += (psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri])).transpose();
        DDCostDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdW[Ri] * dWdDes);
    }

    return DDCostDDesDDes;
}
MatrixXd adjointAdjointHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes [ n_elem + 1, n_dvar ]
    MatrixXd dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate ddAreadDesdDes [ n_elem+1, n_dvar, n_dvar ]
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1); // by (nDes, nDes)
    ddAreadDesdDes = evalddAreadDesdDes(x, dx, design);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dCostdW [3*n_elem]
    VectorXd dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdW [ 3*n_elem, 3*n_elem ]
    SparseMatrix<double> ddCostdWdW = evaldCostdWdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdDes [ 3*n_elem, n_dvar ]
    MatrixXd ddCostdWdDes = evalddCostdWdArea(n_elem) * dAreadDes;
    // Evaluate pCostpArea [ n_elem+1 ]
    VectorXd pCostpArea = evaldCostdArea(n_elem);
    // Evaluate dCostdDes [ n_dvar ]
    VectorXd dCostdDes = pCostpArea.transpose() * dAreadDes;
    // Evaluate ddCostdDesdDes [ n_dvar, n_dvar ]
    MatrixXd ddCostdDesdDes = dAreadDes.transpose() * evalddCostdAreadArea(n_elem) * dAreadDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
	MatrixXd dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdArea
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdArea = evalddRdWdArea(area, flo_opts, flow_data);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, n_dvar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdArea[Ri] * dAreadDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;

    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    // *************************************
    // Solve for Adjoint 1 psi(1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
    psi = slusolver1.solve(dCostdW);

    // *************************************
    // Solve for Adjoint 2 gamma (n_dvar Flow Eval)
    // *************************************
    MatrixXd gamma(3 * n_elem, n_dvar);
    gamma = slusolver2.solve(dRdDes);
    
    // *************************************
    // Solve for Adjoint 3 beta (n_dvar Flow Eval)
    // *************************************
    MatrixXd beta(3 * n_elem, n_dvar);
    MatrixXd RHS(3 * n_elem, n_dvar);

    RHS = ddCostdWdDes + ddCostdWdW * gamma;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        RHS += psi(Ri) * ddRdWdDes[Ri];
        RHS += psi(Ri) * ddRdWdW[Ri] * gamma;
    }
    beta = slusolver1.solve(RHS);

    dRdDes = dRdArea * dAreadDes;
    // *************************************
    // Evaluate total derivative DDCostDDesDDes
    // *************************************
    MatrixXd DDCostDDesDDes(n_dvar, n_dvar);
    DDCostDDesDDes.setZero();
    DDCostDDesDDes = ddCostdDesdDes
                   + gamma.transpose() * ddCostdWdDes
                   + beta.transpose() * dRdDes;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
 //     DDCostDDesDDes += psi(Ri) * ddRdDesdDes[Ri];
        DDCostDDesDDes += psi(Ri) * (gamma.transpose() * ddRdWdDes[Ri]);
    }
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDCostDDesDDes += pCostpArea(Si) * ddAreadDesdDes[Si];
        DDCostDDesDDes += psi.dot(dRdArea.col(Si)) * ddAreadDesdDes[Si];
    }

    return DDCostDDesDDes;
}

MatrixXd adjointDirectHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes [ n_elem + 1, n_dvar ]
    MatrixXd dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate ddAreadDesdDes [ n_elem+1, n_dvar, n_dvar ]
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1); // by (nDes, nDes)
    ddAreadDesdDes = evalddAreadDesdDes(x, dx, design);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dCostdW [3*n_elem]
    VectorXd dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdW [ 3*n_elem, 3*n_elem ]
    SparseMatrix<double> ddCostdWdW = evaldCostdWdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdDes [ 3*n_elem, n_dvar ]
    MatrixXd ddCostdWdDes = evalddCostdWdArea(n_elem) * dAreadDes;
    // Evaluate pCostpArea [ n_elem+1 ]
    VectorXd pCostpArea = evaldCostdArea(n_elem);
    // Evaluate dCostdDes [ n_dvar ]
    VectorXd dCostdDes = pCostpArea.transpose() * dAreadDes;
    // Evaluate ddCostdDesdDes [ n_dvar, n_dvar ]
    MatrixXd ddCostdDesdDes = dAreadDes.transpose() * evalddCostdAreadArea(n_elem) * dAreadDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
	MatrixXd dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdArea
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdArea = evalddRdWdArea(area, flo_opts, flow_data);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, n_dvar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdArea[Ri] * dAreadDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW.transpose()
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
    psi = slusolver1.solve(dCostdW);

    // *************************************
    // Evaluate dWdDes (n_dvar Flow Eval)
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if (slusolver2.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    MatrixXd dWdDes(3 * n_elem, n_dvar);
    dWdDes = solve_linear(-dRdW, dRdDes, 0, opt_opts.dwdx_tol);

    // *************************************
    // Solve for dpsidDes (n_dvar Flow Eval)
    // *************************************
    MatrixXd RHS(3 * n_elem, n_dvar);
    MatrixXd dpsidDes(3 * n_elem, n_dvar);
    RHS = ddCostdWdDes + ddCostdWdW * dWdDes;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        RHS += psi(Ri) * ddRdWdDes[Ri];
        RHS += psi(Ri) * ddRdWdW[Ri] * dWdDes;
    }
    dpsidDes = slusolver1.solve(RHS);

    // *************************************
    // Evaluate total derivative DDCostDDesDDes
    // *************************************
    MatrixXd DDCostDDesDDes(n_dvar, n_dvar);
    DDCostDDesDDes.setZero();
    DDCostDDesDDes = ddCostdDesdDes
                   + dWdDes.transpose() * ddCostdWdDes
                   + dRdDes.transpose() * dpsidDes;

    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        DDCostDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri]).transpose();
    }
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDCostDDesDDes += pCostpArea(Si) * ddAreadDesdDes[Si];
        DDCostDDesDDes += psi.dot(dRdArea.col(Si)) * ddAreadDesdDes[Si];
    }

    return DDCostDDesDDes;
}
MatrixXd directDirectHessian(
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes [ n_elem + 1, n_dvar ]
    MatrixXd dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate ddAreadDesdDes [ n_elem+1, n_dvar, n_dvar ]
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1); // by (nDes, nDes)
    ddAreadDesdDes = evalddAreadDesdDes(x, dx, design);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dCostdW [3*n_elem]
    VectorXd dCostdW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdW [ 3*n_elem, 3*n_elem ]
    SparseMatrix<double> ddCostdWdW = evaldCostdWdW(opt_opts, flo_opts, flow_data.W, dx);
    // Evaluate ddCostdWdDes [ 3*n_elem, n_dvar ]
    MatrixXd ddCostdWdDes = evalddCostdWdArea(n_elem) * dAreadDes;
    // Evaluate pCostpArea [ n_elem+1 ]
    VectorXd pCostpArea = evaldCostdArea(n_elem);
    // Evaluate dCostdDes [ n_dvar ]
    VectorXd dCostdDes = pCostpArea.transpose() * dAreadDes;
    // Evaluate ddCostdDesdDes [ n_dvar, n_dvar ]
    MatrixXd ddCostdDesdDes = dAreadDes.transpose() * evalddCostdAreadArea(n_elem) * dAreadDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
	MatrixXd dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdArea
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdArea = evalddRdWdArea(area, flo_opts, flow_data);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, n_dvar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdArea[Ri] * dAreadDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);

    // *************************************
    // Evaluate dWdDes (n_dvar Flow Eval)
    // *************************************
    double errdwddes = -1.0;
    MatrixXd dWdDes(3 * n_elem, n_dvar);
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > factdrdw;
    factdrdw.compute(-dRdW);
    if (factdrdw.info() != 0) std::cout<<"Factorization failed. Error: "<<factdrdw.info()<<std::endl;
    if (opt_opts.exact_hessian == 1) {
        dWdDes = solve_linear(-dRdW, dRdDes, 0, opt_opts.dwdx_tol);
    }
    else if (opt_opts.exact_hessian == -1 || opt_opts.exact_hessian == -3) {
        // Iterative Solution of dWdDes
        dWdDes = solve_linear(-dRdW, dRdDes, 3, opt_opts.dwdx_tol);

        std::cout<<"dWdDes ||r|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()<<std::endl;
        std::cout<<"dWdDes ||r||/||b|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()/dRdDes.norm()<<std::endl;

        // Direct Solution of dWdDes
        MatrixXd realdWdDes(3*n_elem,n_dvar);
        realdWdDes = solve_linear(-dRdW, dRdDes, 0, opt_opts.dwdx_tol);

        std::cout<<"Relative error of approximate dWdDes vs exact dWdDes:"<<std::endl;
        std::cout<<(realdWdDes - dWdDes).norm()/realdWdDes.norm()<<std::endl;
        errdwddes = (realdWdDes - dWdDes).norm()/realdWdDes.norm();
    }
    // *************************************
    // Evaluate ddWdDesdDes (n_dvar * (n_dvar+1) / 2 Flow Eval)
    // *************************************
    std::vector <MatrixXd> ddWdDesdDes(3 * n_elem);
    MatrixXd dummy(n_dvar, n_dvar);
    dummy.setZero();
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        ddWdDesdDes[Wi] = dummy;
    }

    if (opt_opts.exact_hessian >= -1) {
        VectorXd ddWdDesidDesj(3 * n_elem);
        VectorXd RHS(3 * n_elem);
        for (int di = 0; di < n_dvar; di++) {
            for (int dj = di; dj < n_dvar; dj++) {
                // Evaluate RHS
                RHS.setZero();
                for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
                    RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdArea[Ri] * dAreadDes.col(di);
                    RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdArea[Ri] * dAreadDes.col(dj);
                    RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
                }

                for (int Si = 0; Si < n_elem + 1; Si++) {
                    RHS += dRdArea.col(Si) * ddAreadDesdDes[Si](di, dj);
                }

                // Solve
                ddWdDesidDesj = factdrdw.solve(RHS);

                for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
                    ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                    if (di != dj) {
                        ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                    }
                }
            }
        }
    }

    // *************************************
    // Evaluate total derivative DDCostDDesDDes
    // *************************************
    MatrixXd DDCostDDesDDes(n_dvar, n_dvar);
    MatrixXd A(n_dvar, n_dvar);
    A = dWdDes.transpose() * ddCostdWdDes
        + (dWdDes.transpose() * ddCostdWdDes).transpose()
        + dWdDes.transpose() * ddCostdWdW * dWdDes
        + dWdDes.transpose() * ddCostdWdW * dWdDes
        ;
//      + dWdDes.transpose() * ddCostdWdW * dWdDes;
    DDCostDDesDDes = ddCostdDesdDes
                   + dWdDes.transpose() * ddCostdWdDes
                   + (dWdDes.transpose() * ddCostdWdDes).transpose()
                   + dWdDes.transpose() * ddCostdWdW * dWdDes;
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDCostDDesDDes += pCostpArea(Si) * ddAreadDesdDes[Si];
    }
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        DDCostDDesDDes += dCostdW(Wi) * ddWdDesdDes[Wi];
    }

    if (opt_opts.exact_hessian == -1) {
        for (int i = 0; i<n_dvar; i++) {
            for (int j = 0; j<n_dvar; j++) {
                std::cout << i << "\t" << j 
                    << "\t" << std::setprecision(5) 
                    << errdwddes / (DDCostDDesDDes(i,j)/A(i,j))
                    << std::endl;
            }
        }
    }

    return DDCostDDesDDes;
}

std::vector <MatrixXd> evalddWdDesdDes_FD(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;

    struct Design<double> pert_design = design;
    class Flow_data<double> pert_flow = flow_data;
    VectorXd W0(3 * n_elem),
             W1(3 * n_elem),
             W2(3 * n_elem),
             W3(3 * n_elem),
             W4(3 * n_elem);

    std::vector <MatrixXd> dwddesddes(3 * n_elem);
    MatrixXd dummy(n_dvar, n_dvar);
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        dwddesddes[Wi] = dummy;
    }

    std::vector<double> pertArea(n_elem + 1);

    double dhi, dhj;
    double h = 1e-3;

    quasiOneD(x, area, flo_opts, &pert_flow);
    for (int Wi = 0; Wi<3*n_elem; Wi++) {
        W0[Wi] = pert_flow.W[Wi];
    }

    for (int i = 0; i < n_dvar; i++) {
        for (int j = i; j < n_dvar; j++) {
            dhi = pert_design.design_variables[i] * h;
            dhj = pert_design.design_variables[j] * h;
            if (i == j)
            {
                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] += dhi;
                pert_design.design_variables[j] += dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W1[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] += dhi;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W2[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] -= dhi;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W3[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] -= dhi;
                pert_design.design_variables[j] -= dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W4[Wi] = pert_flow.W[Wi];
                }

                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    dwddesddes[Wi](i, j) =
                    (-W1[Wi] + 16*W2[Wi] - 30*W0[Wi] + 16*W3[Wi] - W4[Wi]) / (12 * dhi * dhj);
                }
            }
            else
            {
                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] += dhi;
                pert_design.design_variables[j] += dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W1[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] += dhi;
                pert_design.design_variables[j] -= dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W2[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] -= dhi;
                pert_design.design_variables[j] += dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W3[Wi] = pert_flow.W[Wi];
                }

                pert_design.design_variables = design.design_variables;
                pert_design.design_variables[i] -= dhi;
                pert_design.design_variables[j] -= dhj;
                pertArea = evalS(pert_design, x, dx);
                quasiOneD(x, pertArea, flo_opts, &pert_flow);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W4[Wi] = pert_flow.W[Wi];
                }

                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    dwddesddes[Wi](i, j) = (W1[Wi] - W2[Wi] - W3[Wi] + W4[Wi]) / (4 * dhi * dhj);
                    dwddesddes[Wi](j, i) = dwddesddes[Wi](i, j);
                }

            } // if diag
        }// dj loop
    }// di loop
    return dwddesddes;
}

std::vector <MatrixXd> evalddWdDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design)
{
    const int n_elem = flo_opts.n_elem;
    const int n_dvar = design.n_design_variables;
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dAreadDes [ n_elem + 1, n_dvar ]
    MatrixXd dAreadDes = evaldAreadDes(x, dx, design);
    // Evaluate ddAreadDesdDes [ n_elem+1, n_dvar, n_dvar ]
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1); // by (nDes, nDes)
    ddAreadDesdDes = evalddAreadDesdDes(x, dx, design);
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Evaluate dRdArea
	MatrixXd dRdArea = evaldRdArea(flo_opts, flow_data);
    // Evaluate dRdDes
    MatrixXd dRdDes = dRdArea * dAreadDes;
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdArea
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdArea = evalddRdWdArea(area, flo_opts, flow_data);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, n_dvar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdArea[Ri] * dAreadDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);

    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver;
    slusolver.analyzePattern(-dRdW);
    slusolver.factorize(-dRdW);
    if (slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

    // Solve dWdDes
    MatrixXd dWdDes(3 * n_elem, n_dvar);
    dWdDes = slusolver.solve(dRdDes);

    // Evaluate ddWdDesdDes
    std::vector <MatrixXd> ddWdDesdDes(3 * n_elem);
    MatrixXd dummy(n_dvar, n_dvar);
    VectorXd ddWdDesidDesj(3 * n_elem);
    VectorXd RHS(3 * n_elem);
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        ddWdDesdDes[Wi] = dummy;
    }

    for (int di = 0; di < n_dvar; di++) {
        for (int dj = di; dj < n_dvar; dj++) {
            // Evaluate RHS
            RHS.setZero();
            for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
                RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdArea[Ri] * dAreadDes.col(di);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdArea[Ri] * dAreadDes.col(dj);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
            }

            for (int Si = 0; Si < n_elem + 1; Si++) {
                RHS += dRdArea.col(Si) * ddAreadDesdDes[Si](di, dj);
            }

            // Solve
            ddWdDesidDesj = slusolver.solve(RHS);

            for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
                ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                if (di != dj)
                {
                    ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                }
            }
        }
    }

//  std::vector <MatrixXd> ddWdDesdDes_FD = evalddWdDesdDes_FD(x, dx, area, design);
//  std::cout<<"Comparing AN and FD ddWdDesdDes"<<std::endl;
//  for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
//      if (std::abs(ddWdDesdDes[Wi].norm()>1e-12))
//      {
//          std::cout<<"Wi "<<Wi<<" diff: "<<
//              (ddWdDesdDes_FD[Wi] - ddWdDesdDes[Wi]).norm()/ddWdDesdDes[Wi].norm()<<std::endl;
//      }
//      else
//      {
//          std::cout<<"Wi "<<Wi<<
//              " zero DD norm. FD norm: "<<(ddWdDesdDes_FD[Wi]).norm()<<std::endl;
//      }
//  }

    return ddWdDesdDes;
}
