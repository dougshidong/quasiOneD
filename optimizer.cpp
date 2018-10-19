#include "optimizer.hpp"
#include "structures.hpp"
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<stdlib.h>//exit
#include"quasiOneD.hpp"
#include"fitness.hpp"
#include"grid.hpp"
#include"gradient.hpp"
#include"analyticHessian.hpp"
#include"output.hpp"
#include<time.h>
#include<stdlib.h>     /* srand, rand */

using namespace Eigen;

void test_grad(
	const std::vector<double> &x,
	const std::vector<double> &dx, 
	const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design);
void test_hessian(
	const std::vector<double> &x,
	const std::vector<double> &dx, 
	const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design);

double linesearch_backtrack_unconstrained(
    const double initial_alpha,
    const std::vector<double> x,
    const std::vector<double> dx,
    const VectorXd &pk,
    const VectorXd &gradient,
    const double current_cost,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
    VectorXd* const searchD,
    class Flow_data<double>* const flow_data,
    struct Design<double>* const current_design);

MatrixXd BFGS(
    const MatrixXd &oldH,
    const VectorXd &oldg,
    const VectorXd &currentg,
    const VectorXd &searchD);

double checkCond(MatrixXd H);
MatrixXd invertHessian(MatrixXd H);
LLT<MatrixXd> checkPosDef(MatrixXd H);

VectorXd implicitSmoothing(VectorXd gradient, double epsilon);


void optimizer(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &initial_design)
{
	int n_elem = flo_opts.n_elem;
	int n_dvar = opt_opts.n_design_variables;
	class Flow_data<double> flow_data;

    const int n_face = n_elem + 1;
    const int n_elem_ghost = n_elem + 2;
    flow_data.dt.resize(n_elem_ghost);
    flow_data.W.resize(3*n_elem_ghost);
    flow_data.W_stage.resize(3*n_elem_ghost);
    flow_data.residual.resize(3*n_elem_ghost);
    flow_data.fluxes.resize(3*n_face);

	struct Design<double> current_design = initial_design;
	current_design.design_variables = initial_design.design_variables;

    std::vector<double> area = evalS<double>(current_design, x, dx);

    std::vector<double> gradient_norm_list, timeVec, Herror, svdvalues, svdvaluesreal, Hcond;
    MatrixXd H(n_dvar, n_dvar),
	         H_BFGS(n_dvar, n_dvar),
			 realH(n_dvar, n_dvar);

    VectorXd pk(n_dvar), searchD(n_dvar);

    clock_t tic = clock();

    std::ofstream myfile;
    myfile.open("convergence.dat");
    myfile << " Iteration \t Cost Function \t Gradient Norm \t Average Error \n";

	quasiOneD(x, area, flo_opts, &flow_data);
    double current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);

    VectorXd psi(3 * n_elem);

    VectorXd gradient(n_dvar);
    VectorXd oldGrad(n_dvar); //BFGS
    gradient = getGradient(opt_opts.gradient_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);

	bool testJacobian = true;
	//testJacobian = false;
	if (testJacobian) {
		test_grad(x, dx, area, flo_opts, flow_data, opt_opts, current_design);
		return;//exit(0);
	}
	bool testingGradient = true;
	//testingGradient = false;
	if (testingGradient) {
		test_grad(x, dx, area, flo_opts, flow_data, opt_opts, current_design);
        test_hessian(x, dx, area, flo_opts, flow_data, opt_opts, current_design);
		return;//exit(0);
	}

    // Initialize B
    H.setIdentity();
    H = H * 1.0;
    if (opt_opts.exact_hessian != 0) {
        H = getAnalyticHessian(opt_opts.hessian_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);
        checkCond(H);
        H = invertHessian(H);
    }

    double gradient_norm = 0;
    for (int i = 0; i < n_dvar; i++)
        gradient_norm += pow(gradient[i], 2);
    gradient_norm = sqrt(gradient_norm);
    int it_design = 0;

    // Design Loop
    while(gradient_norm > opt_opts.opt_tol && it_design < opt_opts.opt_maxit)
    {
        it_design++ ;

		printf("Iteration: %d, Cost Function %23.14e Gradient Norm: %23.14e \n", it_design, current_cost, gradient_norm);
		//std::cout<<"Current Design:\n";
		//for (int i = 0; i < n_dvar; i++) {
		//	std::cout<<designVar[i]<<std::endl;
		//}

//      1  =  Steepest Descent
//      2  =  Quasi-Newton (BFGS)
//      3  =  Newton
//      4  =  Truncated Newton with Adjoint-Direct Matrix-Vector Product
        if (opt_opts.descent_type == 1) {
            //int expo = rand() % 5 + 1 - 3;
            //double averageerr = 0;
            //for (int i = 0; i<n_dvar; i++) {
            //    srand (time(NULL));
            //    double fMin = -2.0;
            //    double fMax = 2.0;
            //    double expo = (double)rand() / RAND_MAX;
            //    expo =  fMin + expo * (fMax - fMin);
            //    expo =  0.0;
            //    pk(i) =  -10*gradient(i)*pow(10,expo);
            //    averageerr += fabs(expo)/n_dvar;
            //    //pk(i) =  -gradient(i);
            //}
            //myfile << it_design << "\t" << current_cost <<"\t"<< gradient_norm << "\t" << averageerr << "\n";
            //myfile.flush();
            pk =  -500*gradient;
        } else if (opt_opts.descent_type == 2) {
            if (it_design > 1) {
                H_BFGS = BFGS(H, oldGrad, gradient, searchD);
                H = H_BFGS;
				std::cout<<H.eigenvalues()<<std::endl;
            }
            pk = -H * gradient;
        } else if (opt_opts.descent_type == 3) {
            //H = getAnalyticHessian(x, dx, flow_data.W, area, designVar, hessian_type);
            //realH = getAnalyticHessian(x, dx, flow_data.W, area, designVar, 2);
            //double err = (realH - H).norm()/realH.norm();
            //std::cout<<"Hessian error: "<<err<<std::endl;
            //Herror.push_back(err);

            //checkCond(H);
            //H = invertHessian(H);
            //Hcond.push_back(checkCond(realH.inverse()));

            //pk = -H * gradient;
        }

		printf("%-15s %-15s\n", "Gradient","Search Direction");
		int step = 1;
		if(n_dvar/8>=1) step = n_dvar/8;
		if(step!=1) printf("Only printing 1 out of %d variables\n", step);
		for (int i=0; i<n_dvar; i+=step) {
			printf("%15.5e %15.5e\n", gradient[i], pk[i]);
		}

        double initial_alpha = 1.0;
		current_cost = linesearch_backtrack_unconstrained(
			initial_alpha, x, dx, pk, gradient, current_cost, flo_opts, opt_opts, &searchD, &flow_data, &current_design);

		area = evalS(current_design, x, dx);
        oldGrad = gradient;
		gradient = getGradient(opt_opts.gradient_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);

        gradient_norm = 0;
        for (int i = 0; i < n_dvar; i++)
            gradient_norm += pow(gradient[i], 2);
        gradient_norm = sqrt(gradient_norm);
        gradient_norm_list.push_back(gradient_norm);

		clock_t toc = clock();
        double elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
        timeVec.push_back(elapsed);
        std::cout<<"Time: "<<elapsed<<std::endl;

        std::cout<<"End of Design Iteration: "<<it_design<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    std::cout<<gradient<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for (int i = 0; i < n_dvar; i++)
        std::cout<<current_design.design_variables[i]<<std::endl;


    double final_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
    std::cout<<"Fitness: "<<final_cost<<std::endl;

    outVec(constants.case_name, "OptConv.dat", "w", gradient_norm_list);
    outVec(constants.case_name, "OptTime.dat", "w", timeVec);
    outVec(constants.case_name, "HessianErr.dat", "w", Herror);
    outVec(constants.case_name, "HessianCond.dat", "w", Hcond);
    outVec(constants.case_name, "svd.dat", "w", svdvalues);
    outVec(constants.case_name, "svdreal.dat", "w", svdvaluesreal);

    myfile.close();

    return;
}

double linesearch_backtrack_unconstrained(
    const double initial_alpha,
    const std::vector<double> x,
    const std::vector<double> dx,
    const VectorXd &pk,
    const VectorXd &gradient,
    const double current_cost,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
    VectorXd* const searchD,
    class Flow_data<double>* const flow_data,
    struct Design<double>* const current_design)
{
	int n_dvar = current_design->n_design_variables;
	double alpha = initial_alpha;
    double c1 = 1e-4;
    double c_pk_grad = c1 * gradient.dot(pk);

	// Copy current design
	struct Design<double> new_design = *current_design;
    for (int i = 0; i < n_dvar; i++) {
        new_design.design_variables[i] = current_design->design_variables[i] + alpha * pk[i];
    }

    std::vector<double> tempA = evalS(new_design, x, dx);
	quasiOneD(x, tempA, flo_opts, flow_data);
    double new_cost = evalFitness(dx, flo_opts, flow_data->W, opt_opts);
	while(new_cost > (current_cost + alpha * c_pk_grad))
    {
        alpha = alpha * 0.5;
        printf("Alpha Reduction: %e\n", alpha);
		if (alpha < 1e-14) {
			printf("Error. Can't find step size. Returning with tiny step.\n");
		}

        for (int i = 0; i < n_dvar; i++) {
            new_design.design_variables[i] = current_design->design_variables[i] + alpha * pk[i];
		}
        tempA = evalS(new_design, x, dx);
		quasiOneD(x, tempA, flo_opts, flow_data);
		new_cost = evalFitness(dx, flo_opts, flow_data->W, opt_opts);
        std::cout<<"new_cost: "<<new_cost<<std::endl;
        printf("current_cost + alpha/2.0 * c_pk_grad: %e\n", current_cost + alpha/ 2.0 * c_pk_grad);
    }

    *current_design = new_design;
    for (int i = 0; i < n_dvar; i++) {
        (*searchD)[i] = alpha * pk[i];
    }

    return new_cost;
}


double checkCond(MatrixXd H) {
    JacobiSVD<MatrixXd> svd(H);
    double svdmax = svd.singularValues()(0);
    double svdmin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = svdmax / svdmin;
    std::cout<<"Condition Number of H:"<<std::endl;
    std::cout<<cond<<std::endl;

    return cond;
}

MatrixXd invertHessian(MatrixXd H) {
    LLT<MatrixXd> llt = checkPosDef(H);
    return llt.solve(MatrixXd::Identity(H.rows(), H.rows()));
}

LLT<MatrixXd> checkPosDef(MatrixXd H) {
    LLT<MatrixXd> llt;
    VectorXcd eigval = H.eigenvalues();
    double shift = 1e-5;
    if (eigval.real().minCoeff() < 0) {
        MatrixXd eye(H.rows(),H.rows());
        eye.setIdentity();
        std::cout<<"Matrix is not Positive Semi-Definite"<<std::endl;
        std::cout<<"Eigenvalues:"<<std::endl;
        std::cout<<eigval<<std::endl;
        llt.compute(H + (shift - eigval.real().minCoeff()) * eye);
        checkCond(H + (shift - eigval.real().minCoeff()) * eye);
    }
    else {
        llt.compute(H);
    }
    return llt;
}

MatrixXd BFGS(
    const MatrixXd &oldH,
    const VectorXd &oldg,
    const VectorXd &currentg,
    const VectorXd &searchD)
{
	int n_rows = oldH.rows();
    MatrixXd newH(n_rows, n_rows);
    VectorXd dg(n_rows), dx(n_rows);
    MatrixXd dH(n_rows, n_rows), a(n_rows, n_rows), b(n_rows, n_rows);

    dg = currentg - oldg;
    dx = searchD;
	const double dgdx = dx.dot(dg);

	if(dgdx < 0) {
		printf("Negative curvature. Not updating BFGS \n");
		return oldH;
	}

    //a = ((dx.transpose() * dg + dg.transpose() * oldH * dg)(0) * (dx * dx.transpose()))
    //     / ((dx.transpose() * dg)(0) * (dx.transpose() * dg)(0));
    a = ((dgdx + dg.transpose() * oldH * dg) * (dx * dx.transpose())) / (dgdx*dgdx);
    b = (oldH * dg * dx.transpose() + dx * dg.transpose() * oldH) / (dgdx);

    dH = a - b;

    newH = oldH + dH;

    return newH;
}

VectorXd implicitSmoothing(VectorXd gradient, double epsilon)
{
    int n = gradient.size();
    MatrixXd A(n, n);
    A.setZero();

    for (int i = 0; i < n-1; i++) {
        A(i  , i) = 1.0 + 2.0 * epsilon;
        A(i+1, i) = -epsilon;
        A(i, i+1) = -epsilon;
    }
    A(n-1,n-1) = 1.0 + 2.0 * epsilon;

    LLT<MatrixXd> llt;
    llt.compute(A);
    if (llt.info() != 0)
        std::cout<<"Factorization failed. Error: "<<llt.info()<<std::endl;
    VectorXd smoothGrad = llt.solve(gradient);

    return smoothGrad;
}

void test_grad(
	const std::vector<double> &x,
	const std::vector<double> &dx, 
	const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design)
{
	printf("Testing gradient...\n");
	int n_dvar = design.n_design_variables;
    VectorXd adjoint_gradient(n_dvar),
             direct_gradient(n_dvar),
             cfinite_gradient(n_dvar),
             ffinite_gradient(n_dvar);

    printf("Evaluating adjoint gradient...\n");
    adjoint_gradient = getGradient(1, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);
    printf("Evaluating direct gradient...\n");
    direct_gradient = getGradient(2, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);
    printf("Evaluating central FD gradient...\n");
    cfinite_gradient = getGradient(-3, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);
    printf("Evaluating forward FD gradient...\n");
    ffinite_gradient = getGradient(-1, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, design);

	printf("%-23s %-23s %-23s %-23s %-23s %-23s %-23s %-23s\n", "Adjoint", "Direct-Diff", "Central FD", "ForwardFD", "AD-DM", "AD-FD","DM-FD", "CFD-FFD");
	//printf(" Adjoint                  Direct-Diff             Central FD              AD-DM                   AD-FD                   DM-FD\n");
	for (int i = 0; i < n_dvar; i++) {
		double g1 = adjoint_gradient[i];
		double g2 = direct_gradient[i];
		double g3 = cfinite_gradient[i];
		double g4 = ffinite_gradient[i];
		double r1 = (g1-g2)/(g1+1e-15);
		double r2 = (g1-g3)/(g1+1e-15);
		double r3 = (g2-g3)/(g2+1e-15);
		double r4 = (g3-g4)/(g3+1e-15);
		printf("%23.15e %23.15e %23.15e %23.15e %23.15e %23.15e %23.15e %23.15e\n", g1, g2, g3, g4, r1, r2, r3, r4);
	}
	return;
}
void test_hessian(
	const std::vector<double> &x,
	const std::vector<double> &dx, 
	const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design)
{
	printf("Testing Hessian...\n");
    struct Optimization_options<double> o_opts_copy = opt_opts; // Want to make sure it's the exact Hessian
    o_opts_copy.exact_hessian = 1;

    const int n_dvar = design.n_design_variables;
    MatrixXd DA(n_dvar, n_dvar),
             AA(n_dvar, n_dvar),
             AD(n_dvar, n_dvar),
             DD(n_dvar, n_dvar),
             FDG(n_dvar, n_dvar),
             FD2(n_dvar, n_dvar);
    printf("Evaluating DD...\n");
    DD = getAnalyticHessian(0, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
    printf("Evaluating AD...\n");
    AD = getAnalyticHessian(1, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
    printf("Evaluating AA...\n");
    AA = getAnalyticHessian(2, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
    printf("Evaluating DA...\n");
    DA = getAnalyticHessian(3, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
    printf("Evaluating FDG...\n");
    FDG = getAnalyticHessian(-1, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
    printf("Evaluating FD2...\n");
    FD2 = getAnalyticHessian(-2, o_opts_copy.cost_function, x, dx, area, flo_opts, flow_data, o_opts_copy, design);
	//printf(" Adjoint                  Direct-Diff             Central FD              AD-DM                   AD-FD                   DM-FD\n");

	printf("%-23s %-23s %-23s %-23s %-23s %-23s %-23s\n", "DD", "DD-AD", "DA-AA", "DD-DA", "DD-FDG", "DD-FD2", "FDG-FD2");
	//printf(" Adjoint                  Direct-Diff             Central FD              AD-DM                   AD-FD                   DM-FD\n");
	for (int i = 0; i < n_dvar; i++) {
        for (int j = i; j < n_dvar; j++) {
            double h1 = DD(i,j);
            double h2 = AD(i,j);
            double h3 = AA(i,j);
            double h4 = DA(i,j);
            double h5 = FDG(i,j);
            double h6 = FD2(i,j);
            double r2 = (h1-h2)/(h1+1e-15);
            double r3 = (h1-h3)/(h1+1e-15);
            double r4 = (h1-h4)/(h1+1e-15);
            double r5 = (h1-h5)/(h1+1e-15);
            double r6 = (h1-h6)/(h1+1e-15);
            double r7 = (h5-h6)/(h1+1e-15);
            printf("%23.15e %23.15e %23.15e %23.15e %23.15e %23.15e %23.15e\n", h1, r2, r3, r4, r5, r6, r7);
        }
	}
	return;
}
