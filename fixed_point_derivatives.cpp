#include<math.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<Eigen/Sparse>

#include"boundary_conditions.hpp"
#include"boundary_gradient.hpp"
#include"timestep.hpp"
#include"grid.hpp"

#include<adolc/adolc.h>
#include<adolc/adolc_sparse.h>

using namespace Eigen;

void eval_dGdW_dGdX_adolc(
    const std::vector<double> &x,
	const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
	const struct Design<double> &design,
	SparseMatrix<double> *const dGdW,
	SparseMatrix<double> *const dGdX)
{
	const int n_dvar = design.n_design_variables;
	const int n_elem = flo_opts.n_elem;

    const int tag = 100;
    trace_on(tag);

	int n_indep = 0; // Expect n_face + n_elem+2
	int n_dep   = 0; // Expect n_elem

	const int n_state = 3;
	const int n_indep_w = n_state*(n_elem);
	const int n_indep_dvar = n_dvar;
	const int n_indep_expected = n_indep_w + n_indep_dvar;
	const int n_dep_expected = n_state*n_elem;
	double *indep = myalloc1(n_indep_expected); // freed
	double *dep = myalloc1(n_dep_expected); // freed

    class Flow_data<adouble> aflow_data(n_elem);
	for (int iw_elem = 1; iw_elem < n_elem+1; iw_elem++) {
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = iw_elem*n_state + iw_state;
			indep[n_indep] = flow_data.W[ki];
			aflow_data.W[ki] <<= indep[n_indep];
			n_indep++;
		}
    }
	for (int iw_state = 0; iw_state < n_state; iw_state++) {
		const int iw_elem_first = 0;
		const int ki_first = iw_elem_first*n_state + iw_state;
		aflow_data.W[ki_first] = flow_data.W[ki_first];
		const int iw_elem_last = n_elem+1;
		const int ki_last = iw_elem_last*n_state + iw_state;
		aflow_data.W[ki_last] = flow_data.W[ki_last];
	}

	
	struct Design<adouble> adesign;
	adesign.h					= design.h;
	adesign.t1 					= design.t1;
    adesign.t2 					= design.t2;
	adesign.spline_degree  		= design.spline_degree;
	adesign.parametrization  	= design.parametrization;
	adesign.n_design_variables  = design.n_design_variables;
	adesign.design_variables.resize(n_dvar);

	for (int i = 0; i < n_dvar; i++) {
		indep[n_indep] = design.design_variables[i];
        adesign.design_variables[i] <<= indep[n_indep];
		n_indep++;
    }

    std::vector<double> dx = eval_dx(x);
    std::vector<adouble> aarea = evalS(adesign, x, dx);
	stepInTime(flo_opts, aarea, dx, &aflow_data);

	for (int iw_elem = 1; iw_elem < n_elem+1; iw_elem++) {
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = iw_elem*n_state + iw_state;
			aflow_data.W[ki] >>= dep[n_dep];
			n_dep++;
		}
    }

    trace_off(tag);

	assert(n_dep == n_dep_expected);
	assert(n_indep == n_indep_expected);

    (*dGdW).setZero();
    (*dGdX).setZero();

	bool sparse = true;
	if (sparse) {
		unsigned int *rind  = NULL;        /* row indices    */
		unsigned int *cind  = NULL;        /* column indices */
		double       *values = NULL;       /* values         */
		int nnz;
		int options[4];

		// CFL computation results in non-sparse due to the 'fmax'
		options[0] = 0;          /* sparsity pattern by index domains (default) */ 
		options[1] = 1;          /*                         safe mode (default) */ 
		options[2] = 0;          /*              not required if options[0] = 0 */ 
		options[3] = 0;          /*                column compression (default) */ 
		sparse_jac(tag, n_dep, n_indep, 0, indep, &nnz, &rind, &cind, &values, options);

		for (int i = 0; i < nnz; i++) {
			const unsigned int row = rind[i];
			const unsigned int col = cind[i];
			const double val = values[i];
			//std::cout<<row<<" "<<col<<" "<<val<<std::endl;
			if((int)col < n_indep_w) {
				(*dGdW).insert(row, col) = val;
			} else {
				(*dGdX).insert(row, col-n_indep_w) = val;
			}

		}
		free(rind); rind=NULL;
		free(cind); cind=NULL;
		free(values); values=NULL;

	} else {
		double **jac = myalloc2(n_dep, n_indep); // freed
		jacobian(tag, n_dep, n_indep, indep, jac);
		for (int ir_elem = 0; ir_elem < n_elem; ir_elem++) {
			for (int ir_state = 0; ir_state < n_state; ir_state++) {

				const int row = ir_elem*n_state + ir_state;

				for (int iw_elem = ir_elem-1; iw_elem < ir_elem+2; iw_elem++) {
					if(iw_elem < 0 || iw_elem >= n_elem) continue;
					for (int iw_state = 0; iw_state < n_state; iw_state++) {
						const int col = iw_elem*n_state + iw_state;
						std::cout<<row<<" "<<col<<" "<<jac[row][col]<<std::endl;
						(*dGdW).insert(row, col) = jac[row][col];
					}
				}
			}
		}
	}
	myfree1(indep);
	myfree1(dep);

	return;
}
