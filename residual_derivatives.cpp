#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include"boundary_conditions.hpp"
#include"timestep.hpp"
#include"residuald1.hpp"
#include<adolc/adolc.h>

using namespace Eigen;
SparseMatrix<double> eval_dRdW_dRdX_adolc(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const class Flow_data<double> &flow_data)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;

    const int tag = 100;
    trace_on(tag);

	int n_indep = 0; // Expect n_face + n_elem+2
	int n_dep   = 0; // Expect n_elem

	const int n_state = 3;
	const int n_indep_w = n_state*(n_elem+2);
	const int n_indep_area = n_face;
	const int n_indep_expected = n_indep_w + n_indep_area;
	const int n_dep_expected = n_state*n_elem;
	double *indep = myalloc1(n_indep_expected); // freed
	double *dep = myalloc1(n_dep_expected); // freed

    class Flow_data<adouble> aflow_data(n_elem);
	for (int iw_elem = 0; iw_elem < n_elem+2; iw_elem++) {
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = iw_elem*n_state + iw_state;
			indep[n_indep] = flow_data.W[ki];
			aflow_data.W[ki] <<= indep[n_indep];
			n_indep++;
		}
    }
    std::vector<adouble> aarea(n_face);
	for (int i = 0; i < n_face; i++) {
		indep[n_indep] = area[i];
        aarea[i] <<= indep[n_indep];
		n_indep++;
    }

    getDomainResi(flo_opts, aarea, aflow_data.W, &aflow_data.fluxes, &aflow_data.residual);

	for (int ir_elem = 1; ir_elem < n_elem+1; ir_elem++) {
		for (int ir_state = 0; ir_state < n_state; ir_state++) {
			const int ki = ir_elem*n_state + ir_state;
			aflow_data.residual[ki] >>= dep[n_dep];
			n_dep++;
		}
    }

    trace_off(tag);

	assert(n_dep == n_dep_expected);
	assert(n_indep == n_indep_expected);

	//unsigned int *rind  = NULL;        /* row indices    */
	//unsigned int *cind  = NULL;        /* column indices */
	//double       *values = NULL;       /* values         */
    //int nnz;
    //int options[4];

    //options[0] = 0;          /* sparsity pattern by index domains (default) */ 
    //options[1] = 0;          /*                         safe mode (default) */ 
    //options[2] = 0;          /*              not required if options[0] = 0 */ 
    //options[3] = 0;          /*                column compression (default) */ 
    //sparse_jac(tag, n_dep, n_indep, 0, indep, &nnz, &rind, &cind, &values, options);

	double **jac = myalloc2(n_dep, n_indep); // freed
    jacobian(tag, n_dep, n_indep, indep, jac);

	myfree1(indep);
	myfree1(dep);

	printf("\n");
	printf("\n");
    MatrixXd dRdW_dense(n_resi, n_resi);
	for (int ir_elem = 1; ir_elem < n_elem+1; ir_elem++) {
		for (int ir_state = 0; ir_state < n_state; ir_state++) {

			const int row = (ir_elem-1)*n_state + ir_state;

			for (int iw_elem = 1; iw_elem < n_elem+1; iw_elem++) {
				for (int iw_state = 0; iw_state < n_state; iw_state++) {
					const int col_withBC = iw_elem*n_state + iw_state;
					const int col = (iw_elem-1)*n_state + iw_state;
					dRdW_dense(row, col) = jac[row][col_withBC];
					printf("%d %d %f \n", row, col, jac[row][col_withBC]);
				}
			}
			printf("\n");
		}
	}
	printf("\n");
	printf("\n");
	std::cout<<dRdW_dense<<std::endl<<std::endl;

    SparseMatrix<double> dRdW_FD = evaldRdW_FD(area, flo_opts, flow_data);
	std::cout<<MatrixXd(dRdW_FD)<<std::endl;
	abort();

    SparseMatrix<double> dRdW = dRdW_dense.sparseView();

	myfree2(jac);

    return dRdW;
}

//void dWbcdW_adolc(
//    const Flow_options &flo_opts,
//    const class Flow_data<dreal> &flow_data,
//	MatrixXd *const dWidWi,
//	MatrixXd *const dWidWd,
//	MatrixXd *const dWodWo,
//	MatrixXd *const dWodWd)
//{
//    const int tag = 101;
//    trace_on(tag);
//
//	int n_indep = 0; // Expect n_face + n_elem+2
//	int n_dep   = 0; // Expect n_elem
//
//	const int n_state = 3;
//	const int n_indep_expected = n_state*2;
//	const int n_dep_expected = n_state;
//	double *indep = myalloc1(n_indep_expected); // freed
//	double *dep = myalloc1(n_dep_expected); // freed
//
//    class Flow_data<adouble> aflow_data(2);
//	for (int iw_elem = 0; iw_elem < 2; iw_elem++) {
//		for (int iw_state = 0; iw_state < n_state; iw_state++) {
//			const int ki = iw_elem*n_state + iw_state;
//			indep[n_indep] = flow_data.W[ki];
//			aflow_data.W[ki] <<= indep[n_indep];
//			n_indep++;
//		}
//    }
//	for (int iw_elem = 0; iw_elem < 2; iw_elem++) {
//		for (int iw_state = 0; iw_state < n_state; iw_state++) {
//			const int ki = (n_elem+iw_elem)*n_state + iw_state;
//			indep[n_indep] = flow_data.W[ki];
//			const int ki2 = (2+iw_elem)*n_state + iw_state;
//			aflow_data.W[ki2] <<= indep[n_indep];
//			n_indep++;
//		}
//    }
//
//    std::vector<adouble> aarea(n_face);
//	for (int i = 0; i < n_face; i++) {
//		indep[n_indep] = area[i];
//        aarea[i] <<= indep[n_indep];
//		n_indep++;
//    }
//
//	const int first_cell = 1;
//	const int last_cell = n_elem;
//	const adouble dt1 = flow_data.dt[first_cell];
//	const adouble dt2 = flow_data.dt[last_cell];
//	const adouble dx  = 1.0/n_elem;
//
//    getDomainResi(flo_opts, aarea, aflow_data.W, &aflow_data.fluxes, &aflow_data.residual);
//
//	for (int ir_elem = 1; ir_elem < n_elem+1; ir_elem++) {
//		for (int ir_state = 0; ir_state < n_state; ir_state++) {
//			const int ki = ir_elem*n_state + ir_state;
//			aflow_data.residual[ki] >>= dep[n_dep];
//			n_dep++;
//		}
//    }
//
//    trace_off(tag);
//
//	assert(n_dep == n_dep_expected);
//	assert(n_indep == n_indep_expected);
//
//	//unsigned int *rind  = NULL;        /* row indices    */
//	//unsigned int *cind  = NULL;        /* column indices */
//	//double       *values = NULL;       /* values         */
//    //int nnz;
//    //int options[4];
//
//    //options[0] = 0;          /* sparsity pattern by index domains (default) */ 
//    //options[1] = 0;          /*                         safe mode (default) */ 
//    //options[2] = 0;          /*              not required if options[0] = 0 */ 
//    //options[3] = 0;          /*                column compression (default) */ 
//    //sparse_jac(tag, n_dep, n_indep, 0, indep, &nnz, &rind, &cind, &values, options);
//
//	double **jac = myalloc2(n_dep, n_indep); // freed
//    jacobian(tag, n_dep, n_indep, indep, jac);
//
//	myfree1(indep);
//	myfree1(dep);
//
//	printf("\n");
//	printf("\n");
//    MatrixXd dRdW_dense(n_resi, n_resi);
//	for (int ir_elem = 1; ir_elem < n_elem+1; ir_elem++) {
//		for (int ir_state = 0; ir_state < n_state; ir_state++) {
//
//			const int row = (ir_elem-1)*n_state + ir_state;
//
//			for (int iw_elem = 1; iw_elem < n_elem+1; iw_elem++) {
//				for (int iw_state = 0; iw_state < n_state; iw_state++) {
//					const int col_withBC = iw_elem*n_state + iw_state;
//					const int col = (iw_elem-1)*n_state + iw_state;
//					dRdW_dense(row, col) = jac[row][col_withBC];
//					printf("%d %d %f \n", row, col, jac[row][col_withBC]);
//				}
//			}
//			printf("\n");
//		}
//	}
//	printf("\n");
//	printf("\n");
//	std::cout<<dRdW_dense<<std::endl<<std::endl;
//
//    SparseMatrix<double> dRdW_FD = evaldRdW_FD(area, flo_opts, flow_data);
//	std::cout<<MatrixXd(dRdW_FD)<<std::endl;
//	abort();
//
//    SparseMatrix<double> dRdW = dRdW_dense.sparseView();
//
//	myfree2(jac);
//
//    return dRdW;
//}
