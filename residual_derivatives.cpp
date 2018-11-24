#include<math.h>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<adolc/adolc.h>
#include<adolc/adolc_sparse.h>
#include"solve_linear.hpp"
#include"boundary_conditions.hpp"
#include"boundary_gradient.hpp"
#include"timestep.hpp"
#include"residuald1.hpp"

using namespace Eigen;

void dWbcdW_adolc(
    const Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
	Matrix3d *const dWidWd,
	Matrix3d *const dWodWd);

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
	const int n_indep_w = n_state*(n_elem);
	const int n_indep_area = n_face;
	const int n_indep_expected = n_indep_w + n_indep_area;
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
    std::vector<adouble> aarea(n_face);
	for (int i = 0; i < n_face; i++) {
		indep[n_indep] = area[i];
        aarea[i] <<= indep[n_indep];
		n_indep++;
    }

	//inletBC(flo_opts, &aflow_data);
	//outletBC(flo_opts, &aflow_data);
    getDomainResi(flo_opts, aarea, &aflow_data);

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

    SparseMatrix<double> dRdW(n_resi, n_resi);
	const int n_stencil = 3;
	dRdW.reserve(n_elem*n_state*n_state*n_stencil);

	bool sparse = true;
	if (sparse) {
		unsigned int *rind  = NULL;        /* row indices    */
		unsigned int *cind  = NULL;        /* column indices */
		double       *values = NULL;       /* values         */
		int nnz;
		int options[4];

		options[0] = 0;          /* sparsity pattern by index domains (default) */ 
		options[1] = 0;          /*                         safe mode (default) */ 
		options[2] = 0;          /*              not required if options[0] = 0 */ 
		options[3] = 0;          /*                column compression (default) */ 
		sparse_jac(tag, n_dep, n_indep, 0, indep, &nnz, &rind, &cind, &values, options);

		for (int i = 0; i < nnz; i++) {
			const unsigned int row = rind[i];
			const unsigned int col = cind[i];
			const double val = values[i];
			if((int)col < n_indep_w) {
				dRdW.insert(row, col) = val;
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
						dRdW.insert(row, col) = jac[row][col];
					}
				}
			}
		}
	}
	myfree1(indep);
	myfree1(dep);

    //Matrix3d dWidWd, dWodWd, dRdW_block;
    //dWbcdW_adolc(flo_opts, flow_data, &dWidWd, &dWodWd);

    //for (int row = 0; row < n_state; row++) {
    //    for (int col = 0; col < n_state; col++) {
    //        dRdW_block(row, col) = jac[row][col];
    //    }
    //}
    //dRdW_block = dRdW_block * dWidWd;
    //for (int row = 0; row < n_state; row++) {
    //    for (int col = 0; col < n_state; col++) {
    //        dRdW.coeffRef(row, col) += dRdW_block(row, col);
    //    }
    //}

    //const int row_offset = n_state*(n_elem-1);
    //int col_offset = n_state*(n_elem+1);
    //for (int row = 0; row < n_state; row++) {
    //    for (int col = 0; col < n_state; col++) {
    //        dRdW_block(row, col) = jac[row+row_offset][col+col_offset];
    //    }
    //}

    //dRdW_block = dRdW_block * dWodWd;
    //col_offset = n_state*(n_elem-1);
    //for (int row = 0; row < n_state; row++) {
    //    for (int col = 0; col < n_state; col++) {
    //        dRdW.coeffRef(row+row_offset, col+col_offset) += dRdW_block(row, col);
    //    }
    //}


	//myfree2(jac);

    return dRdW;
}

void dWbcdW_adolc(
    const Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
	Matrix3d *const dWidWd,
	Matrix3d *const dWodWd)
{
	const int n_elem = flo_opts.n_elem;

    const int tag = 101;
    trace_on(tag);

	int n_indep = 0; // Expect n_face + n_elem+2
	int n_dep   = 0; // Expect n_elem

	const int n_state = 3;
	const int n_indep_expected = n_state*4;
	const int n_dep_expected = n_state*2;

	double *indep = myalloc1(n_indep_expected); // freed
	double *dep = myalloc1(n_dep_expected); // freed

    class Flow_data<adouble> aflow_data(2);
	for (int iw_elem = 0; iw_elem < 2; iw_elem++) {
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = iw_elem*n_state + iw_state;
			indep[n_indep] = flow_data.W[ki];
			aflow_data.W[ki] <<= indep[n_indep];
			n_indep++;
		}
    }
	for (int iw_elem = 0; iw_elem < 2; iw_elem++) {
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = (n_elem+iw_elem)*n_state + iw_state;
			indep[n_indep] = flow_data.W[ki];
			const int ki2 = (2+iw_elem)*n_state + iw_state;
			aflow_data.W[ki2] <<= indep[n_indep];
			n_indep++;
		}
    }

	const int first_cell = 1;
	const int last_cell = n_elem;
	const adouble dt1 = flow_data.dt[first_cell];
	const adouble dt2 = flow_data.dt[last_cell];
	const adouble dx  = 1.0/n_elem;

    for (int i = 0; i < 1; i++) {
        inletBC(flo_opts, &aflow_data);
        outletBC(flo_opts, &aflow_data);
    }

	{  // Inlet
        const int iw_elem = 0;
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			const int ki = iw_elem*n_state + iw_state;
			aflow_data.W[ki] >>= dep[n_dep];
			n_dep++;
		}
    }
	{  // Outlet
        const int iw_elem = 1;
		for (int iw_state = 0; iw_state < n_state; iw_state++) {
			//const int ki = (n_elem+iw_elem)*n_state + iw_state;
			const int ki2 = (2+iw_elem)*n_state + iw_state;
			aflow_data.W[ki2] >>= dep[n_dep];
			n_dep++;
		}
    }

    trace_off(tag);

	assert(n_dep == n_dep_expected);
	assert(n_indep == n_indep_expected);

	double **jac = myalloc2(n_dep, n_indep); // freed
    jacobian(tag, n_dep, n_indep, indep, jac);

	myfree1(indep);
	myfree1(dep);

    Matrix3d dBidWi, dBidWd, dBodWd, dBodWo;
    int row_offset = 0*n_state;
    int col_offset = 0*n_state;
    for (int row = 0; row < n_state; row++) {
        for (int col = 0; col < n_state; col++) {
            dBidWi(row, col) = jac[row+row_offset][col+col_offset];
        }
    }
    row_offset = 0*n_state;
    col_offset = 1*n_state;
    for (int row = 0; row < n_state; row++) {
        for (int col = 0; col < n_state; col++) {
            dBidWd(row, col) = jac[row+row_offset][col+col_offset];
        }
    }
    row_offset = 1*n_state;
    col_offset = 2*n_state;
    for (int row = 0; row < n_state; row++) {
        for (int col = 0; col < n_state; col++) {
            dBodWd(row, col) = jac[row+row_offset][col+col_offset];
        }
    }
    row_offset = 1*n_state;
    col_offset = 3*n_state;
    for (int row = 0; row < n_state; row++) {
        for (int col = 0; col < n_state; col++) {
            dBodWo(row, col) = jac[row+row_offset][col+col_offset];
        }
    }

    Matrix3d identity = MatrixXd::Identity(3,3);
    if ((identity-dBidWi).norm() == 0) {
        (*dWidWd).setZero();
    } else {
        (*dWidWd) = solve_dense_linear((identity - dBidWi), dBidWd, 0, 1e-16);
    }
    (*dWidWd) = dBidWd;


    //Vector3cd eigenvalues = (identity - dBodWo).eigenvalues();
    //double min_eig = eigenvalues.real().array().abs().minCoeff();
    //if (min_eig <= 1e-12) {
    //if ((identity-dBodWo).norm() <= 1e-12) {
    if ((identity-dBodWo).norm() == 0) {
        (*dWodWd).setZero();
    } else {
        (*dWodWd) = solve_dense_linear((identity - dBodWo), dBodWd, 0, 1e-16);
    }
	(*dWodWd) = dBodWd;
    //std::cout<<"(identity - dBodWo)"<<std::endl;
    //std::cout<<(identity - dBodWo)<<std::endl;
    //std::cout<<"(identity - dBodWo).eigenvalues())"<<std::endl;
    //std::cout<<(identity - dBodWo).eigenvalues()<<std::endl;
    //std::cout<<(dBodWd).eigenvalues()<<std::endl;

    //std::vector<double> dBidWi_AN(9), dBidWd_AN(9), dBodWd_AN(9), dBodWo_AN(9);
	//dRdW_BC_inlet(flo_opts, flow_data.W, dBidWi_AN, dBidWd_AN);
	//dRdW_BC_outlet(flo_opts, flow_data.W, dBodWd_AN, dBodWo_AN);

    //Matrix3d dBidWi_e, dBidWd_e, dBodWd_e, dBodWo_e, dWdW_e;
    //for (int row = 0; row < 3; row++) {
    //    for (int col = 0; col < 3; col++) {
    //        const int k = row * 3 + col;
    //        dBidWi_e(row,col) = dBidWi_AN.at(k);
    //        dBidWd_e(row,col) = dBidWd_AN.at(k);
    //        dBodWo_e(row,col) = dBodWo_AN.at(k);
    //        dBodWd_e(row,col) = dBodWd_AN.at(k);
    //    }
    //}
    //if ((identity-dBidWi_e).norm() == 0) {
    //    dWdW_e.setZero();
    //} else {
    //    dWdW_e = ((identity - dBidWi_e).inverse())*dBidWd_e;
    //}

    //std::cout<<"dWi"<<std::endl<<std::endl;
    //std::cout<<*dWidWd<<std::endl<<std::endl;
    //std::cout<<dWdW_e<<std::endl<<std::endl;
    //std::cout<<(*dWidWd-dWdW_e)<<std::endl<<std::endl;

    //std::cout<<dBidWd<<std::endl<<std::endl;
    //std::cout<<dBidWd_e<<std::endl<<std::endl;
    //std::cout<<(dBidWd-dBidWd_e)<<std::endl<<std::endl;

    //std::cout<<dBidWi<<std::endl<<std::endl;
    //std::cout<<dBidWi_e<<std::endl<<std::endl;
    //std::cout<<(dBidWi-dBidWi_e)<<std::endl<<std::endl;


    //if ((identity-dBodWo_e).norm() == 0) {
    //    dWdW_e.setZero();
    //} else {
    //    dWdW_e = ((identity - dBodWo_e).inverse())*dBodWd_e;
    //}
    //std::cout<<"dWodWd"<<std::endl<<std::endl;
    //std::cout<<*dWodWd<<std::endl<<std::endl;
    //std::cout<<dWdW_e<<std::endl<<std::endl;
    //std::cout<<(*dWodWd-dWdW_e)<<std::endl<<std::endl;

    //std::cout<<"dBodWd"<<std::endl<<std::endl;
    //std::cout<<dBodWd<<std::endl<<std::endl;
    //std::cout<<dBodWd_e<<std::endl<<std::endl;
    //std::cout<<(dBodWd-dBodWd_e)<<std::endl<<std::endl;

    //std::cout<<"dBodWo"<<std::endl<<std::endl;
    //std::cout<<dBodWo<<std::endl<<std::endl;
    //std::cout<<dBodWo_e<<std::endl<<std::endl;
    //std::cout<<(dBodWo-dBodWo_e)<<std::endl<<std::endl;
    

}
