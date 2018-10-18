#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include"boundary_conditions.h"
#include"timestep.h"
#include<adolc/adolc.h>

using namespace Eigen;
SparseMatrix<double> dRdW_adolc(
	const struct Flow_options<double> &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    const struct Flow_data<double> &flow_data)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;

    SparseMatrix<double> dRdW(n_resi, n_resi);
    dRdW.reserve(9*n_resi);

    const int tag = 100;
    trace_on(tag);

    struct Flow_data<adouble> aflow_data;
	for (int i = 0; i < n_elem+2; i++) {
        aflow_data.W[i] <<= flow_data.W[i];
    }

    getDomainResi(flo_opts, area, aflow_data.W, &flow_data.fluxes, &flow_data.residual);

	for (int i = 1; i < n_elem+1; i++) {
        aflow_data.residual[i] <<= flow_data.residual[i];
    }

    trace_off(tag);

    return dRdW;
}
