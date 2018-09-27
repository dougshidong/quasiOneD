#include<Eigen/Core>
#include<vector>
#include<iostream>
#include"structures.h"
#include"grid.h"
#include"spline.h"

using namespace Eigen;

MatrixXd evaldAreadDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design)
{
	const int n_elem = x.size();
	const int n_dvar = design.n_design_variables;
    MatrixXd dAreadDes(n_elem + 1, n_dvar);
    if (design.parametrization == 0) 
    {
    //    dAreadDes.setZero();
    //    for (int Si = 1; Si < n_elem; Si++)
    //    {
    //        dAreadDes(Si, Si - 1) = 1;
    //    }
    }
    if (design.parametrization == 1)
    {
    //    double d1 = design.design_variables[0];
    //    double d2 = design.design_variables[1];
    //    double d3 = design.design_variables[2];
    //    double xh;
    //    for (int i = 0; i < n_elem + 1; i++)
    //    {
    //        if (i == 0 || i == n_elem)
    //        {
    //            dAreadDes(i, 0) = 0;
    //            dAreadDes(i, 1) = 0;
    //            dAreadDes(i, 2) = 0;
    //        }
    //        else
    //        {
    //            xh = fabs(x[i] - dx[i] / 2.0);
    //            dAreadDes(i, 0) = - pow(sin(PI * pow(xh, d2)), d3);
    //            dAreadDes(i, 1) = - d1 * d3 * PI * pow(xh, d2)
    //                              * cos(PI * pow(xh, d2)) * log(xh)
    //                              * pow(sin(PI * pow(xh, d2)), d3 - 1);
    //            dAreadDes(i, 2) = - d1 * log(sin(PI * pow(xh, d2)))
    //                              * pow(sin(PI * pow(xh, d2)), d3);
    //        }
    //    }
    }
    if (design.parametrization == 2)
    {
		std::vector<double> control_points;
		control_points.push_back(1); // Clamped
		control_points.insert(control_points.end(), design.design_variables.begin(), design.design_variables.end());
		control_points.push_back(1); // Clamped
		const int n_control_pts = n_dvar+2;
		const int spline_degree = design.spline_degree;

		MatrixXd dAreadSpline = evalSplineDerivative(n_control_pts, spline_degree, control_points, x, dx);

		const int block_start_i = 0;
		const int block_start_j = 1;
		const int i_size = n_elem + 1;
		const int j_size = n_dvar;
		dAreadDes = dAreadSpline.block(block_start_i, block_start_j, i_size, j_size); // Do not return the endpoints
    }
    return dAreadDes;
}

MatrixXd evalddAreadDesdDes_FD(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design,
    const int Si)
{
    const int n_elem = x.size();
    double pertdi, pertdj;
    double pert = 1e-2;
    std::vector<double> area(n_elem + 1);
    std::vector<double> S1(n_elem+1), S2(n_elem+1), S3(n_elem+1), S4(n_elem+1);
    struct Design pert_design = design;
    area = evalS(design, x, dx);

    MatrixXd ddAreadDesdDes(3,3);
    for (int di = 0; di < 3; di++)
    {
        for (int dj = 0; dj < 3; dj++)
        {
            pertdi = design.design_variables[di] * pert;
            pertdj = design.design_variables[dj] * pert;

            for (int m = 0; m<3; m++)
            {
                pert_design.design_variables[m] = design.design_variables[m];
            }

            if (di != dj)
            {
                // ddFluxdWdW1 = ddFlux(i)/dW(i-1)dW(i-1)
                // R1
                pert_design.design_variables[di] = design.design_variables[di] + pertdi;
                pert_design.design_variables[dj] = design.design_variables[dj] + pertdj;

                S1 = evalS(pert_design, x, dx);

                // R2
                pert_design.design_variables[di] = design.design_variables[di] + pertdi;
                pert_design.design_variables[dj] = design.design_variables[dj] - pertdj;

                S2 = evalS(pert_design, x, dx);

                // R3
                pert_design.design_variables[di] = design.design_variables[di] - pertdi;
                pert_design.design_variables[dj] = design.design_variables[dj] + pertdj;

                S3 = evalS(pert_design, x, dx);

                // R4
                pert_design.design_variables[di] = design.design_variables[di] - pertdi;
                pert_design.design_variables[dj] = design.design_variables[dj] - pertdj;

                S4 = evalS(pert_design, x, dx);


                ddAreadDesdDes(di, dj) = 
                    (S1[Si] - S2[Si] - S3[Si] + S4[Si]) / (4.0 * pertdi * pertdj);

                pert_design.design_variables[di] = design.design_variables[di];
                pert_design.design_variables[dj] = design.design_variables[dj];
            }
            else
            {
                // R1
                pert_design.design_variables[di] = design.design_variables[di] + 2.0 * pertdi;

                S1 = evalS(pert_design, x, dx);

                // R2
                pert_design.design_variables[di] = design.design_variables[di] + pertdi;

                S2 = evalS(pert_design, x, dx);
                // R3
                pert_design.design_variables[di] = design.design_variables[di] - pertdi;

                S3 = evalS(pert_design, x, dx);
                // R4
                pert_design.design_variables[di] = design.design_variables[di] - 2.0 * pertdi;
                S4 = evalS(pert_design, x, dx);

                ddAreadDesdDes(di, dj) =
                                (-S1[Si] + 16*S2[Si] - 30*area[Si] + 16*S3[Si] - S4[Si])
                                             / (12 * pertdi * pertdj);
            }
        }
    }
    return ddAreadDesdDes;

}
std::vector <MatrixXd> evalddAreadDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design)
{
    const int PI = atan(1.0) * 4.0;
	const int n_elem = x.size();
    const int n_dvar = design.n_design_variables;
    std::vector <MatrixXd> ddAreadDesdDes(n_elem + 1);
    MatrixXd ddAreaidDesdDes(n_dvar, n_dvar);
    ddAreaidDesdDes.setZero();
    if (design.parametrization == 0 || design.parametrization == 2) {
        for (int Si = 0; Si < n_elem + 1; Si++) {
            ddAreadDesdDes[Si] = ddAreaidDesdDes;
        }
    } else if (design.parametrization == 1) {
        double d1 = design.design_variables[0];
        double d2 = design.design_variables[1];
        double d3 = design.design_variables[2];
        double xh;
        double xd2, spxd2, cpxd2 ;
        for (int Si = 0; Si < n_elem + 1; Si++)
        {
            if (Si == 0 || Si == n_elem)
            {
                ddAreadDesdDes[Si] = Matrix3d::Zero();
            }
            else
            {
                xh = fabs(x[Si] - dx[Si] / 2.0);
                xd2 = pow(xh, d2);
                spxd2 = sin(PI * xd2);
                cpxd2 = cos(PI * xd2);
                // Diagonal
                ddAreaidDesdDes(0, 0) = 0.0;
                ddAreaidDesdDes(1, 1) = -0.5 * d1 * d3 * PI * xd2 * pow(log(xh), 2.0)
                                    * pow(spxd2, d3-2.0)
                                    * ((d3-2.0) * PI * xd2 
                                       + d3 * PI * xd2  * cos(2.0 * PI * xd2)
                                       + sin(2.0 * PI * xd2));
                ddAreaidDesdDes(2, 2) = -d1 * pow(log(spxd2), 2.0) * pow(spxd2, d3);

                // Off-Diagonal
                ddAreaidDesdDes(0, 1) = -d3 * PI * xd2 * cpxd2 * log(xh) * pow(spxd2, d3-1.0);
                ddAreaidDesdDes(1, 0) = ddAreaidDesdDes(0, 1);

                ddAreaidDesdDes(0, 2) = -log(spxd2) * pow(spxd2, d3);
                ddAreaidDesdDes(2, 0) = ddAreaidDesdDes(0, 2);

                ddAreaidDesdDes(1, 2) = -d1 * PI * xd2 * cpxd2 * log(xh)
                                    * (1.0 + d3 * log(spxd2))
                                    * pow(spxd2, d3-1.0);
                ddAreaidDesdDes(2, 1) = ddAreaidDesdDes(1, 2);
                
                ddAreadDesdDes[Si] = ddAreaidDesdDes;
            }
        }
    }
    return ddAreadDesdDes;
}

