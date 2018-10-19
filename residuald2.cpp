#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"convert.hpp"
#include"flux.hpp"
#include"quasiOneD.hpp"
#include"residuald1.hpp"
#include"boundary_conditions.hpp"
#include"boundary_gradient.hpp"
#include"hessianInlet.hpp"
#include"hessianOutlet.hpp"

using namespace Eigen;

std::vector <MatrixXd> evalddFdWdW(
    const double gam,
    const double rho,
    const double u,
    const double p);

void evalddFluxdWdW(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3);

//void evalddFluxdWdW_Scalar_FD(
//    std::vector <MatrixXd> &ddFluxdWdW1,
//    int ii, int jj,
//    std::vector<double> W);


std::vector <MatrixXd> evalddQdWdW(const double gam, const std::vector<double> &W);

std::vector <SparseMatrix<double> > evalddRdWdW_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data)
{
    const int n_elem = flo_opts.n_elem;
    const double dt = 1.0, dx = 1.0;
    const double gam = flo_opts.gam;
    std::vector <SparseMatrix<double> > ddRdWdW_FD(3 * n_elem);

    SparseMatrix<double> ddRidWdW_FD(3 * n_elem, 3 * n_elem);
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdW_FD[Ri] = ddRidWdW_FD;
    }
    class Flow_data<double> pert_flow = flow_data;

    std::vector<double> Q(3 * n_elem, 0);
    double h = 1e-3;
    for (int Ri = 1; Ri < n_elem - 1; Ri++) {
        for (int k = 0; k < 3; k++) {
            const int Rik = Ri * 3 + k;
            const int Rikp = (Ri + 1) * 3 + k;

            for (int Wi = (Ri - 1) * 3; Wi <= (Ri + 1) * 3 + 2; Wi++) {
                const double pertWi = flow_data.W[Wi] * h;
                for (int Wj = (Ri - 1) * 3; Wj <= (Ri + 1) * 3 + 2; Wj++) {
                    const double pertWj = flow_data.W[Wj] * h;
                    if (Wi != Wj) {
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        // R1
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                        WtoQ(gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi1 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R2
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                        WtoQ(gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi2 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R3
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                        WtoQ(gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi3 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R4
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                        WtoQ(gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi4 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                     / (4 * pertWi * pertWj);

                        // Symmetry
//                      ddRdWdW_FD[Rik].insert(Wj, Wi) = ddRdWdW_FD[Rik].coeffRef(Wi, Wj);
                        // Reset pert_flow.W
                        pert_flow.W[Wi] = flow_data.W[Wi];
                        pert_flow.W[Wj] = flow_data.W[Wj];
                    } else {
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }

                        WtoQ(flo_opts.gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi0 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R1
                        pert_flow.W[Wi] = flow_data.W[Wi] + 2.0 * pertWi;

                        WtoQ(flo_opts.gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi1 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R2
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;

                        WtoQ(flo_opts.gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi2 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R3
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;

                        WtoQ(flo_opts.gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi3 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        // R4
                        pert_flow.W[Wi] = flow_data.W[Wi] - 2.0 * pertWi;

                        WtoQ(flo_opts.gam, area, pert_flow.W, Q);
                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                        const double Resi4 = pert_flow.fluxes[Rikp] * area[Ri + 1] - pert_flow.fluxes[Rik] * area[Ri] - Q[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                                                     / (12 * pertWi * pertWj);
                        // Reset pert_flow.W
                        pert_flow.W[Wi] = flow_data.W[Wi];
                        pert_flow.W[Wj] = flow_data.W[Wj];
                    }
                }// Wj Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    //Inlet
    for (int Ri = 0; Ri < 1; Ri++) {
        for (int k = 0; k < 3; k++) {
            const int Rik = Ri * 3 + k;

            for (int Wi = Ri * 3; Wi <= (Ri + 1) * 3 + 2; Wi++) {
                const double pertWi = flow_data.W[Wi] * h;
                for (int Wj = Ri * 3; Wj <= (Ri + 1) * 3 + 2; Wj++) {
                    const double pertWj = flow_data.W[Wj] * h;
                    if (Wi != Wj) {
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        // R1
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi1 = pert_flow.residual[Rik];

                        // R2
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi2 = pert_flow.residual[Rik];

                        // R3
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi3 = pert_flow.residual[Rik];

                        // R4
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi4 = pert_flow.residual[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                     / (4 * pertWi * pertWj);
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                    } else {
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi0 = pert_flow.residual[Rik];

                        // R1
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] + 2.0 * pertWi;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi1 = pert_flow.residual[Rik];

                        // R2
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi2 = pert_flow.residual[Rik];

                        // R3
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi3 = pert_flow.residual[Rik];

                        // R4
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                        pert_flow.W[Wi] = flow_data.W[Wi] - 2.0 * pertWi;

                        inletBC(flo_opts, dt, dx, &pert_flow);
                        const double Resi4 = pert_flow.residual[Rik];

                        ddRdWdW_FD[Rik].insert(Wi, Wj) =
                            (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                            / (12 * pertWi * pertWj);
                        // Reset pert_flow.W
                        for (int m = 0; m < 3 * n_elem; m++) {
                            pert_flow.W[m] = flow_data.W[m];
                        }
                    }
                }// Wj Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    //Outlet
    int Ri = n_elem-1;
    for (int Rk = 0; Rk < 3; Rk++) {
        const int Rik = Ri * 3 + Rk;

        for (int Wi = (Ri - 1) * 3; Wi <= Ri * 3 + 2; Wi++) {
            const double pertWi = flow_data.W[Wi] * h;
            for (int Wj = (Ri - 1) * 3; Wj <= Ri * 3 + 2; Wj++) {
                const double pertWj = flow_data.W[Wj] * h;
                if (Wi != Wj) {
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    // R1
                    pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                    pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi1 = pert_flow.residual[Rik];

                    // R2
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
                    pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi2 = pert_flow.residual[Rik];

                    // R3
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                    pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi3 = pert_flow.residual[Rik];

                    // R4
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
                    pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi4 = pert_flow.residual[Rik];

                    ddRdWdW_FD[Rik].insert(Wi, Wj) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                 / (4 * pertWi * pertWj);
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                } else {
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi0 = pert_flow.residual[Rik];

                    // R1
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] + 2.0 * pertWi;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi1 = pert_flow.residual[Rik];

                    // R2
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi2 = pert_flow.residual[Rik];

                    // R3
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi3 = pert_flow.residual[Rik];

                    // R4
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    pert_flow.W[Wi] = flow_data.W[Wi] - 2.0 * pertWi;

                    outletBC(flo_opts, dt, dx, &pert_flow);
                    const double Resi4 = pert_flow.residual[Rik];

                    ddRdWdW_FD[Rik].insert(Wi, Wj) =
                        (-Resi1 + 16*Resi2 - 30*Resi0 + 16*Resi3 - Resi4)
                        / (12 * pertWi * pertWj);
                    // Reset pert_flow.W
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                } // If not diagonal
            }// Wj Loop
        }// Wi Loop
    }// Rk Loop
    return ddRdWdW_FD;
}

// Calculates Residual Hessian
std::vector < SparseMatrix<double> > evalddRdWdW(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data)
{
    const int n_elem = flo_opts.n_elem;
    std::vector <SparseMatrix<double> > ddRdWdW(3 * n_elem);
    SparseMatrix<double> dummySMat(3 * n_elem, 3 * n_elem);

    for (int Ri = 0; Ri < n_elem; Ri++) {
        for (int Rk = 0; Rk < 3; Rk++) {
            ddRdWdW[Ri * 3 + Rk] = dummySMat;
        }
    }
    // ********************************
    // BC Residual Hessian
    // ********************************
    // Inlet and Outlet Residual Hessian
    // ddRindWdW contains 3 matrices: dr1/dWdW, dr2/dWdW, dr3/dWdW
    // 6 state variables affect the inlet/outlet residuals
    std::vector <MatrixXd> ddRindWdW(3), ddRindWdW_FD(3);
    std::vector <MatrixXd> ddRoutdWdW(3), ddRoutdWdW_FD(3);
    MatrixXd dummyMat(6, 6);
    for (int Rk = 0; Rk < 3; Rk++) {
        ddRindWdW[Rk] = dummyMat;
        ddRoutdWdW[Rk] = dummyMat;
        ddRindWdW_FD[Rk] = dummyMat;
        ddRoutdWdW_FD[Rk] = dummyMat;
    }
    HessianInlet(flo_opts, flow_data.W, ddRindWdW);
    HessianOutlet(flo_opts, flow_data.W, ddRoutdWdW);
//  HessianBCprim_FD(W, ddRindWdW_FD, ddRoutdWdW_FD);
//  HessianBC_FD(W, ddRindWdW, ddRoutdWdW);


    for (int Rk = 0; Rk < 3; Rk++) {
      for (int row = 0; row < 6; row++) {
        for (int col = 0; col < 6; col++) {
          // Inlet
          ddRdWdW[Rk].insert(row, col) = -ddRindWdW[Rk](row, col);
          // Outlet
          int rowi = (n_elem - 2) * 3 + row;
          int coli = (n_elem - 2) * 3 + col;
          ddRdWdW[(n_elem - 1) * 3 + Rk].insert(rowi, coli) = -ddRoutdWdW[Rk](row, col);
        }
      }
    }

    // ********************************
    // Domain Residual Hessian
    // Validated with FD
    // ********************************
    // Evaluate ddFluxdWdW
    // ddFluxdWdW1 = Flux(i) / W(i-1, i-1)
    // ddFluxdWdW2 = Flux(i) / W(i  , i  )
    // ddFluxdWdW3 = Flux(i) / W(i-1, i  )
    // ddFluxdWdW4 = TRANSPOSE(ddFluxdWdW3)
    std::vector <MatrixXd> ddFluxdWdW1(3 * (n_elem + 1)),
                           ddFluxdWdW2(3 * (n_elem + 1)),
                           ddFluxdWdW3(3 * (n_elem + 1));
    evalddFluxdWdW(flo_opts, flow_data.W, ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3);

    // Evaluate ddQdWdW
    std::vector <MatrixXd> ddQdWdW(3 * n_elem);
    ddQdWdW = evalddQdWdW(flo_opts.gam, flow_data.W);

    int Wi1, Wi2;
    for (int Ri = 1; Ri < n_elem - 1; Ri++) {
        for (int Rk = 0; Rk < 3; Rk++) {
            const int Rik = Ri * 3 + Rk;
            const int Rikp = (Ri + 1) * 3 + Rk;

            // R(i) = Flux(i + 1) * area(i + 1) - Flux(i) * area(i) - Q(i)

            // ***************
            // Diagonal Terms
            // ***************

            //ddR(i) / dW(i-1)dW(i-1)
            // = - ddFlux(i) / dW(i-1)dW(i-1) * area(i)
            Wi1 = Ri - 1;
            Wi2 = Ri - 1;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                                -ddFluxdWdW1[Rik].coeffRef(row, col) * area[Ri];
                }
            }

            //ddR(i) / dW(i)dW(i)
            // = ddFlux(i + 1) / dW(i)dW(i) * area(i + 1)
            //   - ddFlux(i) / dW(i)dW(i) * area(i)
            //   - ddQ(i) / dW(i)dW(i) * dArea
            Wi1 = Ri;
            Wi2 = Ri;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                                ddFluxdWdW1[Rikp].coeffRef(row, col) * area[Ri + 1]
                                - ddFluxdWdW2[Rik].coeffRef(row, col) * area[Ri]
                                - ddQdWdW[Rik].coeffRef(row, col) * (area[Ri + 1] - area[Ri]);
                }
            }

            //ddR(i) / dW(i+1)dW(i+1)
            // = ddFlux(i+1) / dW(i+1)dW(i+1) * area(i+1)
            Wi1 = Ri + 1;
            Wi2 = Ri + 1;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                                ddFluxdWdW2[Rikp].coeffRef(row, col) * area[Ri + 1];
                }
            }

            // *******************
            // Off-Diagonal Terms
            // *******************

            // R(i) = Flux(i + 1) * area(i + 1) - Flux(i) * area(i) - Q(i)

            //ddR(i) / dW(i-1)dW(i)
            // = ddFlux(i) / dW(i-1)dW(i) * area(i)
            Wi1 = Ri - 1;
            Wi2 = Ri;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                        -ddFluxdWdW3[Rik].coeffRef(row, col) * area[Ri];
                }
            }
            //ddR(i) / dW(i)dW(i-1)
            // = (ddFlux(i) / dW(i-1)dW(i) * area(i)).transpose()
            Wi1 = Ri;
            Wi2 = Ri - 1;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                        -ddFluxdWdW3[Rik].coeffRef(col, row) * area[Ri];
                }
            }

            //ddR(i) / dW(i)dW(i+1)
            // = ddFlux(i) / dW(i)dW(i+1) * area(i+1)
            Wi1 = Ri;
            Wi2 = Ri + 1;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                        ddFluxdWdW3[Rikp].coeffRef(row, col) * area[Ri + 1];
                }
            }
            //ddR(i) / dW(i+1)dW(i) Symmetry
            // = (ddFlux(i) / dW(i)dW(i+1) * area(i+1)).transpose()
            Wi1 = Ri + 1;
            Wi2 = Ri;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddRdWdW[Rik].insert(Wi1 * 3 + row, Wi2 * 3 + col) =
                        ddFluxdWdW3[Rikp].coeffRef(col, row) * area[Ri + 1];
                }
            }

            //ddR(i) / dW(i-1)dW(i+1) = 0
            //ddR(i) / dW(i+1)dW(i-1) = 0
        }

    }
    return ddRdWdW;
}

std::vector <SparseMatrix<double> > evalddRdWdArea_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data)
{
    const int n_elem = flo_opts.n_elem;
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);

    SparseMatrix<double> ddRidWdArea(3 * n_elem, n_elem + 1);
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdArea[Ri] = ddRidWdArea;
    }

    class Flow_data<double> pert_flow = flow_data;
    std::vector<double> pertArea(n_elem + 1, 0);
    std::vector<double> Q(3 * n_elem, 0);
    double h = 1e-4;
    for (int Ri = 1; Ri < n_elem - 1; Ri++) {
        for (int k = 0; k < 3; k++) {
            const int Rik = Ri * 3 + k;
            const int Rikp = (Ri + 1) * 3 + k;

            for (int Wi = (Ri - 1) * 3; Wi <= (Ri + 1) * 3 + 2; Wi++) {
                const double dh_W = flow_data.W[Wi] * h;
                for (int Si = Ri; Si <= Ri+1; Si++) {
                    const double dh_area = area[Si] * h;
                    for (int m = 0; m < 3 * n_elem; m++) {
                        pert_flow.W[m] = flow_data.W[m];
                    }
                    for (int m = 0; m < n_elem + 1; m++) {
                        pertArea[m] = area[m];
                    }
                    // R1
                    pert_flow.W[Wi] = flow_data.W[Wi];
                    pertArea[Si] = area[Si];

                    pert_flow.W[Wi] = flow_data.W[Wi] + dh_W;
                    pertArea[Si] = area[Si] + dh_area;

                    WtoQ(flo_opts.gam, pertArea, pert_flow.W, Q);
                    getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                    const double Resi1 = pert_flow.fluxes[Rikp] * pertArea[Ri + 1] - pert_flow.fluxes[Rik] * pertArea[Ri] - Q[Rik];

                    // R2
                    pert_flow.W[Wi] = flow_data.W[Wi];
                    pertArea[Si] = area[Si];

                    pert_flow.W[Wi] = flow_data.W[Wi] + dh_W;
                    pertArea[Si] = area[Si] - dh_area;

                    WtoQ(flo_opts.gam, pertArea, pert_flow.W, Q);
                    getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                    const double Resi2 = pert_flow.fluxes[Rikp] * pertArea[Ri + 1] - pert_flow.fluxes[Rik] * pertArea[Ri] - Q[Rik];

                    // R3
                    pert_flow.W[Wi] = flow_data.W[Wi];
                    pertArea[Si] = area[Si];

                    pert_flow.W[Wi] = flow_data.W[Wi] - dh_W;
                    pertArea[Si] = area[Si] + dh_area;

                    WtoQ(flo_opts.gam, pertArea, pert_flow.W, Q);
                    getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                    const double Resi3 = pert_flow.fluxes[Rikp] * pertArea[Ri + 1] - pert_flow.fluxes[Rik] * pertArea[Ri] - Q[Rik];

                    // R4
                    pert_flow.W[Wi] = flow_data.W[Wi];
                    pertArea[Si] = area[Si];

                    pert_flow.W[Wi] = flow_data.W[Wi] - dh_W;
                    pertArea[Si] = area[Si] - dh_area;

                    WtoQ(flo_opts.gam, pertArea, pert_flow.W, Q);
                    getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);

                    const double Resi4 = pert_flow.fluxes[Rikp] * pertArea[Ri + 1] - pert_flow.fluxes[Rik] * pertArea[Ri] - Q[Rik];

                    ddRdWdArea[Rik].insert(Wi, Si) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                 / (4 * dh_W * dh_area);
                }// Si Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    return ddRdWdArea;
}

std::vector <SparseMatrix<double> > evalddRdWdArea(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const class Flow_data<double> &flow_data)
{
    const int n_elem = flo_opts.n_elem;
    const double gam = flo_opts.gam;
    // Allocate ddRdWdArea Sparse Matrix
    std::vector <SparseMatrix<double> > ddRdWdArea(3 * n_elem);
    SparseMatrix<double> ddRidWdArea(3 * n_elem, n_elem + 1);
    ddRidWdArea.reserve(9 * 2); // 9 State Variables and 2 Areas affect the Residual
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdArea[Ri] = ddRidWdArea;
    }

    // Get Jacobians and Fluxes
    std::vector<double> Ap(n_elem * 3 * 3, 0), An(n_elem * 3 * 3, 0);
    if (flo_opts.flux_scheme == 0) {
        dFluxdW_scalard(flo_opts, flow_data.W, Ap, An);
    } else {
        abort();
    }
    // Evaluate dpdW
    std::vector<double> dpdW(3 * n_elem, 0);
    dpdW = evaldpdW(gam, flow_data.W);

    int Wi, Si;
    double val;
    for (int Ri = 1; Ri < n_elem - 1; Ri++) {
        for (int Rk = 0; Rk < 3; Rk++) {
            const int Rik = Ri * 3 + Rk;
            // Positive Jacobian of Left Incoming Flux
            Si = Ri;
            Wi = Ri - 1;
            for (int Wk = 0; Wk < 3; Wk++) {
                const int Wik = Wi * 3 + Wk;
                val = -Ap[Wi * 9 + Rk * 3 + Wk];
                ddRdWdArea[Rik].insert(Wik, Si) = val;
            }
            // Negative Jacobian of Left Incoming Flux
            Si = Ri;
            Wi = Ri;
            for (int Wk = 0; Wk < 3; Wk++) {
                const int Wik = Wi * 3 + Wk;
                val = -An[Wi * 9 + Rk * 3 + Wk];
                if (Rk == 1) val -= -dpdW[Wi * 3 + Wk];
                ddRdWdArea[Rik].insert(Wik, Si) = val;
            }
            // Positive Jacobian of Right Incoming Flux
            Si = Ri + 1;
            Wi = Ri;
            for (int Wk = 0; Wk < 3; Wk++) {
                const int Wik = Wi * 3 + Wk;
                val = Ap[Wi * 9 + Rk * 3 + Wk];
                if (Rk == 1) val -= dpdW[Wi * 3 + Wk];
                ddRdWdArea[Rik].insert(Wik, Si) = val;
            }
            // Negative Jacobian of Right Incoming Flux
            Si = Ri + 1;
            Wi = Ri + 1;
            for (int Wk = 0; Wk < 3; Wk++) {
                const int Wik = Wi * 3 + Wk;
                val = An[Wi * 9 + Rk * 3 + Wk];
                ddRdWdArea[Rik].insert(Wik, Si) = val;
            }
        }// Rk Loop
    }// Ri Loop
    return ddRdWdArea;
}

std::vector <MatrixXd> evalddQdWdW(const double gam, const std::vector<double> &W)
{
    const int n_elem = W.size()/3;
    std::vector <MatrixXd> ddQdWdW(3 * n_elem);
    MatrixXd ddQidWdW(3, 3);
    ddQidWdW.setZero();
    for (int i = 0; i < n_elem; i++) {
        ddQdWdW[i * 3 + 0] = ddQidWdW;
        ddQdWdW[i * 3 + 2] = ddQidWdW;
    }
    for (int i = 0; i < n_elem; i++) {
        const double rho = W[i * 3 + 0];
        const double u = W[i * 3 + 1] / rho;
        ddQidWdW(0, 0) = u * u * (1.0 - gam) / rho;
        ddQidWdW(1, 1) = (1.0 - gam) / rho;
        ddQidWdW(0, 1) = u * (gam - 1.0) / rho;
        ddQidWdW(1, 0) = ddQidWdW(0, 1);
        // (:,3) = 0
        // (3,:) = 0
        ddQdWdW[i * 3 + 1] = ddQidWdW;
    }
    return ddQdWdW;
}
std::vector <MatrixXd> evalddFdWdW(
    const double gam,
    const double rho,
    const double u,
    const double p)
{
    std::vector <MatrixXd> ddFdWdW(3);
    MatrixXd ddfidWdW(3, 3);
    // F1
    ddfidWdW.setZero();
    ddFdWdW[0] = ddfidWdW;
    // F2
    ddfidWdW(0, 0) = -(u * u * (gam - 3.0) / rho);
    ddfidWdW(1, 1) = (3.0 - gam) / rho;
    ddfidWdW(0, 1) = (u * (gam - 3.0) / rho);
    ddfidWdW(1, 0) = ddfidWdW(0, 1);
    ddFdWdW[1] = ddfidWdW;
    // F3
    ddfidWdW(0, 0) = ( 2.0 * p * u * gam
                    + rho * pow(u, 3.0)
                    * (- 2.0 * gam * gam + 5.0 * gam - 3.0) )
                    / (pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 1) = 3.0 * u * (1.0 - gam) / rho;
    ddfidWdW(2, 2) = 0;

    ddfidWdW(0, 1) = (- 2.0 * p * gam
                    + rho * pow(u, 2.0)
                    * (5.0 * gam * gam - 11.0 * gam + 6.0) )
                    / (2.0 * pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 0) = ddfidWdW(0, 1);

    ddfidWdW(0, 2) = -u * gam / rho;
    ddfidWdW(2, 0) = ddfidWdW(0, 2);

    ddfidWdW(1, 2) = gam / rho;
    ddfidWdW(2, 1) = ddfidWdW(1, 2);

    ddFdWdW[2] = ddfidWdW;

    return ddFdWdW;
}

MatrixXd evalddlambdadWdW(
    const double gam,
    const double rho,
    const double rho_u,
    const double e)
{
    const double u = rho_u / rho;
    const double p = get_p(gam, rho, rho_u, e);
    const double c = sqrt(p * gam / rho);

    MatrixXd ddlambdadWdWp(3, 3);

    ddlambdadWdWp(0, 0) = ( 6.0 * p * gam
                          + rho * u * (8.0 * c - u * gam * (gam - 1.0)) )
                          / (16.0 * c * pow(rho, 3.0));
    ddlambdadWdWp(1, 1) = -c * (gam - 1.0) / (4.0 * p);
    ddlambdadWdWp(2, 2) = c * (1.0 - gam) / (8.0 * p * p);

    ddlambdadWdWp(0, 1) = 0.25 * (c * u * (gam - 1.0) / p - 2.0 / rho);
    ddlambdadWdWp(0, 2) = - ( gam * gam * ( 2.0 * p + rho * u * u * (gam - 1.0) ) )
                          / ( 16.0 * pow(c * rho, 3) );

    ddlambdadWdWp(1, 0) = (-4.0 + u * gam * (gam - 1) / c)
                          / (8.0 * rho * rho);
    ddlambdadWdWp(1, 2) = c * u * (gam - 1.0) / (8.0 * p * p);

    ddlambdadWdWp(2, 0) = - gam * (gam - 1.0)
                          / (8.0 * c * rho * rho);
    ddlambdadWdWp(2, 1) = 0.0;

    return ddlambdadWdWp * eval_dWpdW(gam, rho, rho_u);
}

//void evalddFluxdWdW_Scalar_FD(
//    std::vector <MatrixXd> &ddFluxdWdW1,
//    int ii, int jj,
//    std::vector<double> W)
//{
//    int Fluxik;
//    int Wi1, Wi2, Wi, Wj;
//    double h = 1e-3;
//    std::vector<double> Q(3 * n_elem, 0);
//
//    MatrixXd dummy(3,3);
//    dummy.setZero();
//    for (int k = 0; k < 3 * (n_elem + 1); k++) {
//        ddFluxdWdW1[k] = dummy;
//    }
//    for (int Fluxi = 1; Fluxi < n_elem; Fluxi++) {
//        Wi1 = Fluxi + ii;
//        Wi2 = Fluxi + jj;
//        for (int Fluxk = 0; Fluxk < 3; Fluxk++) {
//            Fluxik = Fluxi * 3 + Fluxk;
//            for (int row = 0; row < 3; row++) {
//                for (int col = 0; col < 3; col++) {
//                    Wi = Wi1 * 3 + row;
//                    Wj = Wi2 * 3 + col;
//                    const double pertWi = flow_data.W[Wi] * h;
//                    const double pertWj = flow_data.W[Wj] * h;
//
//                    for (int m = 0; m < 3 * n_elem; m++) {
//                        pert_flow.W[m] = flow_data.W[m];
//                    }
//
//                    if (row != col || ii!=jj) {
//                        // ddFluxdWdW1 = ddFlux(i)/dW(i-1)dW(i-1)
//                        // R1
//                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
//                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux1 = pert_flow.fluxes[Fluxik];
//
//                        // R2
//                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
//                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux2 = pert_flow.fluxes[Fluxik];
//
//                        // R3
//                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
//                        pert_flow.W[Wj] = flow_data.W[Wj] + pertWj;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux3 = pert_flow.fluxes[Fluxik];
//
//                        // R4
//                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
//                        pert_flow.W[Wj] = flow_data.W[Wj] - pertWj;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux4 = pert_flow.fluxes[Fluxik];
//
//                        ddFluxdWdW1[Fluxik](row, col) =
//                                    (Flux1 - Flux2 - Flux3 + Flux4) / (4.0 * pertWi * pertWj);
//
//                        // Symmetry
//                        //ddFluxdWdW1[Fluxik](col, row) = ddFluxdWdW1[Fluxik](row, col);
//                        // Reset pert_flow.W
//                        pert_flow.W[Wi] = flow_data.W[Wi];
//                        pert_flow.W[Wj] = flow_data.W[Wj];
//                    } else {
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux0 = pert_flow.fluxes[Fluxik];
//
//                        // R1
//                        pert_flow.W[Wi] = flow_data.W[Wi] + 2.0 * pertWi;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux1 = pert_flow.fluxes[Fluxik];
//
//                        // R2
//                        pert_flow.W[Wi] = flow_data.W[Wi] + pertWi;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux2 = pert_flow.fluxes[Fluxik];
//
//                        // R3
//                        pert_flow.W[Wi] = flow_data.W[Wi] - pertWi;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux3 = pert_flow.fluxes[Fluxik];
//
//                        // R4
//                        pert_flow.W[Wi] = flow_data.W[Wi] - 2.0 * pertWi;
//                        getFlux<double>(flo_opts, pert_flow.W, &pert_flow.fluxes);
//                        const double Flux4 = pert_flow.fluxes[Fluxik];
//
//                        ddFluxdWdW1[Fluxik](row, col) =
//                                        (-Flux1 + 16*Flux2 - 30*Flux0 + 16*Flux3 - Flux4)
//                                                     / (12 * pertWi * pertWj);
//                        // Reset pert_flow.W
//                        pert_flow.W[Wi] = flow_data.W[Wi];
//                        pert_flow.W[Wj] = flow_data.W[Wj];
//                    }
//                } // col loop
//            } // row loop
//        }// Fluxk
//    }//Fluxi
//
//}
void evalddFluxdWdW_scalard(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3)
{
    const int n_elem = flo_opts.n_elem;
    const double gam = flo_opts.gam;
    const double eps = flo_opts.scalar_d_eps;
    std::vector <MatrixXd> ddFdWdW(3 * n_elem);
    std::vector <MatrixXd> ddFidWdW(3);

    std::vector <Vector3d> dlambdadW(n_elem);

    std::vector <MatrixXd> ddlambdadWdW(n_elem);
    for (int Wi = 0; Wi < n_elem; Wi++) {
        const double rho = W[Wi*3+0];
        const double u = W[Wi*3+1]/rho;
        const double p = get_p(gam, W[Wi*3+0], W[Wi*3+1], W[Wi*3+2]);
        // Evaluate Convective Hessian
        ddFidWdW = evalddFdWdW(gam, rho, u, p);
        for (int k = 0; k < 3; k++) {
            ddFdWdW[Wi * 3 + k] = ddFidWdW[k];
        }
        // Evaluate dlambdadW
        dlambdadW[Wi] = Map<Vector3d>(evaldlambdadW(gam, W[Wi*3+0],W[Wi*3+1],W[Wi*3+2]).data());

        // Evaluate ddlambdadWdW
        ddlambdadWdW[Wi] = evalddlambdadWdW(gam, W[Wi*3+0],W[Wi*3+1],W[Wi*3+2]);
    }

    int Fluxkim, Fluxki;
    MatrixXd ddFluxidWdW(3, 3);
    for (int Fluxi = 1; Fluxi < n_elem; Fluxi++) {
        for (int Fluxk = 0; Fluxk < 3; Fluxk++) {
            Fluxkim = (Fluxi - 1) * 3 + Fluxk;
            Fluxki = Fluxi * 3 + Fluxk;
            // ddFluxdWdW1 = F(i+1/2) / W(i, i)
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxkim](row, col)
                        - 0.5 * eps * ddlambdadWdW[Fluxi - 1](row, col)
                        * (W[Fluxki] - W[Fluxkim]);

                    if (row == Fluxk)
                    {
                        ddFluxidWdW(row, col) += 0.5 * eps * dlambdadW[Fluxi - 1](col);
                    }
                    if (col == Fluxk)
                    {
                        ddFluxidWdW(row, col) += 0.5 * eps * dlambdadW[Fluxi - 1](row);
                    }
                }
            }
            ddFluxdWdW1[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW2 = F(i+1/2) / W(i + 1, i + 1)
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxki](row, col)
                        - 0.5 * eps * ddlambdadWdW[Fluxi](row, col)
                        * (W[Fluxki] - W[Fluxkim]);
                    if (row == Fluxk)
                    {
                        ddFluxidWdW(row, col) -= 0.5 * eps * dlambdadW[Fluxi](col);
                    }
                    if (col == Fluxk)
                    {
                        ddFluxidWdW(row, col) -= 0.5 * eps * dlambdadW[Fluxi](row);
                    }
                }
            }
            ddFluxdWdW2[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW3 = F(i+1/2) / W(i, i + 1)
            ddFluxidWdW.setZero();
            for (int row = 0; row < 3; row++) {
                ddFluxidWdW(row, Fluxk) -= 0.5 * eps * dlambdadW[Fluxi - 1](row);
            }
            for (int col = 0; col < 3; col++) {
                // Important to have += since one of the entries is non-zero
                ddFluxidWdW(Fluxk, col) += 0.5 * eps * dlambdadW[Fluxi](col);
            }
            ddFluxdWdW3[Fluxki] = ddFluxidWdW;
        }
    }

}

void evalddFluxdWdW(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3)
{
    if (flo_opts.flux_scheme == 0) {
        evalddFluxdWdW_scalard(flo_opts, W, ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3);
    } else {
        abort();
    }

}
