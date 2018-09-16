#include <math.h>
#include <vector>
#include "mesh.h"
#include "globals.h"
#include<Eigen/Dense>

using namespace Eigen;
Mesh::Mesh(int n_elements, double x_start = 0, double x_end = 1) {
	x_0 = x_start;
	x_n = x_end;
	n_elem = n_elements;
	n_face = n_elem + 1;

    x.resize(n_face);
    area.resize(n_face);

    dx.resize(n_elem);
	// Build uniform x distribution at interface/flux locations
    double dx_uniform = (x_n - x_0) / n_elem;
    for(int i = 0; i < n_face; i++) {
        x[i] = i*dx_uniform;
	}
    for(int i = 0; i < n_elem; i++) {
        dx[i] = x[i+1] - x[i];
	}
}

std::vector<double> sine_parametrization(std::vector<double> geom, std::vector<double> x);
std::vector <double> eval_spline(std::vector <double> control_pts, std::vector <double> x);

std::vector<double> Mesh::set_area(std::vector<double> geom, int param) {
    int n_geom = geom.size();
    int n_area = area.size();
    if(param == 0) {
		if(n_geom != n_area) abort();
		area = geom;
        for(int i = 0; i < n_face; i++) {
            area[i] = geom[i - 1];
		}

    } else if(param == 1) {
		if(n_geom != 3) abort();
        area = sine_parametrization(geom, x);

    } else if(param == 2) {
		if(n_geom != n_control_pts) abort();

        std::vector<double> control_pts(n_control_pts);
        for(int i = 1; i < n_control_pts - 1; i++) {
            control_pts[i] = geom[i - 1];
        }
        control_pts[0] = 1;
        control_pts[n_control_pts - 1] = 1;
        area = eval_spline(control_pts, x);
    }
    return area;
}

std::vector<double> sine_parametrization(std::vector<double> geom, std::vector<double> x) {
    int n_area = x.size();
	std::vector<double> area(n_area);
	// Sets the area as 1-a*(sin(pi*x^b))^c
    for(int i = 0; i < n_area; i++) {
        area[i] = 1 - geom[0] * pow(sin(M_PI * pow(x[i], geom[1])), geom[2]);
    }
    area[0] = 1;
    area[n_area-1] = 1;

    return area;
}

std::vector <double> get_knots(int n_knots);
double eval_bij(double x, int i, int j, std::vector <double> t);

std::vector <double> eval_spline(std::vector <double> control_pts, std::vector <double> x) {
    // Get Knot Vector
    int n_knots = n_control_pts + spline_degree + 1;
    std::vector <double> knots(n_knots);
    knots = get_knots(n_knots);

    // Define Area at i+1/2
    int n_face = x.size();
    std::vector <double> area(n_face, 0);
    for(int i = 0; i < n_face; i++) {
        for(int j = 0; j < n_control_pts; j++) {
            area[i] += control_pts[j] * eval_bij(x[i], j, spline_degree, knots);
        }
    }
    area[0] = 1;
    area[n_face-1] = 1;

    return area;
}

MatrixXd eval_dArea_dControl(std::vector <double> x, std::vector <double> dx) {
    MatrixXd dArea_dControl(nx + 1, nDesVar);

    // Get Knot Vector
    int n_knots = n_control_pts + spline_degree + 1;
    std::vector <double> knots(n_knots);
    knots = get_knots(n_knots);

    dArea_dControl.setZero();
    for(int i = 0; i < nx + 1; i++) {
        for(int j = 1; j < n_control_pts - 1; j++) {// Not including the inlet/outlet 
            dArea_dControl(i, j - 1) += eval_bij(x[i], j, spline_degree, knots);
        }
    }
    return dArea_dControl;
}

std::vector <double> fit_bspline_to_area(std::vector <double> x, std::vector <double> area, int n_control_pts, int spline_degree) {
    // Get Knot Vector
    int n_knots = n_control_pts + spline_degree + 1;
    std::vector <double> knots(n_knots);
    knots = get_knots(n_knots);

	int n_area = area.size();
    VectorXd control_pts_solution(n_control_pts);
	VectorXd rhs(n_area);
    MatrixXd A(n_area, n_control_pts);
    
    A.setZero();
    for(int i = 0; i < n_area; i++) {
        rhs(i) = area[i];
        for(int j = 0; j < n_control_pts; j++) {
            A(i, j) = eval_bij(x[i], j, spline_degree, knots);
        }
    }
    control_pts_solution = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(rhs);
    
    std::vector <double> control_pts(n_control_pts);
    for(int j = 0; j < n_control_pts; j++) {
        control_pts[j] = control_pts_solution(j);
    }
    control_pts[0] = area[0];
    control_pts[n_control_pts - 1] = area[nx];

    return control_pts;
}

double eval_bij(double x, int i, int j, std::vector <double> t) {
    if(j==0) {
        if(t[i] <= x && x < t[i+1]) return 1;
        else return 0;
    }
    double h = eval_bij(x, i,   j-1, t);
    double k = eval_bij(x, i+1, j-1, t);
    double bij = 0;
    if(h!=0) bij += (x        - t[i]) / (t[i+j]   - t[i]  ) * h;
    if(k!=0) bij += (t[i+j+1] - x   ) / (t[i+j+1] - t[i+1]) * k;
    return bij;
}

std::vector <double> get_knots(int n_knots) {
    std::vector <double> knots(n_knots);
    int n_outer = 2 * (spline_degree + 1);
    int n_inner = n_knots - n_outer;
    double eps = 2e-15; // Allow Spline Definition at End Point
    // Clamped Open-Ended
    for(int iknot = 0; iknot < spline_degree + 1; iknot++)
    {
        knots[iknot] = a_geom;
        knots[n_knots - iknot - 1] = b_geom + eps;
    }
    // Uniform Knot Vector
    double knot_dx = (b_geom + eps - a_geom) / (n_inner + 1);
    for(int iknot = 1; iknot < n_inner+1; iknot++)
    {
        knots[iknot + spline_degree] = iknot * knot_dx;
    }
    return knots;
}
