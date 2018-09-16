#ifndef mesh_h
#define mesh_h

#include <vector>
class Mesh {
public:
	double x_0, x_n;
	int n_elem, n_face;
    std::vector<double> x, dx, area;

	Mesh(int n_elements, double a, double b);
	std::vector<double> set_area(std::vector<double> geom, int param);
};
#endif
