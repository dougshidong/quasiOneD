#ifndef optimizer_h
#define optimizer_h

void design(int nx, int descentType, int gradientType, int fitnessFun,
	      std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar);

#endif
