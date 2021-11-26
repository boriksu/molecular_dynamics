#include "md.h"

/* сила, производная от потенциала*/
double force_LD(double r)
{
	double x;

	//cout << "HELLO \n";

	if (r > RCUT)
		return 0;

	if (r < RMIN)
		return force_LD(RMIN);

	x = SIGMA / r;
	//return -48 * EPS / SIGMA * (pow(x, 13) - 0.5 * pow(x, 7));
	return -4 * EPS / SIGMA * (pow(x, 13) - pow(x, 7));
}

/* потенциал Леннарда-Джонса */
double potential_LD(double r)
{
	double x;

	//cout << "HI \n";

	if (r > RCUT)
		return 0;

	if (r < RMIN)
		return potential_LD(RMIN);
	x = SIGMA / r;
	return 4 * EPS * (pow(x, 12) - 0.5 * pow(x, 6));
}

/* алгоритм Верле */
//vector Verle_R(vector r, vector dr, vector f, double m, double dt)
//{
//	return r + (dr + (f / (2 * m)) * dt * dt);
//}
//
//vector Verle_V(vector dr, double dt)
//{
//	return dr / (2 * dt);
//}