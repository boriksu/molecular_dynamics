#ifndef mdh
#define mdh

#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <memory.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "class_vector.h"

#define PRINT_VEC(vec) (cout << "(" << vec.getX() << ";" << vec.getY() << ";" << vec.getZ() << ")" << endl);
#define POTENTIAL_ENERGY_PERCENT 0.01
#define ENERGY_PERCENT 0.1

using namespace std;
constexpr double EPS = 1;;
constexpr double SIGMA = 1;;
constexpr double RCUT = 3;;
constexpr double RMIN = 0.00001;;

typedef double (*force_type)(double);
typedef double (*potential_type)(double);
//typedef Vec3<double> vector;
//typedef vector(*alg_R)(vector, vector, vector, double, double);
//typedef vector(*alg_V)(vector, double);

void test();
double force_LD(double r);
double potential_LD(double r);
//vector Verle_R(vector r, vector dr, vector f, double m, double dt);
//vector Verle_V(vector dr, double dt);

class ssystem
{
	public:
		ssystem(force_type Force, potential_type Potential) : force(Force), potential(Potential)
		{
			dr = new vector[N];
			v = new vector[N];
			f = new vector[N];
			r = new vector[N];
			m = new double[N];
			time = 0;
			energy = 0;
			kinetic_energy = 0;
			initial_energy = 0;
			out.open("output.xmol");
			term_out.open("output.001");
			L_FREE_MOTION = (pow(Lx * Ly * Lz, 1.0 / 3) / (2 * N));
		}

		~ssystem()
		{
			delete[] dr;
			delete[] v;
			delete[] f;
			delete[] r;
			delete[] m;
			out.close();
			term_out.close();
		}

		ofstream out;									// вывод в два файла
		ofstream term_out;

		void new_frame();
		void print_to_file(int i);
		void calc_P();
		void calc_forces();
		void integrate();
		void initialize();
		void calc_stats();
		vector NIM(vector r1, vector r2, double Dx, double Dy, double Dz) const;
		double NIM_fix(double coord, double l) const;
		void PBC(vector& dist, double Dx, double Dy, double Dz) const;
		double correctCoord(double coord, double left, double right) const;
		void measure_all();
		void output_data();
		vector Verle_R(vector r, vector dr, vector f, double m, double dt);
		vector Verle_V(vector dr, double dt);

	private:
		double L_FREE_MOTION;
		const int N = 432;								// кол-во частиц
		const double Lx = 276, Ly = 276, Lz = 276;		// размер системы 10е-11 m
		const double dt = 0.001;						// шаг моделировани€
		const double T0 = 1.38 * 111.84;				// начальна€ температура
		//bool v0_rand = true;							// начальна€ скорость частиц есть
		bool isInitialised = false;						// инициализирована ли система
		vector *r;										// координаты частиц
		vector *dr;										// рассто€ние до предыдущих координат
		vector *v;										// скорости
		vector *f;										// силы
		double *m;										// массы
		double energy;									// полна€ энерги€ системы
		double initial_energy = 0;						// начальна€ полна€ энерги€ системы
		double kinetic_energy;
		double potential_energy;
		double averV;									// средн€€ скорость частиц
		double sqrAverV;								// среднеквадратическа€ скорость частиц
		double T;										// температура в данный момент времени
		double msys = 0;								// полна€ масса системы
		vector Psys;									// момент импульса системы
		force_type force;								// указатели на функции
		//alg_R сalc_R;									// указатели на функции
		//alg_V сalc_V;									// указатели на функции
		potential_type potential;						// указатели на функции
		double time;									// текущий момент времени
		double termalization_time = 5;					// врем€ на термализацию системы 10е-12 с
};

#endif