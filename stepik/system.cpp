#include "md.h"


/* ������� ��� ������ */
void ssystem::initialize()
{
	srand(std::time(NULL));          
	int K = ceil(pow(N, 1.0 / 3));   // ���������� �������, ��� ������������� �������
	double dLx = Lx / K;
	double dLy = Ly / K;
	double dLz = Lz / K;
	int counter = 0;
	new_frame();


	for (int i = 0; i < N; i++)								/* ������� ��������� ��������� �������� */
	{
		v[i] = vector((rand() % 10000 - 5000) / 10000.0,
			(rand() % 10000 - 5000) / 10000.0,
			(rand() % 10000 - 5000) / 10000.0);
		m[i] = 1;
	}
	calc_stats();											// ���� �����������
	cout << "init T = " << T << endl;
	for (int i = 0; i < N; i++)								// ��������������� ��������� - ��������� ���������
	{
		v[i] = vector(v[i].getX() * pow(( 1 + ((T0 / T) - 1)),0.5),
			v[i].getY() * pow((1 + ((T0 / T) - 1)), 0.5),
			v[i].getZ() * pow((1 + ((T0 / T) - 1)), 0.5));
	}
	calc_P();												// ������� ������� ������ ���� �������
	for (int i = 0; i < N; i++)								// �������� ������� �� ���������, ����� ����� ���� ��������� �����������
		v[i] -= Psys * (1 / msys);

	for (int i = 0; i < K; i++)								// ������������ ������
		for (int j = 0; j < K; j++)
			for (int k = 0; k < K; k++)
				if (counter < N)
				{
					r[counter].set((i + 1.0 / 2) * dLx, (j + 1.0 / 2) * dLy,
						(k + 1.0 / 2) * dLz);
					dr[counter].set(v[counter].getX() * 2 * dt,		// ������ 2? 6:50
									v[counter].getY() * 2 * dt,
									v[counter].getZ() * 2 * dt);
					//assert(dr[counter].abs() < L_FREE_MOTION);	// ������� ��������� ������� ������� ����������
					print_to_file(counter);
					counter++;
				}
	isInitialised = true;
	calc_stats();
	initial_energy = energy;
}


void ssystem::calc_stats()
{
	double avgV = 0;
	double sqrAvgV = 0;
	double kinEnergy = 0;
	double potEnergy = 0;
	double E = 0;
	const double cor_ratio = sqrt(3 * 3.1416 / 8);
	double ratio = 0;

#pragma omp parallel for reduction (+ : avgV, sqrAvgV, kinEnergy, potEnergy)
	for (int i = 0; i < N; i++)
	{
		if (isInitialised)
		{
			//cout << "\ninside\n";
			avgV += v[i].abs();
			sqrAvgV += v[i].abs2();
			kinEnergy += m[i] * v[i].abs2() / 2;
			for (int j = 0; j < N; j++)
				if (i != j)
				{
					vector dr = NIM(r[i], r[j], Lx, Ly, Lz);
					double rij = dr.abs();
					potEnergy += potential(rij);
				}
		}
		else
			kinEnergy += m[i] * v[i].abs2() / 2;
	}
	avgV /= N;
	sqrAvgV = sqrt(sqrAvgV / N);
	potEnergy /= 2;
	T = kinEnergy / (3.0 * N / 2);
	E = potEnergy + kinEnergy;

	if (isInitialised)					// �� ������� �� ������� ������� �������,
	{									// ���� � 2 � ������, �� �������� ������� ������������ � ��� ������ ����������
		ratio = sqrAvgV / avgV;			// 10:40
		//assert(fabs(ratio - cor_ratio) / cor_ratio < 0.05);
		if (initial_energy != 0)
			assert(fabs((E - initial_energy) / initial_energy) < ENERGY_PERCENT);
	}
	averV = avgV;
	sqrAverV = sqrAvgV;
	kinetic_energy = kinEnergy;
	potential_energy = potEnergy;
	energy = E;
}

void ssystem::measure_all()
{
	calc_stats();
	cout << "avg v: " << averV << endl;
	calc_P();
	cout << "Mass center impulse = " << Psys.abs() << endl;
	for (int i = 0; i < N; i++)
		v[i] -= Psys * (1 / msys);
	cout << "Kinetic Energy = " << kinetic_energy << endl;
	cout << "Full Energy = " << energy << endl;
}

void ssystem::output_data()
{
	cout << "\nModelling is finished, measuring final data\n" << endl;
	measure_all();
}


/* ������� ��� ������ .xmol - format file */
void ssystem::new_frame()
{
	out << N << endl;
	out << "****** time = " << time << " *******" << endl;
	if (time > termalization_time)
	{
		term_out << N << endl;
		term_out << "****** time = " << time << " *******" << endl;

	}
}

void ssystem::print_to_file(int i)
{
	out << "Ar " << setw(10) << r[i].getX() << " " << setw(10)
		<< r[i].getY() << " " << setw(10) << r[i].getZ() << " " << setw(10)
		<< v[i].getX() << " " << setw(10) << v[i].getY() << " " << setw(10)
		<< v[i].getZ() << setw(10) << v[i].abs() << endl;
	if (time > termalization_time)
	{
		term_out << "Ar " << setw(10) << r[i].getX() << " " << setw(10)
			<< r[i].getY() << " " << setw(10) << r[i].getZ() << " "
			<< setw(10) << v[i].getX() << " " << setw(10)
			<< v[i].getY() << " " << setw(10) << v[i].getZ() << endl;
	}
}

/* ������� ������� �������� ������ ���� ������� */
void ssystem::calc_P() {
	Psys = vector(0, 0, 0);
	if (msys == 0)
		for (int i = 0; i < N; i++)
		{
			Psys += v[i] * m[i];
			msys += m[i];
		}
	else
		for (int i = 0; i < N; i++)
			Psys += v[i] * m[i];
}

/* ������� ������� ��� */
void ssystem::calc_forces()
{
	double ff;
	double _rij;

#pragma omp parallel for private (ff, _rij)
	for (int i = 0; i < N; i++)
	{
		f[i] = vector(0, 0, 0);
		for (int j = 0; j < N; j++)
		{
			//cout << "f[" << i << "] : ";
			//f[i].print();
			if (i != j)
			{
				vector rij = NIM(r[i], r[j], Lx, Ly, Lz);
				//cout << "rij ";
				//rij.print();
				_rij = rij.abs();
				//cout << "rij abs" << _rij << endl;
				ff = force(_rij);
				if (ff != 0)
					cout << "force " << ff << " | ";
				vector dr = rij;
				//cout << "dr ";
				//dr.print();
				dr.normalise();
				//cout << "dr norm ";
				//dr.print();
				f[i] += dr * ff;
				//cout << "after f[" << i << "] : ";
				//f[i].print();
			}
			//if (j == 5)
			//	exit(0);
		}
		//cout << "F [" << i << "] = " << setw(10) << f[i].getX() << " " << setw(10)
		//	<< f[i].getY() << " " << setw(10) << f[i].getZ() << endl;
		//exit(0);
	}
}


/* ������� �������������� ��������� �������� */
void ssystem::integrate()
{
	vector rBuf;
	time += dt;
	new_frame();
#pragma omp parallel for private (rRuf)
	for (int i = 0; i < N; i++)
	{
		rBuf = r[i];									// ���������� ���������� � ������� ������ �������
		r[i] = Verle_R(r[i], dr[i], f[i], m[i], dt);		// ��������� ����� ���������� �� �����
		dr[i] = r[i] - rBuf;

		//assert(dr[i].abs() < L_FREE_MOTION);

		v[i] = Verle_V(dr[i], dt);
		PBC(r[i], Lx, Ly, Lz);							// ������������� ��������� �������, ���� ������� �������� �� �������
	}
	for (int i = 0; i < N; i++)
		print_to_file(i);

}

vector ssystem::Verle_R(vector r, vector dr, vector f, double m, double dt)
{
	return r + (dr + (f / (2 * m)) * dt * dt);
}

vector ssystem::Verle_V(vector dr, double dt)
{
	return dr / (2 * dt);
}

/* ����� ��������� ������� */
vector ssystem::NIM(vector r1, vector r2, double Dx, double Dy, double Dz) const
{
	double x = -(r1.getX() - r2.getX());
	double y = -(r1.getY() - r2.getY());
	double z = -(r1.getZ() - r2.getZ());

	vector dist = vector(0, 0, 0);

	dist.setX(NIM_fix(x, Dx));
	dist.setY(NIM_fix(y, Dy));
	dist.setZ(NIM_fix(z, Dz));

	//assert(dist.abs() < vector(Dx, Dy, Dz).abs());

	return dist;
}

double ssystem::NIM_fix(double coord, double l) const				// �������� �������
{
	if (coord >= l / 2.0)
		coord = l - coord;
	else if (coord <= l / 2.0)
		coord = coord + l;

	return coord;
}

void ssystem::PBC(vector &dist, double Dx, double Dy, double Dz) const
{
	double aDx = 0;
	double aDy = 0;
	double aDz = 0;

	dist.setX(correctCoord(dist.getX(), aDx, Dx));			// ������� ����������, ����� � ����� �������
	dist.setY(correctCoord(dist.getY(), aDy, Dy));
	dist.setZ(correctCoord(dist.getZ(), aDz, Dz));
}

double ssystem::correctCoord(double coord, double left, double right) const
{
	double l = right - left;			// ����� ����� ������ ������������� ����� ����� ���� ���
	double d;

	if (coord >= right)
	{
		coord = coord - l;
		//d = coord - left;
		//coord = coord - l * floor(d / l);
		// ���������������� ����� - �� �����, � ����� �����, ����� ������� ������ � ������
	}
	else if (coord < left)
	{
		d = left - coord;
		coord = right = l * (d / l);
		//coord = roght - l * (d / l - floor(d / l));
		// ���������������� �����
	}

	return coord;
}