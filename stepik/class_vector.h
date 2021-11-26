#pragma once
#ifndef class_vector
#define class_vector

#include <iostream>
using namespace std;
class vector
{
private:
	double X, Y, Z;
public:
	vector()
	{
		X = 0.0;
		Y = 0.0;
		Z = 0.0;
	};
	vector(double x, double y, double z) {

		X = x;
		Y = y;
		Z = z;
	};

	~vector()
	{
		//cout << "end vector" << endl;
	}

	void print() { cout << "X = " << X << " Y = " << Y << " Z = " << Z << endl; }

	double getX() { return X; }
	double getY() { return Y; }
	double getZ() { return Z; }
	void set(double x, double y, double z) {
		X = x;
		Y = y;
		Z = z;
	}
	void setX(double x) { X = x; }
	void setY(double y) { Y = y; }
	void setZ(double z) { Z = z; }
	double abs() {
		double value = 0.0;
		value = X * X + Y * Y + Z * Z;
		value = sqrt(value);
		return value;
	}

	double abs2() {									// possible
		double value = 0.0;
		value = X * X + Y * Y + Z * Z;
		//value = sqrt(value);
		return value;
	}

	void normalise()
	{
		double value = X * X + Y * Y + Z * Z;
		value = sqrt(value);
		double inv_length = (1 / value);
		X *= inv_length;
		Y *= inv_length;
		Z *= inv_length;
	}

	//friend vector operator/ (vector A, double t);

	vector operator/ (double t)
	{
		this->X /= t; this->Y /= t; this->Z /= t; return *this;
	}

	vector operator* ( double t)
	{
		this->X *= t; this->Y *= t; this->Z *= t; return *this;
	}

	vector operator+ ( vector B)
	{
		vector T;
		T.X = this->X + B.X; T.Y = this->Y + B.Y; T.Z = this->Z + B.Z;
		return T;
	}

	vector operator- (vector B)
	{
		vector T;
		T.X = this->X - B.X; T.Y = this->Y - B.Y; T.Z = this->Z - B.Z;
		return T;
	}

	vector operator+= (vector B)
	{
		vector T;
		T.X = this->X + B.X; T.Y = this->Y + B.Y; T.Z = this->Z + B.Z;
		return T;
	}

	vector operator-= (vector B)
	{
		vector T;
		T.X = this->X - B.X; T.Y = this->Y - B.Y; T.Z = this->Z - B.Z;
		return T;
	}

	////friend vector operator* (vector A, double t);
	//friend vector operator+ (vector A, vector B);
	//friend vector operator- (vector A, vector B);
	//friend vector operator-= (vector A, vector B);
	//friend vector operator+= (vector A, vector B);

	//void input() { cin >> X >> Y; }
	//void output() { cout << “(“ << X << “, ” << Y << “)”; }
	//void set(float x, float y;) { X = x; Y = y; }
};

//vector operator/ (vector A, double t)
//{
//	A.X /= t; A.Y /= t; A.Z /= t; return A;
//}

//vector operator* (vector A, double t)
//{
//	A.X *= t; A.Y *= t; A.Z *= t; return A;
//}

//vector operator+ (vector A, vector B)
//{
//	vector T;
//	T.X = A.X + B.X; T.Y = A.Y + B.Y; T.Z = A.Z + B.Z;
//	return T;
//}

//vector operator- (vector A, vector B)
//{
//	vector T;
//	T.X = A.X - B.X; T.Y = A.Y - B.Y; T.Z = A.Z - B.Z;
//	return T;
//}

//vector operator-= (vector A, vector B)
//{
//	vector T;
//	T.X = A.X - B.X; T.Y = A.Y - B.Y; T.Z = A.Z - B.Z;
//	return T;
//}
//
//vector operator+= (vector A, vector B)
//{
//	vector T;
//	T.X = A.X + B.X; T.Y = A.Y + B.Y; T.Z = A.Z + B.Z;
//	return T;
//}

#endif