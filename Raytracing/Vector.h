#pragma once
class Vector {

public:
	double coords[3];
	
	Vector(double x = 0, double y = 0, double z = 0);

	double operator[](int i) const;
	double& operator[](int i);

	double getNorm2() const;
	double getNorm() const;
	void normalize();
	Vector& operator+=(const Vector& B);
};

Vector operator+(const Vector& A, const Vector& B);
Vector operator-(const Vector& A, const Vector& B);
Vector operator-(const Vector& A);
Vector operator*(double a, const Vector& B);
Vector operator*(const Vector& A, double b);
Vector operator/(const Vector& A, double b);
Vector operator*(const Vector& A, const Vector& B);

Vector randomCos(const Vector& N);
Vector cross(const Vector& A, const Vector& B);

double dot(const Vector& A, const Vector& B);