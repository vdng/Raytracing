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
};

Vector operator+(const Vector& A, const Vector& B);
Vector operator-(const Vector& A, const Vector& B);
Vector operator*(double a, const Vector& B);
Vector operator*(const Vector& A, double b);

double dot(const Vector& A, const Vector& B);