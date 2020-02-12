#pragma once

class Vector {
public:
	//double coords[3];
	
	Vector(double x = 0, double y = 0, double z = 0);

	Vector reflect(const Vector& N) const;
	void normalize();

	double operator[](int i) const;
	double& operator[](int i);
	Vector& operator+=(const Vector& B);

	double getNorm2() const;

private:
	double coords[3];
};

Vector operator+(const Vector& A, const Vector& B);
Vector operator-(const Vector& A, const Vector& B);
Vector operator-(const Vector& A);
Vector operator*(double a, const Vector& B);
Vector operator*(const Vector& A, double b);
Vector operator/(const Vector& A, double b);
Vector operator*(const Vector& A, const Vector& B);

Vector randomCos(const Vector& N);
Vector randomPhong(const Vector& radius, double phongExponent);
Vector cross(const Vector& A, const Vector& B);

double dot(const Vector& A, const Vector& B);