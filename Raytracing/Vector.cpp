#include "Vector.h"
#define _USE_MATH_DEFINES
#include <math.h>       /* sqrt */
#include <random>

std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0, 1);

Vector::Vector(double x, double y, double z) {
	coords[0] = x;
	coords[1] = y;
	coords[2] = z;
}

double Vector::operator[](int i) const {
	return coords[i];
}

double& Vector::operator[](int i) {
	return coords[i];
}

double Vector::getNorm2() const {
	return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
}

double Vector::getNorm() const {
	return sqrt(getNorm2());
}

void Vector::normalize() {
	double n = sqrt(getNorm2());
	coords[0] /= n;
	coords[1] /= n;
	coords[2] /= n;
}

Vector& Vector::operator+=(const Vector& B)
{
	coords[0] += B[0];
	coords[1] += B[1];
	coords[2] += B[2];
	return *this;
}

Vector randomCos(const Vector& N) {
	Vector T1;
	if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2]))
		T1 = Vector(0, -N[2], N[1]);
	else if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2]))
		T1 = Vector(-N[2], 0, N[1]);
	else
		T1 = Vector(-N[1], N[0], 0);
	T1.normalize();
	Vector T2 = cross(N, T1);
	T2.normalize();

	double r1 = distrib(engine);
	double r2 = distrib(engine);

	double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
	double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
	double z = sqrt(r2);

	Vector R = x * T1 + y * T2 + z * N;
	return R;
}

Vector operator+(const Vector& A, const Vector& B) {
	return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
};

Vector operator-(const Vector& A) {
	return Vector(-A[0], -A[1], -A[2]);
}

Vector operator-(const Vector& A, const Vector& B) {
	return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}

Vector operator*(double a, const Vector& B) {
	return Vector(a * B[0], a * B[1], a * B[2]);
}

Vector operator*(const Vector& A, double b) {
	return Vector(A[0] * b, A[1] * b, A[2] * b);
}

Vector operator/(const Vector& A, double b) {
	return Vector(A[0] / b, A[1] / b, A[2] / b);
}

Vector operator*(const Vector& A, const Vector& B)
{
	return Vector(A[0] * B[0], A[1] * B[1], A[2] * B[2]);
}

Vector cross(const Vector& A, const Vector& B) {
	return Vector(A[1] * B[2] - B[1] * A[2], A[2] * B[0] - B[2] * A[0], A[0] * B[1] - B[0] * A[1]);
}

double dot(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}