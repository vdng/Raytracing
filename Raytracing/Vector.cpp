#include "Vector.h"
#include <math.h>       /* sqrt */

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

double dot(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}