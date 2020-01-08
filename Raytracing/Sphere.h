#pragma once

#include "Vector.h"
#include "Ray.h"

class Sphere {
	Vector rho;
	Vector O;
	double R; 

public:
	Sphere(const Vector& O, double R, const Vector& rho);

	double intersect(const Ray& r, Vector& P, Vector& N);
	Vector intensity(const Ray& r, Vector& P, Vector& N, const Vector& L, const double intensiteL);
};