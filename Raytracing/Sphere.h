#pragma once

#include "Vector.h"
#include "Ray.h"

class Sphere {
	Vector rho;
	Vector O;
	double R; 
	bool mirror;
	bool transparent;
	double refractiveIndex;

public:
	Sphere(const Vector& O, double R, const Vector& rho, bool mirror = false, bool transparent = false, double refractiveIndex = 1.5);

	double intersect(const Ray& r, Vector& P, Vector& N);

	bool is_mirror();
	bool is_transparent();
	double get_refractiveIndex();
	Vector get_albedo();
	double get_rayon();
};