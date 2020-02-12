#pragma once

#include "Vector.h"
#include "Ray.h"

enum class SphereType {
	normal,
	mirror,
	transparent,
	light
};

class Sphere {
private:
	Vector rho;
	Vector O;
	double R; 
	SphereType sphereType;
	double refractiveIndex;
	double phongExponent;
	double ks;

public:
	Sphere(const Vector& O, double R, const Vector& rho, SphereType sphereType = SphereType::normal, 
		double ks = 0, double phongExponent = 1000., double refractiveIndex = 1.5);

	double intersect(const Ray& r, Vector& P, Vector& N);

	SphereType get_sphereType();
	double get_refractiveIndex();
	Vector get_albedo();
	double get_rayon();
	Vector get_center();
	double get_phongExponent();
	double get_ks();
};