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
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, SphereType sphereType = SphereType::normal, 
		double ks = 0, double phongExponent = 1000., double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t);

	// Getters
	SphereType get_sphereType();
	double get_refractiveIndex();
	Vector get_albedo();
	double get_radius();
	Vector get_center();
	double get_phongExponent();
	double get_ks();

private:
	Vector albedo;
	Vector center;
	double radius;
	SphereType sphereType;
	double refractiveIndex;
	double phongExponent;
	double ks;

};