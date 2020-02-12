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

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

	// Getters
	SphereType get_sphereType() const;
	double get_refractiveIndex() const;
	Vector get_albedo() const;
	double get_radius() const;
	Vector get_center() const;
	double get_phongExponent() const;
	double get_ks() const;

private:
	Vector albedo;
	Vector center;
	double radius;
	SphereType sphereType;
	double refractiveIndex;
	double phongExponent;
	double ks;

};