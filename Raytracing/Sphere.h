#pragma once

#include "Vector.h"
#include "Ray.h"
#include "Geometry.h"

class Sphere :
	public Geometry
{
public:
	Sphere(const Vector& center, double radius, 
		const Vector& albedo, Material material, double ks = 0, double phongExponent = 1000., 
		double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

	// Getters
	double get_radius() const;
	Vector get_center() const;

private:
	Vector center;
	double radius;

};