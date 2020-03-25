#pragma once

#include "Vector.h"
#include "Ray.h"
#include "Object.h"

class Sphere :
	public Object
{
public:
	Sphere(
		const Vector& center, double radius,
		const Vector& albedo = Vector(0.5, 0.5, 0.5), Material material = Material::normal,
		double ks = 0, double phongExponent = 1000.,
		double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

	// Getters
	double get_radius() const;
	Vector get_center() const;

private:
	Vector center;
	double radius;
};