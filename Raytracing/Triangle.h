#pragma once
#include "Object.h"
class Triangle :
	public Object
{
public:
	Triangle(
		const Vector& A, const Vector& B, const Vector& C,
		const Vector& albedo = Vector(0.5, 0.5, 0.5), Material material = Material::normal, 
		double ks = 0, double phongExponent = 1000.,
		double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

private:
	Vector vertices[3];
};

