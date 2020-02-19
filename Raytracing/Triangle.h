#pragma once
#include "Geometry.h"
class Triangle :
	public Geometry
{
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, 
		const Vector& albedo, Material material, double ks = 0, double phongExponent = 1000.,
		double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

private:
	Vector vertices[3];
};

