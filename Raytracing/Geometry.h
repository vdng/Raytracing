#pragma once
#include "Ray.h"
#include "Vector.h"

enum class Material
{
	normal,
	mirror,
	transparent,
	light
};

class Geometry
{
public:
	Geometry(const Vector& albedo, Material material, double ks, double phongExponent, double refractiveIndex);
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

	Material get_material() const;
	Vector get_albedo() const;
	double get_phongExponent() const;
	double get_ks() const;
	double get_refractiveIndex() const;

private:
	Vector albedo;
	Material material;
	double refractiveIndex;
	double phongExponent;
	double ks;
};

