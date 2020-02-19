#include "Geometry.h"

Geometry::Geometry(const Vector& albedo, Material material, double ks, double phongExponent, double refractiveIndex) :
	albedo(albedo),
	material(material),
	ks(ks),
	phongExponent(phongExponent),
	refractiveIndex(refractiveIndex)
{}

bool Geometry::intersect(const Ray & r, Vector & P, Vector & N, double& t) const
{
	return false;
}

Material Geometry::get_material() const
{
	return material;
}

double Geometry::get_refractiveIndex() const
{
	return refractiveIndex;
}

Vector Geometry::get_albedo() const
{
	return albedo;
}

double Geometry::get_phongExponent() const
{
	return phongExponent;
}

double Geometry::get_ks() const
{
	return ks;
}
