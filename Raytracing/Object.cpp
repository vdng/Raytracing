#include "Object.h"
#include "Triangle.h"
#include <vector>
#include <map>
#include <string>
#include <algorithm>

Object::Object(
    const Vector& albedo, Material material,
    double ks, double phongExponent, double refractiveIndex) :
    albedo(albedo),
    material(material),
    ks(ks),
    phongExponent(phongExponent),
    refractiveIndex(refractiveIndex)
{}

bool Object::intersect(const Ray& r, Vector& P, Vector& N, double& t) const
{
    return false;
}

Material Object::get_material() const
{
    return material;
}

double Object::get_refractiveIndex() const
{
    return refractiveIndex;
}

Vector Object::get_albedo() const
{
    return albedo;
}

double Object::get_phongExponent() const
{
    return phongExponent;
}

double Object::get_ks() const
{
    return ks;
}
