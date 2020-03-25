#include "Sphere.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

Sphere::Sphere(
	const Vector& center, double radius,
	const Vector& albedo, Material material,
	double ks, double phongExponent,
	double refractiveIndex) :
	Object(albedo, material, ks, phongExponent, refractiveIndex),
	center(center),
	radius(radius)
{};


bool Sphere::intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
	Vector rayDirection = r.get_direction();
	Vector rayOrigin = r.get_origin();

	double a = 1;
	double b = 2 * dot(rayDirection, rayOrigin - center);
	double c = (rayOrigin - center).getNorm2() - radius * radius;

	double delta = b * b - 4 * a * c;
	if (delta < 0) return false;

	double sqrtDelta = sqrt(delta);
	double t1 = (-b - sqrtDelta) / (2 * a);
	double t2 = (-b + sqrtDelta) / (2 * a);
	if (t2 < 0) return false;

	t = t1 < 0 ? t2 : t1;
	P = rayOrigin + t * rayDirection;
	N = P - center; N.normalize();

	return true;
}

// Getters
double Sphere::get_radius() const
{
	return radius;
}

Vector Sphere::get_center() const
{
	return center;
}
