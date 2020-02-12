#include "Sphere.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

Sphere::Sphere(const Vector& O, double R, const Vector& rho, SphereType sphereType, double ks, double phongExponent, double n) :
	O(O),
	R(R),
	rho(rho),
	sphereType(sphereType),
	ks(ks),
	phongExponent(phongExponent),
	refractiveIndex(n)
{};


double Sphere::intersect(const Ray& r, Vector& P, Vector& N) {
	double a = 1;
	double b = 2 * dot(r.u, r.C - O);
	double c = (r.C - O).getNorm2() - R * R;

	double delta = b * b - 4 * a * c;
	if (delta < 0) return -1;

	double sqrtDelta = sqrt(delta);
	double t1 = (-b - sqrtDelta) / (2 * a);
	double t2 = (-b + sqrtDelta) / (2 * a);

	if (t2 < 0) return -1;
	double t = t1 < 0 ? t2 : t1;

	P = r.C + t * r.u;
	N = P - O;
	N.normalize();

	return t;
}

SphereType Sphere::get_sphereType()
{
	return sphereType;
}

double Sphere::get_refractiveIndex()
{
	return refractiveIndex;
}

Vector Sphere::get_albedo()
{
	return rho;
}

double Sphere::get_rayon()
{
	return R;
}

Vector Sphere::get_center()
{
	return O;
}

double Sphere::get_phongExponent()
{
	return phongExponent;
}

double Sphere::get_ks()
{
	return ks;
}
