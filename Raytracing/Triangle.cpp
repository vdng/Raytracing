#include "Triangle.h"
#include <cstdlib>

Triangle::Triangle(const Vector& A, const Vector& B, const Vector& C, 
	const Vector& albedo, Material material, double ks, double phongExponent, double refractiveIndex) :
	Geometry(albedo, material, ks, phongExponent, refractiveIndex)
{
	vertices[0] = A;
	vertices[1] = B;
	vertices[2] = C;
}

bool Triangle::intersect(const Ray& r, Vector& P, Vector& N, double& t) const
{
	Vector rayDirection = r.get_direction();
	Vector rayOrigin = r.get_origin();

	N = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]); N.normalize();

	double denom = dot(rayDirection, N);
	if (std::abs(denom) < 1e-12) return false; // rayon parallèle
	
	t = dot(vertices[2] - rayOrigin, N) / denom;
	if (t < 0) return false; // Intersection derrière la caméra

	P = rayOrigin + rayDirection * t;
	double APAB = dot(P - vertices[0], vertices[1] - vertices[0]);
	double ACAB = dot(vertices[2] - vertices[0], vertices[1] - vertices[0]);
	double ABAB = dot(vertices[1] - vertices[0], vertices[1] - vertices[0]);
	double APAC = dot(P - vertices[0], vertices[2] - vertices[0]);
	double ACAC = dot(vertices[2] - vertices[0], vertices[2] - vertices[0]);
	double det = ABAB * ACAC - ACAB * ACAB;
	double beta = (APAB * ACAC - APAC * ACAB) / det;
	double gamma = (ABAB * APAC - ACAB * APAB) / det;
	if (beta < 0 || gamma < 0 || beta + gamma > 1) return false;

	return true;
}
