#include "Scene.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

Scene::Scene(std::vector<Sphere*> spheres, Vector camera, double fov, Vector light, double intensiteL) : 
	spheres(spheres), camera(camera), fov(fov), light(light), intensiteL(intensiteL) {};

bool Scene::intersect(const Ray& r, Vector& P, Vector& N, int& idx)
{
	bool has_intersect = false;	
	Vector Plocal, Nlocal;
	double t = std::numeric_limits<double>::max();
	for (size_t i = 0; i < spheres.size(); i++)
	{
		double tlocal = spheres[i]->intersect(r, Plocal, Nlocal);
		if (0 <= tlocal && tlocal < t)
		{
			has_intersect = true;
			idx = i;
			t = tlocal;
			P = Plocal;
			N = Nlocal;
		}
	}
	return has_intersect;
}

Vector Scene::getColor(const Ray& r) {
	Vector P, N;
	int idx;

	bool has_intersect = intersect(r, P, N, idx);

	Vector PL = light - P;
	double distlight = PL.getNorm();
	PL.normalize();

	Vector P_prime, N_prime, I;
	Ray r_prime(P + 1e-12 * N, PL);
	
	int idx_prime;
	bool has_intersect_prime = intersect(r_prime, P_prime, N_prime, idx_prime);
	if (has_intersect_prime) 
	{
		if (distlight < (P - P_prime).getNorm())
		{
			I = spheres[idx]->intensity(r, P, N, light, intensiteL);
		}
	}
	else
	{
		I = spheres[idx]->intensity(r, P, N, light, intensiteL);
	}
	return I;
}


Sphere Scene::operator[](int i) const {
	return *spheres[i];
}

Sphere& Scene::operator[](int i) {
	return *spheres[i];
}