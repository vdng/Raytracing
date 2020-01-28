#include "Scene.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

Scene::Scene(std::vector<Sphere*> spheres, Vector camera, double fov, Vector light, double intensiteL, double refractiveIndex) :
	spheres(spheres), camera(camera), fov(fov), light(light), intensiteL(intensiteL), refractiveIndex(refractiveIndex) {};

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

Vector Scene::getColor(const Ray& r, int numRebound) {
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

	if (has_intersect)
	{
		if (spheres[idx]->is_mirror() && numRebound > 0)
		{
			Vector R = r.u - 2 * dot(r.u, N) * N;
			Ray rray(P + 1e-4 * N, R);
			I = getColor(rray, numRebound - 1);
		}

		else if (spheres[idx]->is_transparent())
		{
			double n1 = refractiveIndex;
			double n2 = spheres[idx]->get_refractiveIndex();
			Vector N_transparent(N);

			if (dot(r.u, N) > 0) // Le rayon sort de la sphère
			{
				n1 = n2;
				n2 = refractiveIndex;
				N_transparent = -N;
			}

			double radical = 1. - (n1 / n2) * (n1 / n2) * (1 - dot(r.u, N_transparent) * dot(r.u, N_transparent));

			if (radical > 0)
			{
				Vector R = (n1 / n2) * (r.u - dot(r.u, N_transparent) * N_transparent) - sqrt(radical) * N_transparent;
				Ray rray(P - 1e-4 * N_transparent, R);
				I = getColor(rray, numRebound);
			}
		}
		
		else
		{
			if (has_intersect_prime && distlight > (P - P_prime).getNorm())
				I = Vector(0,0,0);
			else
				I = spheres[idx]->intensity(r, P, N, light, intensiteL);
		}
	}
	return I;
}


Sphere Scene::operator[](int i) const {
	return *spheres[i];
}

Sphere& Scene::operator[](int i) {
	return *spheres[i];
}