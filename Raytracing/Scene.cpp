#include "Scene.h"
#include "Vector.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <random>

Scene::Scene(): 
	intensiteL(3 * 1e9), 
	fov(M_PI / 3.), 
	refractiveIndex(1.), 
	light(Vector(-10, 20, 40), 2, Vector(1., 1., 1.)) 
{};

bool Scene::intersect(const Ray& r, Vector& P, Vector& N, int& idx)
{
	bool has_intersect = false;	
	Vector Plocal, Nlocal;
	double t = std::numeric_limits<double>::max();
	for (size_t i = 0; i < spheres.size(); i++)
	{
		double tlocal = spheres[i].intersect(r, Plocal, Nlocal);
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

	if (numRebound == 0)
	{
		return Vector(0, 0, 0);
	}

	Vector P, N, I(0., 0., 0.);
	int idx;

	bool has_intersect = intersect(r, P, N, idx);

	if (has_intersect)
	{
		// Lumière étendue bruitée
		//if (idx == 0) 
		//{
		//	return light.get_albedo() * intensiteL / (4 * M_PI * light.get_rayon() * light.get_rayon());
		//}

		if (spheres[idx].is_mirror())
		{
			Vector R = r.u - 2 * dot(r.u, N) * N;
			Ray rray(P + 1e-4 * N, R);
			I = getColor(rray, numRebound - 1);
		}

		else if (spheres[idx].is_transparent())
		{
			double n1 = refractiveIndex;
			double n2 = spheres[idx].get_refractiveIndex();
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
				Ray tray(P - 1e-4 * N_transparent, R);
				I = getColor(tray, numRebound - 1);
			}
		}
		
		else  // Ni mirroir, ni transparent
		{	
			// Contribution de l'éclairage directe

			//Vector PL = light - P;
			//double distlight2 = PL.getNorm2();
			//PL.normalize();

			//Vector P_prime, N_prime;
			//Ray r_prime(P + 1e-12 * N, PL);

			//int idx_prime;
			//bool has_intersect_prime = intersect(r_prime, P_prime, N_prime, idx_prime);

			//if (has_intersect_prime && distlight2 > (P - P_prime).getNorm2())
			//	I = Vector(0,0,0);
			//else
			//{
			//	Vector PL = light - P;
			//	PL.normalize();
			//	I = intensiteL / M_PI * spheres[idx].get_albedo() * (std::max(0., dot(N, PL)) / (light - P).getNorm2());
			//}

			Vector axeOP = (P - light.get_center()); axeOP.normalize();
			Vector randomDirection = randomCos(axeOP);
			Vector randomPoint = randomDirection * light.get_rayon() + light.get_center();
			Vector wi = (randomPoint - P); wi.normalize();
			double d_light2 = (randomPoint - P).getNorm2();

			Ray r_light(P + 1e-4 * N, wi);
			Vector P_light, N_light;
			int sphere_idx_light;
			bool hasIntersectLight = intersect(r_light, P_light, N_light, sphere_idx_light);

			if (hasIntersectLight && d_light2 * 0.99 > (P - P_light).getNorm2())
			{
				I = Vector(0., 0., 0.);
			}
			else
			{
				I = intensiteL / (4 * M_PI * d_light2) * std::max(0., dot(N, wi) / dot(axeOP, randomDirection)) * spheres[idx].get_albedo();
			}

			// Contribution indirecte
			randomDirection = randomCos(N);
			Ray randomRay(P + 1e-4 * N, randomDirection);

			I += getColor(randomRay, numRebound - 1) * spheres[idx].get_albedo();

		}
	}
	return I;
}

void Scene::set_camera(Vector C)
{
	camera = C;
}

void Scene::set_fov(double f)
{
	fov = f;
}

void Scene::addSphere(const Sphere& sphere)
{
	spheres.push_back(sphere);
}

void Scene::set_light(Sphere L)
{
	light = L;
}

void Scene::set_intensiteL(double intensite)
{
	intensiteL = intensite;
}

void Scene::set_refractiveIndex(double refraction)
{
	refractiveIndex = refraction;
}

Vector Scene::get_camera()
{
	return camera;
}

double Scene::get_fov()
{
	return fov;
}


Sphere Scene::operator[](int i) const {
	return spheres[i];
}

Sphere& Scene::operator[](int i) {
	return spheres[i];
}