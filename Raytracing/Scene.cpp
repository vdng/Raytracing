#include "Scene.h"
#include "Vector.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <random>

extern std::default_random_engine engine;
extern std::uniform_real_distribution<double> distrib;

Scene::Scene(Sphere* light, double totalIntensity):
	lightIntensity(totalIntensity * 4. * M_PI / (4 * M_PI * M_PI * light->get_radius() * light->get_radius())),
	fov(M_PI / 3.), 
	refractiveIndex(1.), 
	light(light) 
{};

void Scene::addSphere(const Sphere& sphere)
{
	objects.push_back(&sphere);
}

void Scene::addTriangle(const Triangle& triangle)
{
	objects.push_back(&triangle);
}


bool Scene::intersect(const Ray& r, Vector& P, Vector& N, int& idx)
{
	bool has_intersect = false;	
	Vector P_local, N_local;
	double t = std::numeric_limits<double>::max();
	for (size_t i = 0; i < objects.size(); i++)
	{
		double t_local;
		if (objects[i]->intersect(r, P_local, N_local, t_local))
		{
			has_intersect = true;
			if (t_local < t)
			{
				idx = i;
				t = t_local;
				P = P_local;
				N = N_local;
			}
		}
	}
	return has_intersect;
}

Vector Scene::getColor(const Ray& r, int numRebound, bool showLights)
{
	if (numRebound == 0) return Vector(0., 0., 0.);

	Vector P, N, I(0., 0., 0.);
	int idx;

	bool has_intersect = intersect(r, P, N, idx);
	Vector rayDirection = r.get_direction();
	// Vector rayOrigin = r.get_origin();

	if (has_intersect)
	{
		const Geometry* currentObject = objects[idx];

		// Lumière étendue bruitée
		//if (idx == 0) 
		//{
		//	return light.get_albedo() * lightIntensity / (4 * M_PI * light.get_radius() * light.get_radius());
		//}

		Material material = currentObject->get_material();

		if (material == Material::light)
		{
			I = showLights ? (light->get_albedo() * lightIntensity) : Vector(0., 0., 0.);
		}

		else if (material == Material::mirror)
		{
			Vector R = rayDirection.reflect(N);
			Ray rray(P + 1e-4 * N, R);
			I = getColor(rray, numRebound - 1);
		}

		else if (material == Material::transparent)
		{
			double n1 = refractiveIndex;
			double n2 = currentObject->get_refractiveIndex();
			Vector N_transparent(N);

			if (dot(rayDirection, N) > 0) // Le rayon sort de la sphère
			{
				n1 = n2;
				n2 = refractiveIndex;
				N_transparent = -N;
			}

			double radical = 1. - (n1 / n2) * (n1 / n2) * (1 - dot(rayDirection, N_transparent) * dot(rayDirection, N_transparent));

			if (radical > 0)
			{
				Vector R = (n1 / n2) * (rayDirection - dot(rayDirection, N_transparent) * N_transparent) - sqrt(radical) * N_transparent;
				Ray tray(P - 1e-4 * N_transparent, R);
				I = getColor(tray, numRebound - 1);
			}
		}
		
		else // if (material == Material::normal)
		{	
			// Contribution de l'éclairage directe
			Vector axeOP = (P - light->get_center()); axeOP.normalize();
			Vector randomDirection = randomCos(axeOP);
			Vector randomPoint = randomDirection * light->get_radius() + light->get_center();
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
				Vector BRDF = currentObject->get_albedo() / M_PI * (1. - currentObject->get_ks()) + currentObject->get_ks() * phongBRDF(wi, rayDirection, N, currentObject->get_phongExponent()) * currentObject->get_albedo();
				double J = 1. * dot(randomDirection, -wi) / d_light2;
				double proba = dot(axeOP, randomDirection) / (M_PI * light->get_radius() * light->get_radius());
				I = lightIntensity * std::max(0., dot(N, wi)) * J * BRDF / proba;
			}

			// Contribution indirecte

			double p = 1 - currentObject->get_ks();
			bool sampleDiffuse;
			Vector R = rayDirection.reflect(N);

			if (distrib(engine) < p) {
				sampleDiffuse = true;
				randomDirection = randomCos(N);
			}
			else 
			{ 
				sampleDiffuse = false;
				randomDirection = randomPhong(rayDirection.reflect(N), currentObject->get_phongExponent());
				if (dot(randomDirection, N) < 0) return Vector(0., 0., 0.);
				if (dot(randomDirection, R) < 0) return Vector(0., 0., 0.);
			}

			Ray randomRay(P + 1e-4 * N, randomDirection);

			double phongProba = (currentObject->get_phongExponent() + 1) / M_PI * std::pow(dot(R, randomDirection), currentObject->get_phongExponent());
			double proba = p * dot(N, randomDirection) / M_PI + (1. - p) * phongProba;

			if (sampleDiffuse)
				I += getColor(randomRay, numRebound - 1, false) * currentObject->get_albedo() * dot(N, randomDirection) / M_PI / proba;
			else
				I += getColor(randomRay, numRebound - 1, false) * dot(N, randomDirection) * phongBRDF(randomDirection, rayDirection, N, currentObject->get_phongExponent()) * currentObject->get_ks() * currentObject->get_albedo() / proba;
		}
	}
	return I;
}

// Setters
// =======

void Scene::set_camera(Vector C)
{
	camera = C;
}

void Scene::set_fov(double f)
{
	fov = f;
}

void Scene::set_light(Sphere* L)
{
	light = L;
}

void Scene::set_lightIntensity(double intensite)
{
	lightIntensity = intensite;
}

void Scene::set_refractiveIndex(double refraction)
{
	refractiveIndex = refraction;
}


// Getters
// =======

Vector Scene::get_camera() const
{
	return camera;
}

double Scene::get_fov() const
{
	return fov;
}


//const Sphere* Scene::operator[](int i) const {
//	return spheres[i];
//}
//
//const Sphere& Scene::operator[](int i) {
//	return *spheres[i];
//}

double Scene::phongBRDF(const Vector& wi, const Vector& wo, const Vector& N, double phongExponent)
{
	Vector reflected_wo = wo.reflect(N);
	double lobe = std::pow(dot(reflected_wo, wi), phongExponent) * (phongExponent + 2) / (2. * M_PI);
	return lobe;
}