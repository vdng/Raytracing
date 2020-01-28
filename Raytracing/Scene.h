#pragma once

#include <vector>
#include "Vector.h"
#include "Sphere.h"

class Scene {
	std::vector<Sphere*> spheres;
	Vector camera;
	double fov;
	Vector light;	// source lumineuse
	double intensiteL;	// intensité de la source lumineuse
	double refractiveIndex;

public:
	Scene(std::vector<Sphere*> spheres, Vector Camera, double fov, Vector Light, double intensiteL, double refractiveIndex = 1.);

	bool intersect(const Ray& r, Vector& P, Vector& N, int& idx);
	Vector getColor(const Ray& r, int numRebound);


	Sphere operator[](int i) const;

	Sphere& operator[](int i);
};
