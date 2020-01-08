#pragma once

#include <vector>
#include "Vector.h"
#include "Sphere.h"

class Scene {
	std::vector<Sphere*> spheres;

public:
	Scene(std::vector<Sphere*> spheres);

	int intersect(const Ray& r, Vector& P, Vector& N);

	Sphere operator[](int i) const;

	Sphere& operator[](int i);
};
