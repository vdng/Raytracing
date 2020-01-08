#include "Scene.h"

Scene::Scene(std::vector<Sphere*> spheres) : spheres(spheres) {};

int Scene::intersect(const Ray& r, Vector& P, Vector& N)
{
	int idx = -1;	// indice de la sphere la plus proche ; -1 si aucune intersection
	Vector Plocal, Nlocal;
	double t = std::numeric_limits<double>::max();
	for (size_t i = 0; i < spheres.size(); i++)
	{
		double tlocal = spheres[i]->intersect(r, Plocal, Nlocal);
		if (0 <= tlocal && tlocal < t)
		{
			idx = i;
			t = tlocal;
			P = Plocal;
			N = Nlocal;
		}
	}
	return idx;
}

Sphere Scene::operator[](int i) const {
	return *spheres[i];
}

Sphere& Scene::operator[](int i) {
	return *spheres[i];
}