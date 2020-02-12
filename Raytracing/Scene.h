#pragma once

#include <vector>
#include "Vector.h"
#include "Sphere.h"

class Scene {
private:
	std::vector<Sphere> spheres;
	Vector camera;
	double fov;
	Sphere light;	// source lumineuse
	double intensiteL;	// intensité de la source lumineuse
	double refractiveIndex;

	double phongBRDF(const Vector& wi, const Vector& wo, const Vector& N, double phongExponent);

public:
	Scene(Sphere light, double totalIntensity);

	bool intersect(const Ray& r, Vector& P, Vector& N, int& idx);
	Vector getColor(const Ray& r, int numRebound, bool showLights = true);

	void set_camera(Vector C);
	void set_fov(double f);
	void set_light(Sphere L);
	void set_intensiteL(double intensite);
	void set_refractiveIndex(double refraction);

	Vector get_camera();
	double get_fov();

	void addSphere(const Sphere& sphere);

	Sphere operator[](int i) const;
	Sphere& operator[](int i);
};