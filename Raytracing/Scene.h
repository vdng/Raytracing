#pragma once

#include <vector>
#include "Vector.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Geometry.h"

class Scene {
public:
	Scene(Sphere* light, double totalIntensity);

	void addSphere(const Sphere& sphere);
	void addTriangle(const Triangle& triangle);

	bool intersect(const Ray& r, Vector& P, Vector& N, int& idx);
	Vector getColor(const Ray& r, int numRebound, bool showLights = true);

	void set_camera(Vector C);
	void set_fov(double f);
	void set_light(Sphere* L);
	void set_lightIntensity(double intensity);
	void set_refractiveIndex(double refraction);

	Vector get_camera() const;
	double get_fov() const;

	//const Sphere* operator[](int i) const;
	//const Sphere& operator[](int i);

private:
	std::vector<const Geometry*> objects;
	Vector camera;
	double fov;
	Sphere* light;	// source lumineuse
	double lightIntensity;	// intensité de la source lumineuse
	double refractiveIndex;

	double phongBRDF(const Vector& wi, const Vector& wo, const Vector& N, double phongExponent);
};