#pragma once

#include <vector>
#include "Vector.h"
#include "Object.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Geometry.h"

class Scene {
public:
	Scene(Sphere* light, double totalIntensity);

	// Methods
	void addSphere(const Sphere& sphere);
	void addTriangle(const Triangle& triangle);
	void addGeometry(const Geometry& geometry);

	bool intersect(const Ray& r, Vector& P, Vector& N, int& idx);
	Vector getColor(const Ray& r, int numRebound, bool showLights = true);

	// Setters
	void set_camera(Vector C);
	void set_fov(double f);
	void set_light(Sphere* L);
	void set_lightIntensity(double intensity);
	void set_refractiveIndex(double refraction);

	// Getters
	Vector get_camera() const;
	double get_fov() const;

	//const Sphere* operator[](int i) const;
	//const Sphere& operator[](int i);

private:
	std::vector<const Object*> objects;
	Vector camera;
	double fov;
	Sphere* light;	// source lumineuse
	double lightIntensity;	// intensité de la source lumineuse
	double refractiveIndex;

	// Methods
	double phongBRDF(const Vector& wi, const Vector& wo, const Vector& N, double phongExponent);
};