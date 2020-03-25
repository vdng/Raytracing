#pragma once
#include "Ray.h"
#include "Vector.h"
#include "Object.h"
#include "TriangleIndexes.h"
#include <vector>
#include "BBox.h"

class Geometry:
	public Object
{
public:
	Geometry(
		const char* obj, double scaling, const Vector& offset,
		const Vector& albedo = Vector(0.5, 0.5, 0.5), Material material = Material::normal,
		double ks = 0, double phongExponent = 1000.,
		double refractiveIndex = 1.5);

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const;

private:
	BBox bb;

	std::vector<TriangleIndexes> indexes;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
	std::vector<Vector> vertexcolors;

	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> w, h;

	// Methods
	void readObj(const char* obj);
};

