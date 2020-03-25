#pragma once
#include "Vector.h"
#include "Ray.h"

class BBox
{
public:
	BBox();
	BBox(const Vector& bmin, const Vector& bmax);

	bool intersect(const Ray& r) const;

	Vector bmin, bmax;
};

