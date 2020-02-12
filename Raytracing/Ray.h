#pragma once

#include "Vector.h"

class Ray {
public:
	Ray(const Vector& origin, const Vector& direction);

	Vector get_origin() const;
	Vector get_direction() const;

private:
	Vector origin;
	Vector direction;
};