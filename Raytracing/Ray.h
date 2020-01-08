#pragma once

#include "Vector.h"

class Ray {
public:
	Vector C, u;
	
	Ray(const Vector& C, Vector u);
};