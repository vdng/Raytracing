#include "Ray.h"

Ray::Ray(const Vector& origin, const Vector& direction) : 
	origin(origin), 
	direction(direction) 
{}
Vector Ray::get_origin() const
{
	return origin;
}
Vector Ray::get_direction() const
{
	return direction;
}
;
