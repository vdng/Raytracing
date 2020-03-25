#include "BBox.h"
#include <algorithm>

BBox::BBox()
{
}

BBox::BBox(const Vector& bmin, const Vector& bmax):
	bmin(bmin), bmax(bmax)
{}

bool BBox::intersect(const Ray & r) const
{
	Vector rayDirection = r.get_direction();
	Vector rayOrigin = r.get_origin();

	double t1x = (bmin[0] - rayOrigin[0]) / rayDirection[0];
	double t2x = (bmax[0] - rayOrigin[0]) / rayDirection[0];
	double tminx = std::min(t1x, t2x);
	double tmaxx = std::max(t1x, t2x);

	double t1y = (bmin[1] - rayOrigin[1]) / rayDirection[1];
	double t2y = (bmax[1] - rayOrigin[1]) / rayDirection[1];
	double tminy = std::min(t1y, t2y);
	double tmaxy = std::max(t1y, t2y);

	double t1z = (bmin[2] - rayOrigin[2]) / rayDirection[2];
	double t2z = (bmax[2] - rayOrigin[2]) / rayDirection[2];
	double tminz = std::min(t1z, t2z);
	double tmaxz = std::max(t1z, t2z);

	if (std::min(std::min(tmaxx, tmaxy), tmaxz) - std::max(std::min(tminx, tminy), tminz) > 0)
		return true;
	return false;
}

