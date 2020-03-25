#include "TriangleIndexes.h"

TriangleIndexes::TriangleIndexes(
	int vtxi, int vtxj, int vtxk,
	int ni, int nj, int nk,
	int uvi, int uvj, int uvk) :
	vtxi(vtxi), vtxj(vtxj), vtxk(vtxk),
	uvi(uvi), uvj(uvj), uvk(uvk),
	ni(ni), nj(nj), nk(nk) {
};