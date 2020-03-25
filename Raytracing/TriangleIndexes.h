#pragma once
class TriangleIndexes {
public:
	TriangleIndexes(
		int vtxi = -1, int vtxj = -1, int vtxk = -1,
		int ni = -1, int nj = -1, int nk = -1,
		int uvi = -1, int uvj = -1, int uvk = -1);

	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int faceGroup;
};
