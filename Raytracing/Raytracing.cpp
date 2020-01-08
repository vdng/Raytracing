#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <algorithm>    // std::max

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) {
		coords[0] = x;
		coords[1] = y;
		coords[2] = z;
	}

	double operator[](int i) const {
		return coords[i];
	}

	double& operator[](int i) {
		return coords[i];
	}

	double getNorm2() const {
		return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
	}

	void normalize() {
		double n = sqrt(getNorm2());
		coords[0] /= n;
		coords[1] /= n;
		coords[2] /= n;
	}

	double coords[3];
};

Vector operator+(const Vector& A, const Vector& B) {
	return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
};

Vector operator-(const Vector& A, const Vector& B) {
	return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}

Vector operator*(double a, const Vector& B) {
	return Vector(a * B[0], a * B[1], a * B[2]);
}

Vector operator*(const Vector& A, double b) {
	return Vector(A[0] * b, A[1] * b, A[2] * b);
}

double dot(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

class Ray {
public:
	Ray(const Vector& C, Vector u) : C(C), u(u) {};
	Vector C, u;
};

class Sphere {
public:
	Sphere(const Vector& O, double R) : O(0), R(R) {};

	double intersect(const Ray& r) {
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).getNorm2() - R * R;

		double delta = b * b - 4 * a * c;
		if (delta < 0) return -1;

		double t1 = (-b - sqrt(delta)) / (2 * a);
		double t2 = (-b + sqrt(delta)) / (2 * a);

		if (t1 >= 0) return t1;
		else return t2;
	}

	Vector O;
	double R;
};


int main() {
	int W = 512;
	int H = 512;

	Vector O;
	Sphere s(O, 10);
	Vector C(0, 0, 55);

	double fov = M_PI / 3;

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W / 2, i - W / 2, - W / (2 * tan(fov / 2)));
			u = u - C;
			u.normalize();
			Ray r(C, u);

			if (s.intersect(r) >= 0) {
				image[(i * W + j) * 3 + 0] = 255;
				image[(i * W + j) * 3 + 1] = 255;
				image[(i * W + j) * 3 + 2] = 255;
			}
			else {
				image[(i * W + j) * 3 + 0] = 0;
				image[(i * W + j) * 3 + 1] = 0;
				image[(i * W + j) * 3 + 2] = 0;
			}

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}