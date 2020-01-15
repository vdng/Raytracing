#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"

#include <algorithm>    // std::max

int main() {
	int W = 512;
	int H = 512;

	const Vector O;
	const Vector rho(1., 1., 0.05);
	Sphere s(O, 10., rho);

	Sphere splafond(Vector(0., 1000., 0.), 940, Vector(0.5, 0.5, 0.));
	Sphere smurfond(Vector(0., 0., -1000.), 940, Vector(0.1, 0., 0.5));
	Sphere ssol(Vector(0., -1000., 0.), 990, Vector(0.1, 0.1, 0.3));
	Sphere smur1(Vector(-1000., 0., 0.), 940, Vector(0., 0., 1.));
	Sphere smur2(Vector(1000., 0., 0.), 940, Vector(0., 0., 1.));

	std::vector<Sphere*> spheres = {&s, &splafond, &smurfond, &ssol};
	const Vector C(0, 0, 55);  // caméra
	double fov = M_PI / 3;
	const Vector L(-10, 20, 40);	// source lumineuse
	const double intensiteL = 3 * 10e7;	// intensité de la source lumineuse
	Scene scene(spheres, C, fov, L, intensiteL);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W / 2., H / 2. - i, - H / (2 * tan(fov / 2)));
			u.normalize();
			Ray r(C, u);

			Vector P, N;
			int idx;

			bool has_intersect = scene.intersect(r, P, N, idx);

			if (has_intersect) 
			{
				Vector I = scene.getColor(r, P, N, idx);
				image[(i * W + j) * 3 + 0] = std::min(pow(I[0], 0.45), 255.);
				image[(i * W + j) * 3 + 1] = std::min(pow(I[1], 0.45), 255.);
				image[(i * W + j) * 3 + 2] = std::min(pow(I[2], 0.45), 255.);
			}
			else 
			{
				image[(i * W + j) * 3 + 0] = 0;
				image[(i * W + j) * 3 + 1] = 0;
				image[(i * W + j) * 3 + 2] = 0;
			}

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}