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
	const Vector rho(0.2, 0.9, 0.7);
	Sphere s(O, 10, rho);

	Sphere mur1(Vector(0, 1000, 0), 940, Vector(00.5, 0.5, 0));
	Sphere mur2(Vector(0, 0, -1000), 940, Vector(0.1, 0, 0.5));
	Sphere mur3(Vector(0, -1000, 0), 990, Vector(0.7, 0.1, 0.3));
	Sphere mur4(Vector(0, 0, 1000), 940, Vector(1, 0, 1));

	std::vector<Sphere*> spheres = {&s, &mur1, &mur2, &mur3, &mur4};

	Scene scene(spheres);
	
	const Vector C(0, 0, 55);  // caméra

	const Vector L(-10, 20, 40);	// source lumineuse
	const double intensiteL = 3 * 10e7;	// intensité de la source lumineuse

	double fov = M_PI / 3;

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector u(j - W / 2., H / 2. - i, - H / (2 * tan(fov / 2)));
			u.normalize();
			Ray r(C, u);

			Vector P, N;

			int k = scene.intersect(r, P, N);

			if (k >= 0) 
			{
				Vector I = scene[k].intensity(r, P, N, L, intensiteL);
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