#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <random>

#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"

#include <algorithm>    // std::max

extern std::default_random_engine engine;
extern std::uniform_real_distribution<double> distrib;

Ray generateRay(Vector camera, double fov, int W, int H, int i, int j) {
	double r1 = distrib(engine), r2 = distrib(engine);
	double R = sqrt(-2 * log(r1));
	double x = R * cos(2 * M_PI * r2) * 0.5;
	double y = R * sin(2 * M_PI * r2) * 0.5;

	Vector u(j - W / 2. + x - 0.5, H / 2. - i + y - 0.5, -H / (2 * tan(fov / 2))); u.normalize();
	Ray r(camera, u);
	return r;
}

int main() {
	int W = 512;
	int H = 512;
	int nRays = 100;

	Sphere s1(Vector(10., 0., 0.), 10., Vector(1., 1., 0.05), false, false, 1.5);
	Sphere s2(Vector(-10., 0., 0.), 10., Vector(0.05, 1., 1.), false, true, 1.5);

	Sphere slum(Vector(-10, 20, 40), 10., Vector(1., 1., 1.), false, false);

	Sphere splafond(Vector(0., 1000., 0.), 940, Vector(0.5, 0.5, 0.));
	Sphere smurfond(Vector(0., 0., -1000.), 940, Vector(0.1, 0., 0.5));
	Sphere ssol(Vector(0., -1000., 0.), 990, Vector(0.1, 0.1, 0.3));
	Sphere smur1(Vector(-1000., 0., 0.), 940, Vector(0., 1., 0.));
	Sphere smur2(Vector(1000., 0., 0.), 940, Vector(0., 0., 1.));

	Scene scene;
	scene.addSphere(slum);
	scene.addSphere(splafond);
	scene.addSphere(smurfond);
	scene.addSphere(ssol);
	scene.addSphere(smur1);
	scene.addSphere(smur2);
	scene.addSphere(s1);
	scene.addSphere(s2);

	scene.set_camera(Vector(0., 0., 55.));
	scene.set_fov(M_PI / 3.);

	scene.set_light(slum);
	scene.set_intensiteL(3* 1e9);
	scene.set_refractiveIndex(1.);

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			//Vector u(j - W / 2., H / 2. - i, - H / (2 * tan(scene.get_fov() / 2)));
			//u.normalize();
			//Ray r(scene.get_camera(), u);

			int numRebound = 5;
			Vector color(0., 0., 0.);
			for (int k = 0; k < nRays; k++)
			{
				Ray r = generateRay(scene.get_camera(), scene.get_fov(), W, H, i, j);
				color += scene.getColor(r, numRebound) / nRays;
			}

			image[(i * W + j) * 3 + 0] = std::min(pow(color[0], 0.45), 255.);
			image[(i * W + j) * 3 + 1] = std::min(pow(color[1], 0.45), 255.);
			image[(i * W + j) * 3 + 2] = std::min(pow(color[2], 0.45), 255.);

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}