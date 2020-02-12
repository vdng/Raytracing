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

Ray generateRay(Vector cameraPosition, double fov, int W, int H, int i, int j, Vector direction, Vector up, double focusDistance, double aperture) {
	Vector right = cross(direction, up);

	double r1 = distrib(engine), r2 = distrib(engine);
	double R = sqrt(-2 * log(r1));
	double dx = R * cos(2 * M_PI * r2) * 0.5;
	double dy = R * sin(2 * M_PI * r2) * 0.5;

	double dx_aperture = (distrib(engine) - 0.5) * aperture;
	double dy_aperture = (distrib(engine) - 0.5) * aperture;

	Vector u((j - W / 2. + dx - 0.5) * right + (H / 2. - i + dy - 0.5) * up + (H / (2 * tan(fov / 2)) * direction)); u.normalize();
	
	Vector destination = cameraPosition + focusDistance * u;
	Vector origin = cameraPosition + Vector(dx_aperture, dy_aperture, 0.);
	
	Vector v = destination - origin; v.normalize();

	Ray r(origin, v);
	return r;
}

int main() {
	int W = 512;
	int H = 512;
	int nRays = 10;
	int numRebound = 5;
	double focusDistance = 45.;
	double aperture = 0.5;

	Sphere s0(Vector(10., 10., 20.), 7., Vector(1., 1., 1.), SphereType::normal, 0.4, 10);
	Sphere s1(Vector(0., 5., 10.), 7., Vector(0.05, 0.05, 0.05), SphereType::normal, 0.4, 1000);
	Sphere s2(Vector(-10., 0., 0.), 7., Vector(0.05, 1., 1.), SphereType::transparent, 0.4, 10);

	Sphere slum(Vector(-10, 20, 30), 10., Vector(1., 1., 1.), SphereType::light);

	Sphere splafond(Vector(0., 1000., 0.), 940, Vector(0.5, 0.5, 0.), SphereType::normal, 0.3);
	Sphere smurfond(Vector(0., 0., -1000.), 940, Vector(0.1, 0., 0.5), SphereType::normal, 0.3);
	Sphere ssol(Vector(0., -1000., 0.), 990, Vector(0.1, 0.1, 0.3), SphereType::normal, 0.3);
	Sphere smur1(Vector(-1000., 0., 0.), 940, Vector(0., 1., 0.), SphereType::normal, 0.3);
	Sphere smur2(Vector(1000., 0., 0.), 940, Vector(0., 0., 1.), SphereType::normal, 0.3);

	double totalIntensity = 1e9;
	Scene scene(slum, totalIntensity);
	scene.addSphere(slum);
	scene.addSphere(splafond);
	scene.addSphere(smurfond);
	scene.addSphere(ssol);
	scene.addSphere(smur1);
	scene.addSphere(smur2);

	scene.addSphere(s0);
	scene.addSphere(s1);
	scene.addSphere(s2);

	scene.set_camera(Vector(0., 0., 55.));
	scene.set_fov(M_PI / 3.);

	Vector direction(0., 0, -1.); /*direction.normalize();*/
	Vector up(0., 1., 0.); /*direction.normalize();*/

	scene.set_refractiveIndex(1.);

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector color(0., 0., 0.);

			for (int k = 0; k < nRays; k++)
			{
				Ray r = generateRay(scene.get_camera(), scene.get_fov(), W, H, i, j, direction, up, focusDistance, aperture);
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