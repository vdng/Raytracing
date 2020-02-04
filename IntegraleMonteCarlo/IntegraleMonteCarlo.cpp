// IntegraleMonteCarlo.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
// S{f(x}dx = (1/N) * S_i{f(x_i)/p(x_i)}
// Exercice avec f(x) = cos^30(x)
// p : fonction de densité gaussienne

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>


double gaussDensity(double sigma, double x) {
	return 1 / (sigma * sqrt(2 * M_PI)) * exp(-(x * x) / (2 * sigma * sigma));
}

int main()
{
	std::default_random_engine engine;
	std::uniform_real_distribution<double> distrib(0, 1);

	double S = 0;
	int N = 1e7;
	double x, r1(0), r2(0);
	double sigma = 0.25;

	for (int i = 0; i < N; i++)
	{
		r1 = distrib(engine);
		r2 = distrib(engine);
		x = sqrt(-2 * log(r1)) * cos(2*M_PI * r2) * sigma;
		S += pow(cos(x), 30) / gaussDensity(sigma, x);
	}
	S /= N;
	std::cout << S << std::endl;
	return 0;
}
