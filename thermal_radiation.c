#include <stdio.h>
#include <math.h>

#include "constants.h"

double planck(double lambda, double T) {
	/* Returns B as an irradiance! */
	return 2 * PI * H * pow(C, 2) / pow(lambda, 5) * 1 / (exp(H * C / (lambda * K_B * T)) - 1);
}

double planck_int(double lambda_lower, double lambda_upper, double T) {
	/* Returns B as an irradiance! */
	if (lambda_lower == lambda_upper) {
		return 0;
	} else if (lambda_lower > lambda_upper) {
		double swap = lambda_lower;
		lambda_lower = lambda_upper;
		lambda_upper = swap;
	}

	/* Go one step per every half-micrometer */
	int nsteps = (lambda_upper - lambda_lower) / .5E-6 + 1;

	double delta_lambda = (lambda_upper - lambda_lower) / (double) nsteps;
	double sum = 0;
	for (int i=0; i < nsteps; i++) {
		double lambda = lambda_lower + delta_lambda * (i + 0.5);
		sum += planck(lambda, T) * delta_lambda;
	}

	return sum;
}

void schwarzschild(int nlevels, double albedo, double *dtau_layers, double *B_layers, double B_surface, double *Edn_levels, double *Eup_levels) {
	/* Note, B should be an irradiance! */
	double L_0 = 0.;
	int nsteps = 20;

	/* Initialize energy arrays */
	for (int i=0; i < nlevels; i++) {
		Eup_levels[i] = 0;
		Edn_levels[i] = 0;
	}

	double delta_mu = 1. / (double) nsteps;
	for (int n=0; n < nsteps; n++) {
		double mu = delta_mu * (n + 0.5);

		double L = L_0;
		/* Down */
		Edn_levels[0] += 2. * PI * L * mu * delta_mu;
		for (int i=0; i < nlevels-1; i++) {
			L = L * exp(-dtau_layers[i] / mu) + (1 - exp(-dtau_layers[i] / mu)) * B_layers[i] / PI;
			// DEBUG: printf("Ldn %8.f; B %8.2f; B_surface %8.2f, dtau %8.2f\n", L, B_layers[i] / pi, B_surface, dtau_layers[i]);
			Edn_levels[i+1] += 2. * PI * L * mu * delta_mu;
		}
		/* Up */
		L = B_surface / PI * (1. - albedo) + L * albedo;
		Eup_levels[nlevels-1] += 2. * PI * L * mu * delta_mu;
		for (int i=nlevels-2; i >= 0; i--) {
			L = L * exp(-dtau_layers[i] / mu) + (1 - exp(-dtau_layers[i] / mu)) * B_layers[i] / PI;
			// DEBUG: printf("Lup %8.f; B %8.2f; B_surface %8.2f, dtau %8.2f\n", L, B_layers[i] / pi, B_surface, dtau_layers[i]);
			Eup_levels[i] += 2. * PI * L * mu * delta_mu;
		}
	}

}
