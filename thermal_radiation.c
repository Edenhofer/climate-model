#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "eddington_v2.h"

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

void doubling_adding(int nlevels, double albedo_ground, double *r_layers, double *t_layers, double *t_direct_layers, double *s_direct_layers, double *r_direct_layers, double E_direct_toa, double *E_direct_levels, double *Edn_levels, double *Eup_levels) {
	int nlayers = nlevels - 1;

	double reflectivity_layers[nlayers];
	double transmissivity_layers[nlayers];
	double scattering_direct_layers[nlayers];
	double transmissivity_direct_layers[nlayers];

	reflectivity_layers[0] = r_layers[0];
	transmissivity_layers[0] = t_layers[0];
	transmissivity_direct_layers[0] = t_direct_layers[0];
	scattering_direct_layers[0] = s_direct_layers[0];
	E_direct_levels[0] = E_direct_toa;

	for (int i=1; i < nlayers; i++) {
		E_direct_levels[i] = transmissivity_direct_layers[i-1] * E_direct_toa;
		reflectivity_layers[i] = r_layers[i] + (reflectivity_layers[i-1] * t_layers[i] * t_layers[i]) / (1. - reflectivity_layers[i-1] * r_layers[i]);
		transmissivity_layers[i] = (transmissivity_layers[i-1] * t_layers[i]) / (1. - reflectivity_layers[i-1] * r_layers[i]);
		transmissivity_direct_layers[i] = transmissivity_direct_layers[i-1] * t_direct_layers[i];
		scattering_direct_layers[i] = (t_layers[i] * scattering_direct_layers[i-1] + transmissivity_direct_layers[i-1] * r_direct_layers[i] * reflectivity_layers[i-1] * t_layers[i]) / (1. - reflectivity_layers[i-1] * r_layers[i]) + transmissivity_direct_layers[i-1] * s_direct_layers[i];
	}
	E_direct_levels[nlevels-1] = transmissivity_direct_layers[nlayers-1] * E_direct_toa;

	Edn_levels[nlevels-1] = (scattering_direct_layers[nlayers-1] + transmissivity_direct_layers[nlayers-1] * reflectivity_layers[nlayers-1] * albedo_ground) / (1. - reflectivity_layers[nlayers-1] * albedo_ground) * E_direct_toa;
	Eup_levels[nlevels-1] = albedo_ground * (Edn_levels[nlevels-1] + transmissivity_direct_layers[nlayers-1] * E_direct_toa);
	for (int i=nlevels-2; i > 0; i--) {
		Edn_levels[i] = (reflectivity_layers[i-1] * t_layers[i] * Eup_levels[i+1]) / (1. - reflectivity_layers[i-1] * r_layers[i]) + (E_direct_toa * scattering_direct_layers[i-1] + E_direct_levels[i] * r_direct_layers[i] * reflectivity_layers[i-1]) / (1. - reflectivity_layers[i-1] * r_layers[i]);
		Eup_levels[i] = (t_layers[i] * Eup_levels[i+1]) / (1. - reflectivity_layers[i-1] * r_layers[i]) + (E_direct_toa * scattering_direct_layers[i-1] * r_layers[i] + E_direct_levels[i] * r_direct_layers[i]) / (1. - reflectivity_layers[i-1] * r_layers[i]);
	}
	Edn_levels[0] = 0.;
	Eup_levels[0] = t_layers[0] * Eup_levels[1] + r_direct_layers[0] * E_direct_toa;

	/* DEBUG code
	printf("%5s %12s %12s %12s\n", "Level", "E_dir", "Edn_lvl", "Eup_lvl");
	for (int i=0; i < nlevels; i++) {
		printf("%2d %12.4f %12.4f %12.4f\n", i, E_direct_levels[i], Edn_levels[i], Eup_levels[i]);
	}
	*/
}

void doubling_adding_eddington(int nlevels, double albedo_ground, double mu0, double E_direct_toa, double *g_layers, double *omega0_layers, double *dtau, double *E_direct_levels, double *Edn_levels, double *Eup_levels) {
	int nlayers = nlevels - 1;

	double r_layers[nlayers], t_layers[nlayers], t_direct_layers[nlayers], s_direct_layers[nlayers], r_direct_layers[nlayers];

	for (int i=0; i < nlayers; i++) {
		eddington_v2(dtau[i], g_layers[i], omega0_layers[i], mu0, &t_layers[i], &r_layers[i], &r_direct_layers[i], &s_direct_layers[i], &t_direct_layers[i]);
	}
	double E_direct_frac = E_direct_toa * mu0;
	doubling_adding(nlevels, albedo_ground, r_layers, t_layers, t_direct_layers, s_direct_layers, r_direct_layers, E_direct_frac, E_direct_levels, Edn_levels, Eup_levels);
}
