#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "thermal_radiation.h"

int negCompare(const void *a, const void *b) {
	if ((*(double*)b - *(double*)a) < 0)
		return -1;
	else if ((*(double*)b - *(double*)a) > 0)
		return 1;
	else
		return 0;
}

double TToTheta(double T, double p) {
	return T * pow(1000./p, R_A/C_P);
}

double ThetaToT(double Theta, double p) {
	return Theta * pow(p/1000., R_A/C_P);
}

double barometric_PToZ(double p, double T_b, double p_b) {
	return -1 * R_STAR * T_b / (M_AIR * G) * log(p / p_b);
}

void convection(double *temperature, double *pressure_layers, int nlayers) {
	double pot_temperature[nlayers];

	for (int i=0; i < nlayers; i++) {
		pot_temperature[i] = TToTheta(temperature[i], pressure_layers[i]);
	}
	qsort(&pot_temperature, nlayers, sizeof(double), negCompare);
	for (int i=0; i < nlayers; i++) {
		temperature[i] = ThetaToT(pot_temperature[i], pressure_layers[i]);
		// DEBUG: printf("layer %2d :: temperature %8.2fK :: potential temperature %8.2fK\n", i, temperature[i], pot_temperature[i]);
	}
}

void band_flux(int nlevels, int nbands, double albedo, double *lambda_bands, double *dtau_windows, double *temperature_layers, double *total_Edn_levels, double *total_Eup_levels) {
	int nlayers = nlevels - 1;
	double Eup_levels[nlevels];
	double Edn_levels[nlevels];
	double B_layers[nlayers];  /* unit: W * s / (m**2 * r) */

	/* Reset energy fluxes */
	for (int i=0; i < nlevels; i++) {
		total_Eup_levels[i] = 0;
		total_Edn_levels[i] = 0;
	}

	/* Atmospheric window: varying tau; e.g. for three windows
	 * lambda < 1st val: tau varying between given bounds (vertically integrated)
	 * 1st val < lambda < 2nd val: tau varying between given bounds (vertically integrated)
	 * 2nd val < lambda: tau varying between given bounds (vertically integrated)
	 */
	double dtau_layers[nlayers];
	for (int n=0; n < nbands - 1; n++) {
		for (int i=0; i < nlayers; i++) {
			B_layers[i] = planck_int(lambda_bands[n], lambda_bands[n+1], temperature_layers[i]);
			dtau_layers[i] = dtau_windows[n] / (double) nlayers;
		}

		double B_surface = planck_int(lambda_bands[n], lambda_bands[n+1], temperature_layers[nlayers-1]);
		schwarzschild(nlevels, albedo, dtau_layers, B_layers, B_surface, Edn_levels, Eup_levels);
		for (int i=0; i < nlevels; i++) {
			total_Eup_levels[i] += Eup_levels[i];
			total_Edn_levels[i] += Edn_levels[i];
		}
	}
}

void heating(double *temperature, double delta_t, double p0, int nlayers) {
	double delta_p = p0/nlayers;
	temperature[nlayers-1] += E_ABS * delta_t * G / (delta_p * C_P);
}

void cooling(double *temperature, double delta_t, double p0, int nlayers) {
	double delta_p = p0/nlayers;
	double emissivity = 1./3.;
	temperature[nlayers-1] -= emissivity * SIGMA_SB * pow(temperature[nlayers-1], 4) * delta_t * G / (delta_p * C_P);
}

int main() {
	int niterations = 100000;  /* Number of iterations to run the model */
	int nlayers = 20;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */

	double p0 = 1000;  /* unit: hPa */
	double albedo = 0.;

	double pressure_layers[nlayers];
	double temperature_layers[nlayers];  /* unit: Kelvin */
	double z_layers[nlayers];

	double pressure_levels[nlevels];
	double temperature_levels[nlevels];  /* unit: Kelvin */
	double z_levels[nlevels];  /* altitude; unit: m */
	double total_Eup_levels[nlevels];
	double total_Edn_levels[nlevels];

	printf("Initializing arrays...\n");
	for (int i=0; i < nlevels; i++) {
		pressure_levels[i] = p0/nlayers * i;
		printf("level %2d :: pressure %8.2fhPa\n", i, pressure_levels[i]);
	}
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = ( pressure_levels[i] + pressure_levels[i+1] ) / 2;
		temperature_layers[i] = 200. + 20. * (double) i;
	}

	printf("Beginning iterative climate modelling...\n");
	int it = 0;
	while (it < niterations) {
		it++;

		/* Not-yet adaptive time difference between heating steps */
		double delta_t = 3600;
		double delta_p = p0/nlayers;

		//heating(temperature_layers, delta_t, p0, nlayers);
		//cooling(temperature_layers, delta_t, p0, nlayers);

		convection(temperature_layers, pressure_layers, nlayers);

		double lambda_bands[4] = {1E-20, 8E-6, 12E-6, 1E-03};
		double dtau_windows[3] = {1., 0., 1.};
		band_flux(nlevels, 4, albedo, lambda_bands, dtau_windows, temperature_layers, total_Edn_levels, total_Eup_levels);

		if (it%100 == 0) {
			printf("iteration %4d\n", it);
			for (int i=0; i < nlevels; i++) {
				printf("total_Edn_levels: %6.4f, total_Eup_levels: %6.4f, temperature_layers: %6.4f\n", total_Edn_levels[i], total_Eup_levels[i], temperature_layers[i]);
			}
		}

		/* Preemptively adapt the value of Eup at the surface to make the following loop work.
		 * This is necessary as the surface does not transfer heat any lower. Furthermore, apply surface heating.
		 */
		total_Eup_levels[nlevels-1] = total_Edn_levels[nlevels-1] + E_ABS;
		for (int i=0; i < nlayers; i++) {
			double delta_E_abs = (total_Edn_levels[i] - total_Eup_levels[i]) - (total_Edn_levels[i+1] - total_Eup_levels[i+1]);
			/* Convert between hPa and Pa by multiplying delta_p with 100 */
			temperature_layers[i] += delta_E_abs * G / (100 * delta_p * C_P) * delta_t;
		}
	}

	temperature_levels[0] = T_UNIVERSE;
	temperature_levels[nlevels-1] = temperature_layers[nlayers-1];
	for (int i=1; i < nlayers; i++) {
		temperature_levels[i] = (temperature_layers[i] + temperature_layers[i-1]) / 2.;
	}
	for (int i=0; i < nlevels; i++) {
		z_levels[i] = barometric_PToZ(pressure_levels[i], temperature_layers[nlayers-1], p0);
	}
	for (int i=0; i < nlayers; i++) {
		z_layers[i] = (z_levels[i] + z_levels[i+1]) / 2.;
	}

	printf("Modelled results by levels...\n");
	for (int i=0; i < nlevels; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm\n", i, pressure_levels[i], temperature_levels[i], z_levels[i]);
	}
	printf("Modelled results by layers...\n");
	for (int i=0; i < nlayers; i++) {
		printf("layer %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm\n", i, pressure_layers[i], temperature_layers[i], z_layers[i]);
	}

	return 0;
}
