#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"

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
		printf("layer %2d :: temperature %8.2fK :: potential temperature %8.2fK\n", i, temperature[i], pot_temperature[i]);
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
	int niterations = 20;  /* Number of iterations to run the model */
	int nlayers = 10;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double delta_t = 5;  /* time difference between heating steps */

	double p0 = 1000;  /* unit: hPa */
	double pressure_levels[nlevels];
	double temperature_levels[nlevels];  /* unit: Kelvin */
	double z_levels[nlevels];  /* altitude; unit: m */
	double pressure_layers[nlayers];
	double temperature_layers[nlayers];  /* unit: Kelvin */

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

		printf("iteration %4d\n", it);
		heating(temperature_layers, delta_t, p0, nlayers);
		cooling(temperature_layers, delta_t, p0, nlayers);
		convection(temperature_layers, pressure_layers, nlayers);

		for (int i=0; i < nlevels; i++) {
			z_levels[i] = barometric_PToZ(pressure_levels[i], temperature_layers[nlayers-1], p0);
		}
	}

	temperature_levels[0] = T_UNIVERSE;
	temperature_levels[nlevels-1] = temperature_layers[nlayers-1];
	for (int i=1; i < nlayers; i++) {
		temperature_levels[i] = (temperature_layers[i] + temperature_layers[i-1]) / 2.;
	}

	printf("Print resulting arrays...\n");
	for (int i=0; i < nlevels; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm\n", i, pressure_levels[i], temperature_levels[i], z_levels[i]);
	}

	return 0;
}
