#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define R_a 287
#define c_p 1004.  /* unit: hPa */
#define g 9.82  /* m/s**2 */
#define E_abs 235  /* unit: W/m**2 */

int negCompare(const void *a, const void *b) {
	if ((*(double*)b - *(double*)a) < 0)
		return -1;
	else if ((*(double*)b - *(double*)a) > 0)
		return 1;
	else
		return 0;
}

double TToTheta(double T, double p) {
	return T * pow(1000./p, R_a/c_p);
}

double ThetaToT(double Theta, double p) {
	return Theta * pow(p/1000., R_a/c_p);
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
	temperature[nlayers-1] += E_abs * delta_t * g / (delta_p * c_p);
}

int main() {
	int niterations = 20;  /* Number of iterations to run the model */
	int nlayers = 10;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double delta_t = 5 * 60;  /* time difference between heating steps */

	double p0 = 1000;  /* unit: hPa */
	double pressure_layers[nlayers];
	double pressure[nlevels];
	double temperature[nlayers];  /* Kelvin */

	printf("Initializing arrays...\n");
	for (int i=0; i < nlevels; i++) {
		pressure[i] = p0/nlayers * i;
		printf("level %2d :: pressure %8.2fhPa\n", i, pressure[i]);
	}
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = ( pressure[i] + pressure[i+1] ) / 2;
		temperature[i] = 200. + 20. * (double) i;
	}

	printf("Beginning iterative climate modelling...\n");
	int it = 0;
	while (it < niterations) {
		it++;

		printf("iteration %4d\n", it);
		heating(temperature, delta_t, p0, nlayers);
		convection(temperature, pressure_layers, nlayers);
	}

	printf("Print temperature array...\n");
	for (int i=0; i < nlayers; i++) {
		printf("layer %2d :: temperature %8.2fK\n", i, temperature[i]);
	}

	return 0;
}
