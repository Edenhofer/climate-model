#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define R_a 287
#define c_p 1004.

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
		printf("layer %2d :: temperature %8.2fK :: potential temperature %8.2fK\n", i, temperature[i], pot_temperature[i]);
	}
	qsort(&pot_temperature, nlayers, sizeof(double), negCompare);
	for (int i=0; i < nlayers; i++) {
		temperature[i] = ThetaToT(pot_temperature[i], pressure_layers[i]);
	}
}

void heating(double *temperature) {

}

int main() {
	double p0 = 1000;  /* unit: hPa */
	int nlayers = 10;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
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

	convection(temperature, pressure_layers, nlayers);

	printf("Print sorted temperature array...\n");
	for (int i=0; i < nlayers; i++) {
		printf("layer %2d :: temperature %8.2fK\n", i, temperature[i]);
	}

	return 0;
}
