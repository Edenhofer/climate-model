#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define R_a 287
#define c_p 1004.

int negCompare(const void * a, const void * b) {
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

int main() {
	double p0 = 1000;  /* unit: hPa */
	int nlayers = 10;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double pressure_layers[nlayers];
	double pressure[nlevels];
	double temperature[nlayers];  /* Kelvin */
	double pot_temperature[nlevels];

	printf("Initializing pressure array...\n");
	for (int i=0; i < nlevels; i++) {
		pressure[i] = p0/nlayers * i;
		printf("%f\n", pressure[i]);
	}

	printf("Initializing arrays...\n");
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = ( pressure[i] + pressure[i+1] ) / 2;
		temperature[i] = 200. + 20. * (double) i;

		pot_temperature[i] = TToTheta(temperature[i], pressure_layers[i]);
		printf("layer %2d :: temperature %8.2fK :: potential temperature %8.2fK\n", i, temperature[i], pot_temperature[i]);
	}

	qsort(&pot_temperature, nlayers, sizeof(double), negCompare);

	printf("Sorted temperature array...\n");
	for (int i=0; i < nlayers; i++) {
		temperature[i] = ThetaToT(pot_temperature[i], pressure_layers[i]);
		printf("layer %2d :: temperature %8.2fK :: potential temperature %8.2fK\n", i, temperature[i], pot_temperature[i]);
	}

	return 0;
}
