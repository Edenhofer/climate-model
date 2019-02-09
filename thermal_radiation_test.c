#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "radiative_transfer.h"

int main() {
	int nlayers = 10;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */

	/* Hardcode the initialization values for 10 layers */
	double p0 = 1000;  /* unit: hPa */
	double albedo = 0.;

	double pressure_layers[nlayers];
	double temperature_layers[nlayers];  /* unit: Kelvin */
	double z_layers[nlayers];
	double dtau_layers[nlayers];
	double B_layers[nlayers];  /* unit: W * s / (m**2 * r) */

	double pressure_levels[11] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};  /* unit: hPa */
	double temperature_levels[11] = {127.28, 187.09, 212.42, 229.22, 242.03, 252.48, 261.37, 269.13, 276.04, 282.29, 288.00};  /* unit: Kelvin */
	double z_levels[nlevels];  /* altitude; unit: m */
	double Eup_levels[nlevels];
	double Edn_levels[nlevels];

	printf("Initializing arrays...\n");
	for (int i=0; i < nlevels; i++) {
		z_levels[i] = -1 * R_STAR * temperature_levels[nlayers-1] / (M_AIR * G) * log(pressure_levels[i] / p0);
	}
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = (pressure_levels[i] + pressure_levels[i+1]) / 2.;
		temperature_layers[i] = (temperature_levels[i] + temperature_levels[i+1]) / 2.;
		z_layers[i] = (z_levels[i] + z_levels[i+1]) / 2.;
		B_layers[i] = SIGMA_SB * pow(temperature_layers[i], 4);
	}
	double B_surface = SIGMA_SB * pow(temperature_levels[nlevels-1], 4);

	printf("Print levels...\n");
	for (int i=0; i < nlevels; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm\n", i, pressure_levels[i], temperature_levels[i], z_levels[i]);
	}
	printf("Print layers...\n");
	for (int i=0; i < nlayers; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm\n", i, pressure_layers[i], temperature_layers[i], z_layers[i]);
	}

	printf("Grey atmosphere: delta tau equals to 1\n");
	for (int i=0; i < nlayers; i++) {
		dtau_layers[i] = 1. / (double) nlayers;
	}
	schwarzschild(nlevels, albedo, dtau_layers, B_layers, B_surface, Edn_levels, Eup_levels);
	for (int i=0; i < nlevels; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm :: E_dn %8.3f :: E_up %8.3f\n", i, pressure_levels[i], temperature_levels[i], z_levels[i], Edn_levels[i], Eup_levels[i]);
	}

	printf("Grey atmosphere: delta tau equals to 10\n");
	for (int i=0; i < nlayers; i++) {
		dtau_layers[i] = 10. / (double) nlayers;
	}
	schwarzschild(nlevels, albedo, dtau_layers, B_layers, B_surface, Edn_levels, Eup_levels);
	for (int i=0; i < nlevels; i++) {
		printf("level %2d :: pressure %8.2f :: temperature %8.2fK :: altitude %8.2fm :: E_dn %8.3f :: E_up %8.3f\n", i, pressure_levels[i], temperature_levels[i], z_levels[i], Edn_levels[i], Eup_levels[i]);
	}

	printf("Atmospheric window: varying tau\n");
	/* lambda < 8 micron: tau varying between 0 and 10 (vertically integrated)
	 * 8 micron < lambda < 12 micron: tau=0 (vertically integrated)
	 * 12 micron < lambda: tau varying between 0 and 10 (vertically integrated)
	 */
	double lambda_bands[4] = {1E-20, 8E-6, 12E-6, 1E-01};
	for (int k=0; k < 50; k++) {
		double total_dtau_layers = 0.1 * k;
		double dtau_windows[3] = {total_dtau_layers, 0., total_dtau_layers};

		double total_Edn_surface = 0.;
		double total_Eup_TOA = 0.;

		for (int n=0; n < 3; n++) {
			for (int i=0; i < nlayers; i++) {
				B_layers[i] = planck_int(lambda_bands[n], lambda_bands[n+1], temperature_layers[i]);
				dtau_layers[i] = dtau_windows[n] / (double) nlayers;
			}

			double B_surface = planck_int(lambda_bands[n], lambda_bands[n+1], temperature_levels[nlevels-1]);
			schwarzschild(nlevels, albedo, dtau_layers, B_layers, B_surface, Edn_levels, Eup_levels);
			total_Edn_surface += Edn_levels[nlevels-1];
			total_Eup_TOA += Eup_levels[0];
		}

		printf("%2.1f :: %7.3f :: %7.3f\n", total_dtau_layers, total_Edn_surface, total_Eup_TOA);
	}

	return 0;
}
