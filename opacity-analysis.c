#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "thermal_radiation.h"

#define max(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a < _b ? _a : _b; })

void tee(FILE *f, char const *fmt, ...) {
	va_list ap;
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	va_start(ap, fmt);
	vfprintf(f, fmt, ap);
	va_end(ap);

	/* Ensure results are immediately saved */
	fflush(stdout);
	fflush(f);
}

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

int main() {
	int nlayers = 20;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double bound_delta_time = 12 * 60 * 60;  /* Bound on the timestep */
	double delta_temperature_threshold = 1e-4;  /* Temperature threshold for model termination */

	double p0 = 1000;  /* unit: hPa */
	double albedo = 0.;

	double pressure_layers[nlayers];
	double temperature_layers[nlayers];  /* unit: Kelvin */

	double pressure_levels[nlevels];
	double total_Eup_levels[nlevels];
	double total_Edn_levels[nlevels];

	for (int i=0; i < nlevels; i++) {
		pressure_levels[i] = p0/nlayers * i;
	}
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = (pressure_levels[i] + pressure_levels[i+1]) / 2;
		temperature_layers[i] = 288. - ((double) nlayers - (double) i - 1.) * 50./((double) nlayers - 1.);
	}

	/* Open file for logging and make sure it is always closed */
	FILE *fp;
	fp = fopen("opacity-analysis.csv", "w");

	tee(fp, "%8s %30s %30s\n", "Opacity", "Surface-Temperature_Grey[K]", "Surface-Temperature_Window[K]");
	int nopacities = 100;
	for (int outer=0; outer < nopacities; outer++) {
		double temperature_surface_grey;
		double temperature_surface_window;

		double opacity;
		int opacity_log_threshold = (int) ((float) nopacities / 2.);
		double opacity_log_offset = log((float) opacity_log_threshold);
		int actual_outer = outer;
		if (outer%2 == 0) {
			actual_outer = outer+1;
		}
		if (outer < opacity_log_threshold) {
			opacity = log(actual_outer);
		} else {
			opacity = opacity_log_offset + (float) (actual_outer - opacity_log_threshold);
		}

		double lambda_bands[4] = {1E-20, 8E-6, 12E-6, 1E-03};
		double dtau_windows[3];
		if (outer%2 == 0) {
			/* Gray atmosphere */
			dtau_windows[0] = opacity;
			dtau_windows[1] = opacity;
			dtau_windows[2] = opacity;

		} else if (outer%2 == 1) {
			/* Window atmosphere */
			dtau_windows[0] = opacity;
			dtau_windows[1] = 0.;
			dtau_windows[2] = opacity;
		}

		/* Single iteration of the climate model */
		int it = 0;
		double model_t = 0.;
		double bound_delta_temperature = 1.;
		double temperature_sum_prev = 0.;
		double temperature_sum_curr = temperature_sum_prev + delta_temperature_threshold + 1.;
		while (fabs(temperature_sum_curr - temperature_sum_prev) > delta_temperature_threshold) {
			it++;
			temperature_sum_prev = temperature_sum_curr;
			temperature_sum_curr = 0.;

			double delta_p = p0/nlayers;

			convection(temperature_layers, pressure_layers, nlayers);

			band_flux(nlevels, 4, albedo, lambda_bands, dtau_windows, temperature_layers, total_Edn_levels, total_Eup_levels);

			double max_E_net = 1.;
			/* Preemptively adapt the value of Eup at the surface to make the following loop work.
			 * This is necessary as the surface does not transfer heat any lower. Furthermore, apply surface heating.
			 */
			total_Eup_levels[nlevels-1] = total_Edn_levels[nlevels-1] + E_ABS;
			for (int i=0; i < nlayers; i++) {
				double delta_E_abs = (total_Edn_levels[i] - total_Eup_levels[i]) - (total_Edn_levels[i+1] - total_Eup_levels[i+1]);
				max_E_net = max(fabs(delta_E_abs), max_E_net);
			}

			/* Calculate an adaptive time difference and add it to the runtime.
			 * In the calculation convert between hPa and Pa by multiplying delta_p with 100.
			 */
			double delta_t = min(C_P/G * (100 * delta_p)/max_E_net * bound_delta_temperature, bound_delta_time);
			model_t += delta_t;
			/* Actually adapt the temperature after having defined the time step; also see previous loop */
			for (int i=0; i < nlayers; i++) {
				double delta_E_abs = (total_Edn_levels[i] - total_Eup_levels[i]) - (total_Edn_levels[i+1] - total_Eup_levels[i+1]);
				temperature_layers[i] += delta_E_abs * G / (100 * delta_p * C_P) * delta_t;
				temperature_sum_curr += temperature_layers[i];
			}
		}
		if (outer%2 == 0) {
			temperature_surface_grey = temperature_layers[nlayers-1];
		} else if (outer%2 == 1) {
			temperature_surface_window = temperature_layers[nlayers-1];
		}

		if (outer%2 == 1) {
			tee(fp, "%8.2f %30.4f %30.4f\n", (float) opacity, temperature_surface_grey, temperature_surface_window);
		}
	}

	/* Close file */
	fclose(fp);

	return 0;
}
