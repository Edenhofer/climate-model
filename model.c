#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "thermal_radiation.h"
#include "ascii.h"
#include "gnuplot_i.h"
#include "fpda_rrtm_lw.h"

#define max(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a < _b ? _a : _b; })

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

void rrtm_flux(int nlayers, double *pressure_levels, double *temperature_layers, double *h2ovmr, double *o3vmr, double *co2vmr, double *ch4vmr, double *n2ovmr, double *o2vmr, double *cfc11vmr, double *cfc12vmr, double *cfc22vmr, double *ccl4vmr, double albedo, double *total_Edn_levels, double *total_Eup_levels) {
	int nbands;
	double *band_lbound_bands;
	double *band_ubound_bands;
	double **dtau_mol_lw_bandslayers;
	double **wgt_lw_bandslayers;

	cfpda_rrtm_lw(nlayers, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &nbands, &band_lbound_bands, &band_ubound_bands, &wgt_lw_bandslayers, &dtau_mol_lw_bandslayers);

	int nlevels = nlayers+1;
	double Eup_levels[nlevels];
	double Edn_levels[nlevels];
	double B_layers[nlayers];  /* unit: W * s / (m**2 * r) */

	/* Reset energy fluxes */
	for (int i=0; i < nlevels; i++) {
		total_Eup_levels[i] = 0;
		total_Edn_levels[i] = 0;
	}

	for (int ib=0; ib < nbands; ib++) {
		for (int i=0; i < nlayers; i++) {
			B_layers[i] = wgt_lw_bandslayers[ib][i] * planck_int(1e-2/band_ubound_bands[ib], 1e-2/band_lbound_bands[ib], temperature_layers[i]);
		}

		double B_surface = B_layers[nlayers-1];
		schwarzschild(nlevels, albedo, dtau_mol_lw_bandslayers[ib], B_layers, B_surface, Edn_levels, Eup_levels);
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
	int nlayers = 20;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double bound_delta_time = 12 * 60 * 60;  /* Bound on the timestep */
	double delta_temperature_threshold = 1e-4;  /* Temperature threshold for model termination */

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
	double h2ovmr[nlayers], o3vmr[nlayers], co2vmr[nlayers], ch4vmr[nlayers], n2ovmr[nlayers], o2vmr[nlayers], cfc11vmr[nlayers], cfc12vmr[nlayers], cfc22vmr[nlayers], ccl4vmr[nlayers];

	printf("Initializing arrays...\n");
	printf("Level Pressure[hPa]\n");
	for (int i=0; i < nlevels; i++) {
		pressure_levels[i] = p0/nlayers * i;
		printf("%3d %12.4f\n", i, pressure_levels[i]);
	}

	int nlevels_fpda_file;
	double *h2oppm=NULL, *o3ppm=NULL, *discard1=NULL, *discard2=NULL, *discard3=NULL;
	char fpda_filepath[128] = "ascii/fpda.atm";
	int status;
	status = read_5c_file(fpda_filepath, &discard1, &discard2, &discard3, &h2oppm, &o3ppm, &nlevels_fpda_file);
	if (status != 0) {
		fprintf(stderr, "Error while opening file '%s'. Aborting...\n", fpda_filepath);
		return 1;
	}
	if (nlevels != nlevels_fpda_file) {
		fprintf(stderr, "Reading in H2O and O3 concentrations only works for %3d levels! Aborting...\n", nlevels);
		return 1;
	}
	for (int i=0; i < nlayers; i++) {
		h2ovmr[i] = (h2oppm[i] + h2oppm[i+1]) / 2. * 1e-6;
		o3vmr[i] = (o3ppm[i] + o3ppm[i+1]) / 2. * 1e-6;
	}

	printf("Layer Temperature[K] H2O O3 CO2 CH4 N2O O2\n");
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = (pressure_levels[i] + pressure_levels[i+1]) / 2;
		temperature_layers[i] = 288. - ((double) nlayers - (double) i - 1.) * 50./((double) nlayers - 1.);
		co2vmr[i] = 400e-6;
		ch4vmr[i] = 1.7e-6;
		n2ovmr[i] = 320e-9;
		o2vmr[i] = .209;
		cfc11vmr[i] = 0.;
		cfc12vmr[i] = 0.;
		cfc22vmr[i] = 0.;
		ccl4vmr[i] = 0.;

		printf("%3i %7g %7g %7g %7g %7g %7g %7g\n", i, temperature_layers[i], h2ovmr[i], o3vmr[i], co2vmr[i], ch4vmr[i], n2ovmr[i], o2vmr[i]);
	}

	gnuplot_ctrl *g1;
	g1 = gnuplot_init();

	printf("Beginning iterative climate modelling...\n");
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

		//heating(temperature_layers, delta_t, p0, nlayers);
		//cooling(temperature_layers, delta_t, p0, nlayers);

		convection(temperature_layers, pressure_layers, nlayers);

		/* Compute the flux using fixed bands
		 * ```
		 * double lambda_bands[4] = {1E-20, 8E-6, 12E-6, 1E-03};
		 * double dtau_windows[3] = {1.6, 0., 1.6};
		 * band_flux(nlevels, 4, albedo, lambda_bands, dtau_windows, temperature_layers, total_Edn_levels, total_Eup_levels);
		 * ```
		 */

		rrtm_flux(nlayers, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, albedo, total_Edn_levels, total_Eup_levels);

		if (it%100 == 0) {
			printf("Iteration %7d at time %6.2fd (%8.1fs)\n", it, model_t / (60*60*24), model_t);
			printf("total_Edn_levels[W] total_Eup_levels[W] temperature_layers[K]\n");
			for (int i=0; i < nlayers; i++) {
				printf("%12.4f %12.4f %12.4f\n", total_Edn_levels[i], total_Eup_levels[i], temperature_layers[i]);
			}
			printf("%12.4f %12.4f %12.4f\n", total_Edn_levels[nlevels-1], total_Eup_levels[nlevels-1], (double) NAN);

			for (int i=0; i < nlevels; i++) {
				z_levels[i] = barometric_PToZ(pressure_levels[i], temperature_layers[nlayers-1], p0);
			}
			for (int i=0; i < nlayers; i++) {
				z_layers[i] = (z_levels[i] + z_levels[i+1]) / 2.;
			}

			gnuplot_resetplot(g1);  /* Start with new plot rather than plotting into existing one */
			gnuplot_setstyle(g1, "linespoints");  /* Draw lines and points */
			gnuplot_set_xlabel(g1, "temperature [K]");  /* x-axis label */
			gnuplot_set_ylabel(g1, "altitude [m]");  /* y-axis label */

			/* Plot temperature T as function of z and label with temperature */
			gnuplot_plot_xy(g1, temperature_layers, z_layers, nlayers, "Temperature") ;
		}

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

	/* Close plot after receiving an input */
	printf("Iterations completed. Terminate by closing the standard input... ");
	char tmp;
	scanf(&tmp);
	gnuplot_close(g1);
	printf("\n");

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
	printf("Level Pressure[hPa] Temperature[K] Altitude[m]\n");
	for (int i=0; i < nlevels; i++) {
		printf("%3d %12.4f %12.4f %12.4f\n", i, pressure_levels[i], temperature_levels[i], z_levels[i]);
	}
	printf("Modelled results by layers...\n");
	printf("Layer Pressure[hPa] Temperature[K] Altitude[m]\n");
	for (int i=0; i < nlayers; i++) {
		printf("%3d %12.4f %12.4f %12.4f\n", i, pressure_layers[i], temperature_layers[i], z_layers[i]);
	}

	return 0;
}
