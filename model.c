#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include "constants.h"
#include "thermal_radiation.h"
#include "ascii.h"
#include "gnuplot_i.h"
#include "fpda_rrtm_sw.h"
#include "fpda_rrtm_lw.h"

#define ABSTOL 1e-10

#define max(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) \
	({ typeof(a) _a = (a); typeof(b) _b = (b); _a < _b ? _a : _b; })

int areSame(double a, double b) {
    return fabs(a - b) < ABSTOL;
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

void concentration_interpolation(double *x_lines, double *y_lines, int nlines, double *x_layers, double *y_layers, int nlayers) {
	/* Interpolate concentrations at given values for arbitrary layers and convert from ppm to vmr */
	gsl_interp *workspace;
	gsl_interp_accel *accel;

	workspace = gsl_interp_alloc(gsl_interp_linear, nlines);
	accel = gsl_interp_accel_alloc();

	gsl_interp_init(workspace, x_lines, y_lines, nlines);
	for (int i=0; i < nlayers; i++) {
		/* Interpolate between points and multiply with 1e-6 as to translate ppm to vmr */
		y_layers[i] = gsl_interp_eval(workspace, x_lines, y_lines, x_layers[i], accel) * 1e-6;
	}

	gsl_interp_accel_free(accel);
	gsl_interp_free(workspace);
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

void rrtm_lw_flux(int nlayers, double albedo, double *pressure_levels, double *temperature_layers, double *h2ovmr, double *o3vmr, double *co2vmr, double *ch4vmr, double *n2ovmr, double *o2vmr, double *cfc11vmr, double *cfc12vmr, double *cfc22vmr, double *ccl4vmr, double *liquid_water_path_layers, double cloud_frac, double r_cloud, double *wvl_lbound_bands, double *wvl_ubound_bands, double *q_ext_cloud_bands, double *omega0_cloud_bands, int nlevels_cloud_prop_file, double *total_Edn_levels, double *total_Eup_levels) {
	int nbands;
	int nlevels = nlayers+1;
	double *band_lbound_bands;
	double *band_ubound_bands;
	double **dtau_mol_lw_bandslayers;
	double **wgt_lw_bandslayers;

	cfpda_rrtm_lw(nlayers, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &nbands, &band_lbound_bands, &band_ubound_bands, &wgt_lw_bandslayers, &dtau_mol_lw_bandslayers);

	/* Check the output of the RRTM for consistency with the provided cloud property file */
	if (nlevels_cloud_prop_file != nbands) {
		fprintf(stderr, "Number of long wavelength bands returned by the RRTM differs from the one from the input file. Aborting...\n");
		exit(1);
	}
	for (int ib=0; ib < nbands; ib++) {
		if (areSame(band_lbound_bands[ib], wvl_lbound_bands[ib]) || areSame(band_ubound_bands[ib], wvl_ubound_bands[ib])) {
			fprintf(stderr, "Long wavelength bands returned by the RRTM differ from the ones from the input file. Aborting...\n");
			exit(1);
		}
	}

	/* Reset energy fluxes */
	for (int i=0; i < nlevels; i++) {
		total_Eup_levels[i] = 0;
		total_Edn_levels[i] = 0;
	}

	for (int ib=0; ib < nbands; ib++) {
		double dtau_layers[nlayers];
		double B_layers[nlayers];  /* unit: W * s / (m**2 * r) */
		double Eup_cloud_levels[nlevels], Edn_cloud_levels[nlevels];
		double Eup_plain_levels[nlevels], Edn_plain_levels[nlevels];

		for (int i=0; i < nlayers; i++) {
			B_layers[i] = wgt_lw_bandslayers[ib][i] * planck_int(1e-2/band_ubound_bands[ib], 1e-2/band_lbound_bands[ib], temperature_layers[i]);
		}

		for (int i=0; i < nlayers; i++) {
			/* Direct radiation is subject to extinction due to gas absorption and clouds (no scattering) */
			double dtau_ext_cloud = 3 * liquid_water_path_layers[i] * q_ext_cloud_bands[ib] / (4 * r_cloud * RHO_LIQUID);
			dtau_layers[i] = dtau_mol_lw_bandslayers[ib][i] + dtau_ext_cloud * omega0_cloud_bands[ib];
		}

		double B_surface = B_layers[nlayers-1];
		/* Simulate one atmosphere with and one without clouds */
		schwarzschild(nlevels, albedo, dtau_mol_lw_bandslayers[ib], B_layers, B_surface, Edn_plain_levels, Eup_plain_levels);
		schwarzschild(nlevels, albedo, dtau_layers, B_layers, B_surface, Edn_cloud_levels, Eup_cloud_levels);
		for (int i=0; i < nlevels; i++) {
			total_Eup_levels[i] += cloud_frac * Eup_cloud_levels[i] + (1. - cloud_frac) * Eup_plain_levels[i];
			total_Edn_levels[i] += cloud_frac * Edn_cloud_levels[i] + (1. - cloud_frac) * Edn_plain_levels[i];
		}
	}
}

void rrtm_sw_flux(int nlayers, double albedo_ground, double mu0, double *pressure_levels, double *temperature_layers, double *h2ovmr, double *o3vmr, double *co2vmr, double *ch4vmr, double *n2ovmr, double *o2vmr, double *liquid_water_path_layers, double cloud_frac, double r_cloud, double *wvl_lbound_bands, double *wvl_ubound_bands, double *q_ext_cloud_bands, double *omega0_cloud_bands, double *g_cloud_bands, int nlevels_cloud_prop_file, double *total_E_direct_levels, double *total_Edn_levels, double *total_Eup_levels) {
	int nbands;
	int nlevels = nlayers+1;
	double *band_lbound_bands;
	double *band_ubound_bands;
	double *wgt_sw_bandslayers;
	double **dtau_mol_sw_bandslayers;
	double **dtau_ray_sw_bandslayers;

	cfpda_rrtm_sw(nlayers, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &nbands, &band_lbound_bands, &band_ubound_bands, &wgt_sw_bandslayers, &dtau_mol_sw_bandslayers, &dtau_ray_sw_bandslayers);

	/* Check the output of the RRTM for consistency with the provided cloud property file */
	if (nlevels_cloud_prop_file != nbands) {
		fprintf(stderr, "Number of short wavelength bands returned by the RRTM differs from the one from the input file. Aborting...\n");
		exit(1);
	}
	for (int ib=0; ib < nbands; ib++) {
		if (areSame(band_lbound_bands[ib], wvl_lbound_bands[ib]) || areSame(band_ubound_bands[ib], wvl_ubound_bands[ib])) {
			fprintf(stderr, "Short wavelength bands returned by the RRTM differ from the ones from the input file. Aborting...\n");
			exit(1);
		}
	}

	/* Reset energy fluxes */
	for (int i=0; i < nlevels; i++) {
		total_E_direct_levels[i] = 0;
		total_Eup_levels[i] = 0;
		total_Edn_levels[i] = 0;
	}

	for (int ib=0; ib < nbands; ib++) {
		double dtau_cloud_layers[nlayers], g_cloud_layers[nlayers], omega0_cloud_layers[nlayers];
		double dtau_plain_layers[nlayers], g_plain_layers[nlayers], omega0_plain_layers[nlayers];
		double E_direct_cloud_levels[nlevels], Eup_cloud_levels[nlevels], Edn_cloud_levels[nlevels];
		double E_direct_plain_levels[nlevels], Eup_plain_levels[nlevels], Edn_plain_levels[nlevels];

		for (int i=0; i < nlayers; i++) {
			double dtau_ext_cloud = 3 * liquid_water_path_layers[i] * q_ext_cloud_bands[ib] / (4 * r_cloud * RHO_LIQUID);
			double dtau_sca_cloud = omega0_cloud_bands[ib] * dtau_ext_cloud;

			/* Delta-Scaling */
			double f = pow(g_cloud_bands[ib], 2);
			double g_cloud_scaled = (g_cloud_bands[ib] - f) / (1 - f);
			dtau_sca_cloud *= (1. - f);
			dtau_ext_cloud = dtau_sca_cloud + (1 - omega0_cloud_bands[ib]) * dtau_ext_cloud;

			/* Direct radiation is subject to extinction due to gas absorption, rayleigh scattering and clouds */
			dtau_cloud_layers[i] = dtau_mol_sw_bandslayers[ib][i] + dtau_ray_sw_bandslayers[ib][i] + dtau_ext_cloud;
			g_cloud_layers[i] = (dtau_ray_sw_bandslayers[ib][i] * 0. + dtau_sca_cloud * g_cloud_scaled) / (dtau_ray_sw_bandslayers[ib][i] + dtau_sca_cloud);
			omega0_cloud_layers[i] = (dtau_ray_sw_bandslayers[ib][i] + dtau_sca_cloud) / dtau_cloud_layers[i];
			/* Same thing but without any clouds */
			dtau_plain_layers[i] = dtau_mol_sw_bandslayers[ib][i] + dtau_ray_sw_bandslayers[ib][i];
			g_plain_layers[i] = 0.;
			omega0_plain_layers[i] = dtau_ray_sw_bandslayers[ib][i] / dtau_plain_layers[i];
		}

		doubling_adding_eddington(nlevels, albedo_ground, mu0, wgt_sw_bandslayers[ib], g_plain_layers, omega0_plain_layers, dtau_plain_layers, E_direct_plain_levels, Edn_plain_levels, Eup_plain_levels);
		doubling_adding_eddington(nlevels, albedo_ground, mu0, wgt_sw_bandslayers[ib], g_cloud_layers, omega0_cloud_layers, dtau_cloud_layers, E_direct_cloud_levels, Edn_cloud_levels, Eup_cloud_levels);
		for (int i=0; i < nlevels; i++) {
			total_E_direct_levels[i] += cloud_frac * E_direct_cloud_levels[i] + (1. - cloud_frac) * E_direct_plain_levels[i];
			total_Eup_levels[i] += cloud_frac * Eup_cloud_levels[i] + (1. - cloud_frac) * Eup_plain_levels[i];
			total_Edn_levels[i] += cloud_frac * Edn_cloud_levels[i] + (1. - cloud_frac) * Edn_plain_levels[i];
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
	int nlayers = 22;  /* number of layers */
	int nlevels = nlayers + 1;  /* number of levels */
	double bound_delta_time = 12 * 60 * 60;  /* Bound on the timestep */
	int n_temperature_ter = 5;  /* Number of times the temperature has to be within the threshold */
	double delta_temperature_threshold = 1e-4;  /* Temperature threshold for model termination */

	double p0 = 1000;  /* unit: hPa */
	double delta_p = p0/nlayers;

	double cloud_frac;  // Fraction of clouds when mixing two atmospheres, one with and one without clouds
	double r_cloud = 1e-5;  // Characteristic radius of a droplet; unit: m
	double A_g = 0.12;  // Albedo at the ground
	double mu0 = 0.25;  // Integrate the day-night cycle via a clever parametrization

	double pressure_layers[nlayers];
	double temperature_layers[nlayers];  /* unit: Kelvin */
	double z_layers[nlayers];
	double relative_humidity_layers[nlayers];
	double liquid_water_path_layers[nlayers];

	double pressure_levels[nlevels];
	double temperature_levels[nlevels];  /* unit: Kelvin */
	double z_levels[nlevels];  /* altitude; unit: m */
	double h2ovmr[nlayers], o3vmr[nlayers], co2vmr[nlayers], ch4vmr[nlayers], n2ovmr[nlayers], o2vmr[nlayers], cfc11vmr[nlayers], cfc12vmr[nlayers], cfc22vmr[nlayers], ccl4vmr[nlayers];

	printf("Initializing arrays...\n");
	printf("%5s %12s\n", "Level", "p[hPa]");
	for (int i=0; i < nlevels; i++) {
		pressure_levels[i] = p0/nlayers * i;
		printf("%5d %12.4f\n", i, pressure_levels[i]);
	}
	for (int i=0; i < nlayers; i++) {
		pressure_layers[i] = (pressure_levels[i] + pressure_levels[i+1]) / 2;
		temperature_layers[i] = 288. - ((double) nlayers - (double) i - 1.) * 50./((double) nlayers - 1.);
		liquid_water_path_layers[i] = 0.;
	}

	int nlevels_fpda_file, fpda_status;
	double *h2oppm=NULL, *o3ppm=NULL, *discard1=NULL, *discard2=NULL, *discard3=NULL;
	char fpda_filepath[128] = "ascii/fpda.atm";
	fpda_status = read_5c_file(fpda_filepath, &discard1, &discard2, &discard3, &h2oppm, &o3ppm, &nlevels_fpda_file);
	if (fpda_status != 0) {
		fprintf(stderr, "Error while opening file '%s'. Aborting...\n", fpda_filepath);
		return 1;
	}
	/* Interpolate concentrations for arbitrary layers */
	double pressure_given[nlevels_fpda_file];
	for (int i=0; i < nlevels_fpda_file; i++) {
		pressure_given[i] = p0/(nlevels_fpda_file - 1) * i;
	}
	concentration_interpolation(pressure_given, h2oppm, nlevels_fpda_file, pressure_layers, h2ovmr, nlayers);
	concentration_interpolation(pressure_given, o3ppm, nlevels_fpda_file, pressure_layers, o3vmr, nlayers);

	printf("%5s %12s %12s %12s %12s %12s %12s %12s\n", "Layer", "T[K]", "H2O", "O3", "CO2", "CH4", "N2O", "O2");
	for (int i=0; i < nlayers; i++) {
		co2vmr[i] = 400e-6;
		ch4vmr[i] = 1.7e-6;
		n2ovmr[i] = 320e-9;
		o2vmr[i] = .209;
		cfc11vmr[i] = 0.;
		cfc12vmr[i] = 0.;
		cfc22vmr[i] = 0.;
		ccl4vmr[i] = 0.;

		printf("%5i %12g %12g %12g %12g %12g %12g %12g\n", i, temperature_layers[i], h2ovmr[i], o3vmr[i], co2vmr[i], ch4vmr[i], n2ovmr[i], o2vmr[i]);
	}

	/* Compute the relative humidity profile from the standard atmosphere
	 * in order to keep it constant later on by adapting the VMR of H20.
	 */
	for (int i=0; i < nlayers; i++) {
		/* August-Roche-Magnus formula for saturation water vapor pressure */
		double saturation_h2o_vapor_pressure = 6.1094 * exp(17.625 * (temperature_layers[i] - 273.15) / (temperature_layers[i] - 273.15 + 243.04));
		double saturation_h2o_vmr = saturation_h2o_vapor_pressure / pressure_layers[i];
		relative_humidity_layers[i] = 100 * h2ovmr[i] / saturation_h2o_vmr;
	}

	// Create a cloud
	liquid_water_path_layers[10] = 5.6e-3;
	cloud_frac = .8;

	int nlevels_cloud_prop_lw_file, cloud_prop_lw_status;
	double *wvl_lbound_lw_bands=NULL, *wvl_ubound_lw_bands=NULL, *q_ext_cloud_lw_bands=NULL, *omega0_cloud_lw_bands=NULL, *g_cloud_lw_bands=NULL;
	char cloud_prop_lw_filepath[128] = "rrtm/cldprp/rrtm.lw.int";
	cloud_prop_lw_status = read_5c_file(cloud_prop_lw_filepath, &wvl_lbound_lw_bands, &wvl_ubound_lw_bands, &q_ext_cloud_lw_bands, &omega0_cloud_lw_bands, &g_cloud_lw_bands, &nlevels_cloud_prop_lw_file);
	if (cloud_prop_lw_status != 0) {
		fprintf(stderr, "Error while opening file '%s'. Aborting...\n", cloud_prop_lw_filepath);
		return 1;
	}
	int nlevels_cloud_prop_sw_file, cloud_prop_sw_status;
	double *wvl_lbound_sw_bands=NULL, *wvl_ubound_sw_bands=NULL, *q_ext_cloud_sw_bands=NULL, *omega0_cloud_sw_bands=NULL, *g_cloud_sw_bands=NULL;
	char cloud_prop_sw_filepath[128] = "rrtm/cldprp/rrtm.sw.int";
	cloud_prop_sw_status = read_5c_file(cloud_prop_sw_filepath, &wvl_lbound_sw_bands, &wvl_ubound_sw_bands, &q_ext_cloud_sw_bands, &omega0_cloud_sw_bands, &g_cloud_sw_bands, &nlevels_cloud_prop_sw_file);
	if (cloud_prop_sw_status != 0) {
		fprintf(stderr, "Error while opening file '%s'. Aborting...\n", cloud_prop_sw_filepath);
		return 1;
	}

	gnuplot_ctrl *g1;
	g1 = gnuplot_init();

	printf("Beginning iterative climate modelling...\n");
	int it = 0;
	double model_t = 0.;
	double bound_delta_temperature = 1.;
	int n_temperature = 0;
	double temperature_sum_prev = 0.;
	double temperature_sum_curr = temperature_sum_prev + delta_temperature_threshold + 1.;
	while (n_temperature < n_temperature_ter) {
		it++;
		if (fabs(temperature_sum_curr - temperature_sum_prev) < delta_temperature_threshold) {
			n_temperature++;
		} else {
			n_temperature = 0;
		}
		temperature_sum_prev = temperature_sum_curr;
		temperature_sum_curr = 0.;

		//heating(temperature_layers, delta_t, p0, nlayers);
		//cooling(temperature_layers, delta_t, p0, nlayers);

		/* Compute the flux using fixed bands
		 * ```
		 * double lambda_bands[4] = {1E-20, 8E-6, 12E-6, 1E-03};
		 * double dtau_windows[3] = {1.6, 0., 1.6};
		 * band_flux(nlevels, 4, albedo, lambda_bands, dtau_windows, temperature_layers, total_Edn_levels, total_Eup_levels);
		 * ```
		 */

		double total_Eup_levels[nlevels], total_Edn_levels[nlevels];
		double total_Eup_lw_levels[nlevels], total_Edn_lw_levels[nlevels];
		double total_E_direct_levels[nlevels], total_Eup_sw_levels[nlevels], total_Edn_sw_levels[nlevels];

		/* Keep the relative humidity constant by adapting the VMR of H2O */
		for (int i=0; i < nlayers; i++) {
			/* August-Roche-Magnus formula for saturation water vapor pressure */
			double saturation_h2o_vapor_pressure = 6.1094 * exp(17.625 * (temperature_layers[i] - 273.15) / (temperature_layers[i] - 273.15 + 243.04));
			double saturation_h2o_vmr = saturation_h2o_vapor_pressure / pressure_layers[i];
			h2ovmr[i] = relative_humidity_layers[i] * saturation_h2o_vmr / 100;
		}

		/* Compute the energy fluxes using the RRTM */
		rrtm_lw_flux(nlayers, A_g, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, liquid_water_path_layers, cloud_frac, r_cloud, wvl_lbound_lw_bands, wvl_ubound_lw_bands, q_ext_cloud_lw_bands, omega0_cloud_lw_bands, nlevels_cloud_prop_lw_file, total_Edn_lw_levels, total_Eup_lw_levels);
		rrtm_sw_flux(nlayers, A_g, mu0, pressure_levels, temperature_layers, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, liquid_water_path_layers, cloud_frac, r_cloud, wvl_lbound_sw_bands, wvl_ubound_sw_bands, q_ext_cloud_sw_bands, omega0_cloud_sw_bands, g_cloud_sw_bands, nlevels_cloud_prop_sw_file, total_E_direct_levels, total_Edn_sw_levels, total_Eup_sw_levels);
		for (int i=0; i < nlevels; i++) {
			total_Eup_levels[i] = total_Eup_lw_levels[i] + total_Eup_sw_levels[i];
			total_Edn_levels[i] = total_E_direct_levels[i] + total_Edn_lw_levels[i] + total_Edn_sw_levels[i];
		}

		/* Calculate energy differences in order to adapt the temperature */
		double div_E_layers[nlayers];
		for (int i=0; i < nlayers; i++) {
			div_E_layers[i] = (total_Edn_levels[i] - total_Eup_levels[i]) - (total_Edn_levels[i+1] - total_Eup_levels[i+1]);
		}
		div_E_layers[nlayers-1] += total_Edn_levels[nlevels-1] - total_Eup_levels[nlevels-1];
		/* Calculate an adaptive time difference and add it to the runtime */
		double max_E_net = 0.;
		for (int i=0; i < nlayers; i++) {
			max_E_net = max(fabs(div_E_layers[i]), max_E_net);
		}
		/* Convert between hPa and Pa by multiplying delta_p with 100 */
		double delta_t = min(C_P/G * (100 * delta_p)/max_E_net * bound_delta_temperature, bound_delta_time);
		model_t += delta_t;
		/* Actually adapt the temperature after having defined the time step; also see previous loop */
		for (int i=0; i < nlayers; i++) {
			temperature_layers[i] += div_E_layers[i] * G / (100 * delta_p * C_P) * delta_t;
			temperature_sum_curr += temperature_layers[i];
		}

		/* Sort layers by potential temperature, a.k.a. do convection */
		convection(temperature_layers, pressure_layers, nlayers);

		/* Textual and visual output */
		if (it%5 == 0) {
			for (int i=0; i < nlevels; i++) {
				z_levels[i] = barometric_PToZ(pressure_levels[i], temperature_layers[nlayers-1], p0);
			}
			for (int i=0; i < nlayers; i++) {
				z_layers[i] = (z_levels[i] + z_levels[i+1]) / 2.;
			}

			printf("Iteration %7d at time %6.2fd (%8.1fs)\n", it, model_t / (60*60*24), model_t);
			printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "Edn[W]", "Eup[W]", "Edn_lw[W]", "Eup_lw[W]", "E_dir[W]", "Edn_sw[W]", "Eup_sw[W]", "z_lyr[m]", "T_lyr[K]");
			for (int i=0; i < nlayers; i++) {
				printf("%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.1f %12.4f\n", total_Edn_levels[i], total_Eup_levels[i], total_Edn_lw_levels[i], total_Eup_lw_levels[i], total_E_direct_levels[i], total_Edn_sw_levels[i], total_Eup_sw_levels[i], z_layers[i], temperature_layers[i]);
			}
			printf("%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.1f %12.4f\n", total_Edn_levels[nlevels-1], total_Eup_levels[nlevels-1], total_Edn_lw_levels[nlevels-1], total_Eup_lw_levels[nlevels-1], total_E_direct_levels[nlevels-1], total_Edn_sw_levels[nlevels-1], total_Eup_sw_levels[nlevels-1], (double) NAN, (double) NAN);

			gnuplot_resetplot(g1);  /* Start with a new plot rather than plotting into existing one */
			gnuplot_setstyle(g1, "linespoints");  /* Draw lines and points */
			gnuplot_set_xlabel(g1, "temperature [K]");  /* x-axis label */
			gnuplot_set_ylabel(g1, "altitude [m]");  /* y-axis label */

			/* Plot temperature T as function of z and label with temperature */
			gnuplot_plot_xy(g1, temperature_layers, z_layers, nlayers, "Temperature") ;
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

	printf("Modelled results by levels...\n");
	printf("%5s %12s %12s %12s\n", "Level", "p[hPa]", "z[m]", "T[K]");
	for (int i=0; i < nlevels; i++) {
		printf("%5d %12.4f %12.4f %12.4f\n", i, pressure_levels[i], z_levels[i], temperature_levels[i]);
	}
	printf("Modelled results by layers...\n");
	printf("%5s %12s %12s %12s\n", "Layer", "p[hPa]", "z[m]", "T[K]");
	for (int i=0; i < nlayers; i++) {
		printf("%5d %12.4f %12.4f %12.4f\n", i, pressure_layers[i], z_layers[i], temperature_layers[i]);
	}

	return 0;
}
