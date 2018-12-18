#ifndef __C_TR
#define __C_TR 1

double planck(double lambda, double T);
double planck_int(double lambda_lower, double lambda_upper, double T);
void schwarzschild(int nlevels, double albedo, double *dtau_layers, double *B_layers, double B_surface, double *Edn_levels, double *Eup_levels);
void doubling_adding(int nlevels, double albedo_ground, double *r_layers, double *t_layers, double *t_direct_layers, double *s_direct_layers, double *r_direct_layers, double E_direct_toa, double *E_direct_levels, double *Edn_levels, double *Eup_levels);
void doubling_adding_eddington(int nlevels, double albedo_ground, double mu0, double E_direct_toa, double *g_layers, double *omega0_layers, double *dtau, double *E_direct_levels, double *Edn_levels, double *Eup_levels);

#endif
