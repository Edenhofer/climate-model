#ifndef __C_TR
#define __C_TR 1

double planck(double lambda, double T);
double planck_int(double lambda_lower, double lambda_upper, double T);
void schwarzschild(int nlevels, double albedo, double *dtau_layers, double *B_layers, double B_surface, double *Edn_levels, double *Eup_levels);

#endif
