#ifndef _ANALYSE_
#define _ANALYSE_
void analyse(t_non *non);
float calc_participation_ratio(int N,float *H);
float calc_local_participation_ratio(int N,float *H,float min,float max,float *e,float shift);
float calc_spectral_participation_ratio(int N,float *H);
float calc_local_spectral_participation_ratio(int N,float *H,float min,float max,float *e,float shift);
int find_cEig(float *cEig,float *cDOS,float *dip2,float *H,float *e,int N,float min,float max,int counts,float shift);
void find_dipole_mag(t_non *non,float *dip2,int step,FILE *mu_traj,float *H,float *my_xyz);
void calc_densitymatrix(t_non *non,float *rho,float *rho2,float *rho4,float *local_rho,float *spec_rho,float *H,float* e,float *dip2);
void cluster(t_non *non,float *rho);
#endif // _ANALYSE_
