#ifndef _ANALYSE_
#define _ANALYSE_
void analyse(t_non *non);
float calc_participation_ratio(int N,float *H);
float calc_local_participation_ratio(int N,float *H,float min,float max,float *e,float shift);
int find_cEig(float *cEig,float *cDOS,float *dip2,float *H,float *e,int N,float min,float max,int counts,float shift);
void find_dipole_mag(t_non *non,float *dip2,int step,FILE *mu_traj,float *H);
#endif // _ANALYSE_
