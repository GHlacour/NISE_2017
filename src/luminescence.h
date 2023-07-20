#ifndef _LUMINESCENCE_
#define _LUMINESCENCE_
void luminescence(t_non *non);
void calc_LUM(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu);
void bltz_weight(float *mu_eg,float *Hamil_i_e,t_non *non);
void bltz_weight_itime(float *cr,float *Hamiltonian_i,t_non *non);
#endif // _LUMINESCENCE_
