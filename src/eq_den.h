#ifndef _EQ_DEN_
#define _EQ_DEN_
#include "lapack.h"
void eq_den(float *Hamiltonian_i, float *rho_l, int N, t_non *non);
void matrix_on_real_vector(float *mat,float *vr,int N);
void dipole_double_inverse_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi);
void dipole_double_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi);
void diagonalize_real_nonsym(float* K, float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N);

void pe_pa_cr(int t1, int t2, int t3, float polWeight, 
float up_ver3_SE_re_NR, float up_ver3_SE_re_R, float up_ver3_GB_re_NR, float up_ver3_GB_re_R,
 float up_ver3_EA_re_NR, float up_ver3_EA_re_R,
float up_ver3_SE_im_NR , float up_ver3_GB_im_NR,   float up_ver3_EA_im_NR,   
float up_ver3_SE_im_R,   float up_ver3_GB_im_R,    float up_ver3_EA_im_R,
float *re_2DES_NR_pe,float *re_2DES_R_pe, float *im_2DES_NR_pe,  float *im_2DES_R_pe,  t_non *non);

#endif // _EQ_DEN_
