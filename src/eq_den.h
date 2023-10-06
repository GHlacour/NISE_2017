#ifndef _EQ_DEN_
#define _EQ_DEN_
#include "lapack.h"
void eq_den(float *Hamiltonian_i, float *rho_l, int N, t_non *non);
void matrix_on_real_vector(float *mat,float *vr,int N);
void dipole_double_inverse_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi);
void dipole_double_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi);
void diagonalize_real_nonsym(float* K, float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N);
void write_matrix_to_file(char fname[],float *matrix,int N);
void write_ham_to_file(char fname[],float *Hamiltonian_i,int N);

#endif // _EQ_DEN_
