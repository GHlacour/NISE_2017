#ifndef _PROPAGATE_
#define _PROPAGATE_
#include "lapack.h"
void propagate_vector(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display);
void propagate_matrix(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display);

void propagate_vec_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign);
void propagate_t2_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,float **vr,float **vi,int sign);
int propagate_vec_DIA_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign);
void propagate_vec_RK4(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign);
void propagate_vec_RK4_doubles(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,float *Anh);
void propagate_vec_RK4_doubles_ES(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m);
void propagate_vec_coupling_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign);
void propagate_vec_coupling_S_doubles(t_non *non,float *Hamiltonian_i,float *cr,float
*ci,int m,float *Anh);
void propagate_vec_coupling_S_doubles_ES(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m);
void diagonalizeLPD(float *H,float *v,int N);
void build_diag_H(float *Hamiltonian_i,float *H,float *e,int N);
int time_evolution_mat(t_non *non,float *Hamiltonian_i,float *Ur,float *Ui,int *R,int *C,int m);
void propagate_double_sparce(t_non *non,float *Ur,float *Ui,int *R,int *C,float *fr,float *fi,int elements,int m,float *Anh);
void propagate_double_sparce_ES(t_non *non,float *Ur,float *Ui,int *R,int *C,float *fr,float *fi,int elements,int m);
void propagate_double_2DIR_control(t_non* non, float* Hamil_i_e, float** ft1r, float** ft1i,
                                  float* fr, float* fi, float** rightrr, float** rightri,
                                  float* rightnr, float* rightni, int currentSample,
                                  int molPol, int t3, float* Anh);
void propagate_double_2DES_control(t_non* non, float* Hamil_i_e, float** ft1r, float** ft1i,
                                  float* fr, float* fi, float** rightrr, float** rightri,
                                  float* rightnr, float* rightni, int currentSample,
                                  int molPol, int t3);
#endif // _PROPAGATE_
