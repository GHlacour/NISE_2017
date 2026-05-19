#ifndef _PROPSEG_
#define _PROPSEG_

void propagate_vector_segments(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display, int N_i);
void propagate_matrix_segments(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display, int N_i);
void propagate_vec_coupling_S_segments(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, int sign, int N_i);
void propagate_snapshot(float *U_snap_re, float *U_snap_im, float **U_comp_re, float **U_comp_im, float **temp_re, float **temp_im, int N);


#endif // _PROPSEG_