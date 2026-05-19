#ifndef _MCFRET_ /* ignore */
#define _MCFRET_

void mcfret(t_non *non);
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int segments, float *partition_functions);
void mcfret_autodetect(t_non *non, float treshold);
void mcfret_response_function(float *re_S_1,float *im_S_1,t_non *non,int emission, float *ave_vecr);
void mcfret_coupling(float *J,t_non *non);
void mcfret_energy(float *E,t_non *non,int segments, float *ave_vecr,float *energy_cor);
void mcfret_rate(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs,float *re_Emi,float *im_Emi,float *J,t_non *non);
void mcfret_validate(t_non *non);
void mcfret_eigen(t_non *non,float *rate_matrix,float *re_e,float *im_e,float *vl,float *vr,int segments,float *energy_cor);
void mcfret_analyse(float *E,float *rate_matrix,t_non *non,int segments);
void mcfret_response_function_sub(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci);
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,float *rC,float *iC,int *psites,int segments,int si,int sj,int sk,int N);
float trace_rate(float *matrix,int N);
void write_matrix_to_file_float(char fname[],float *matrix,int N);
void read_response_from_file(char fname[],float *re_R,float *im_R,int N,int tmax);
void triangular_on_square(float *T,float *S,int N);
void average_density_matrix(float *ave_den_mat,t_non *non);

void mcfret_propagation_segmented(float *re_S_1,float *im_S_1,t_non *non);
void mcfret_response_function_sub_segments(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci, int *H_indices_si,int N_i);

void hermitian_conjugate(float *A_re, float *A_im, float *hermi_re, float *hermi_im, int N1, int N2);
void mcfret_rate_from_abs(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs, float *rho_0,float *J,t_non *non);

#endif /* _MCFRET_ */
