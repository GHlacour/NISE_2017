#ifndef _MCFRET_ /* ignore */
#define _MCFRET_
void mcfret(t_non *non);
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int segments);
void mcfret_autodetect(t_non *non, float treshold);
void mcfret_response_function(double *re_S_1,double *im_S_1,t_non *non,int emission);
void mcfret_coupling(double *J,t_non *non);
void mcfret_energy(double *E,t_non *non,int segments);
void mcfret_rate(double *rate_matrix,double *coherence_matrix,int segments,double *re_Abs,double *im_Abs,double *re_Emi,double *im_Emi,double *J,t_non *non);
void mcfret_validate(t_non *non);
void mcfret_analyse(double *E,double *rate_matrix,t_non *non,int segments);
void mcfret_response_function_sub(double *re_S_1,double *im_S_1,int t1,t_non *non,float *cr,float *ci);
void segment_matrix_mul(double *rA,double *iA,double *rB,double *iB,double *rC,double *iC,int *psites,int segments,int si,int sj,int sk,int N);
double trace_rate(double *matrix,int N);
void integrate_rate_response(double *rate_response,int T,double *is13,double *isimple);
void write_matrix_to_file(char fname[],float *matrix,int N);
void write_matrix_to_file_double(char fname[],double *matrix,int N);
void read_matrix_from_file(char fname[],double *matrix,int N);
void read_response_from_file(char fname[],double *re_R,double *im_R,int N,int tmax);
void triangular_on_square(float *T,float *S,int N);
#endif /* _MCFRET_ */
