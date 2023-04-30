#ifndef _MCFRET_ /* ignore */
#define _MCFRET_
void mcfret(t_non *non);
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int segments);
void mcfret_autodetect(t_non *non, float treshold);
void mcfret_response_function(float *re_S_1,float *im_S_1,t_non *non,int emission);
void mcfret_coupling(float *J,t_non *non);
void mcfret_rate(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs,float *re_Emi,float *im_Emi,float *J,t_non *non);
void mcfret_validate(t_non *non);
void mcfret_analyse(t_non *non);
void mcfret_response_function_sub(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci);
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,float *rC,float *iC,int *psites,int segments,int si,int sj,int sk,int N);
float trace_rate(float *matrix,int N);
void integrate_rate_response(float *rate_response,int T,float *is13,float *isimple);
float write_matrix_to_file(char fname[],float *matrix,int N);
#endif /* _MCFRET_ */
