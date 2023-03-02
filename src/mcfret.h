#ifndef _MCFRET_
#define _MCFRET_
void mcfret(t_non *non);
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non);
void mcfret_autodetect(t_non *non, float treshold);
void mcfret_response_function(float* Hamiltonian, float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci, float *density_matrix, int emission);
void mcfret_coupling(t_non *non);
void mcfret_rate(t_non *non);
void mcfret_validate(t_non *non);
void mcfret_analyse(t_non *non);

#endif /* _MCFRET_ */
