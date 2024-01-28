#ifndef _PROJECT_
#define _PROJECT_

void projection(float *phi,t_non *non);
void multi_projection(float *phi_in,float *phi_out,t_non *non,int ip);
void multi_projection_Hamiltonian(float *Hamil_i_e, t_non *non);
void multi_projection_Coupling(float *Hamil_i_e, t_non *non);
int project_dim(t_non* non);
void project_degeneracies(t_non* non,int *degen,int segments);


#endif // _PROJECT_
