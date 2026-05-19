#ifndef _PROJECT_
#define _PROJECT_

void projection(float *phi,t_non *non);
void multi_projection(float *phi_in,float *phi_out,t_non *non,int ip);
void multi_projection_Hamiltonian(float *Hamil_i_e, t_non *non);
void multi_projection_Coupling(float *Hamil_i_e, t_non *non);
void zero_coupling(float *Hamil_i_e, t_non *non);
int project_dim(t_non* non);
void project_degeneracies(t_non* non,int *degen,int segments);
int find_H_indices_segment(int *psites, int *H_indices_si,int si, t_non *non);
void isolate_segment_Hamiltonian_triu(float *Hamiltonian_full_triu, float *Hamiltonian_segment_triu, int *H_indices_si, int N_i, t_non *non);
void isolate_segment_Hamiltonian(float *Hamiltonian_full, float *Hamiltonian_segment, int *H_indices_si, int N_i, t_non *non);
void isolate_coupling_block(float *J_full, float *J_ij, int N_i, int N_j, int *H_indices_si, int *H_indices_sj, t_non *non);
int find_max_segment_size(int *psites, t_non *non);


#endif // _PROJECT_
