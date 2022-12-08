#ifndef _PROJECT_
#define _PROJECT_

void projection(float *phi,t_non *non);
void multi_projection(float *phi_in,float *phi_out,t_non *non,int ip);
void multiprojection (float *Hamil_i_e, t_non *non, int ip)
int project_dim(t_non* non);


#endif // _PROJECT_

