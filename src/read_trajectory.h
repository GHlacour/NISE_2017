#ifndef _TRAJECTORY_
#define _TRAJECTORY_

void open_files(t_non *non,FILE **H_traj,FILE **mu_traj,FILE **Cfile);
void read_coupling(t_non *non,FILE *C_traj,FILE *mu_traj,float* Hamil_i_e,float *mu_xyz);
void read_Hamiltonian(t_non *non,float *Hamil_i_e,FILE *H_traj,int pos);
void read_dipole(t_non *non,FILE *mu_traj,float *mu_eg,float *mu_xyz,int x,int pos);

int read_He(t_non *non,float *He,FILE *FH,int pos);
int read_Dia(t_non *non,float *He,FILE *FE,int pos);
int read_A(t_non *non,float *Anh,FILE *FH,int pos);
int read_mue(t_non *non,float *mue,FILE *FH,int pos,int x);
int read_over(t_non *non,float *over,FILE *FH,int pos,int x);
void muread(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj);
void mureadE(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj,float *mu,float *pol);
int read_cluster(t_non *non,int pos,int *cl,FILE *FH);

#endif // _TRAJECTORY_
