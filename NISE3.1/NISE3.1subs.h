#ifndef _NONSUBS_
#define _NONSUBS_

#include "lapack.h"
void copyvec(float *a,float *b,int N);
void clearvec(float *a,int N);
time_t set_time(time_t t0);
time_t log_time(time_t t0,FILE *log);
int Sindex(int a,int b,int N);
int Eindex(int a,int b,int N);
int read_He(t_non *non,float *He,FILE *FH,int pos);
int read_A(t_non *non,float *Anh,FILE *FH,int pos);
int read_mue(t_non *non,float *mue,FILE *FH,int pos,int x);
int read_over(t_non *non,float *over,FILE *FH,int pos,int x);
void muread(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj);
int read_cluster(t_non *non,int pos,int *cl,FILE *FH);
void propagate_vec_DIA(t_non *non,float *H,float *cr,float *ci,int sign);
int propagate_vec_DIA_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign);
void propagate_vec_coupling_S(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign);
void propagate_vec_coupling_S_doubles(t_non *non,float *Hamiltonian_i,float *cr,float 
*ci,int m,float *Anh);
void diagonalizeLPD(float *H,float *v,int N);
void build_diag_H(float *Hamiltonian_i,float *H,float *e,int N);
void generateCS(float *X,float *Y,float *Z);
void projection(float *phi,t_non *non);
int control(t_non *non);
void dipole_double(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
void dipole_double_last(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
int time_evolution_mat(t_non *non,float *Hamiltonian_i,float *Ur,float *Ui,int *R,int *C,int m);
void propagate_double_sparce(t_non *non,float *Ur,float *Ui,int *R,int *C,float *fr,float *fi,int elements,int m,float *Anh);

#endif // _NONSUBS_
