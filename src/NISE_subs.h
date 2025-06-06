#ifndef _NONSUBS_
#define _NONSUBS_

#include "lapack.h"
void** calloc2D(size_t nRows, size_t nCols, size_t size, size_t sizeP);
void free2D(void** arr);
void not_parallel();
void copyvec(float *a,float *b,int N);
void clearvec(float *a,int N);
void clearvec_double(double* a, int N);
void unitmat(float *a,int N);
void vector_on_vector(float *rr,float *ir,float *vr,float *vi,int N);
void matrix_on_vector(float *c,float *vr,float *vi,int N);
void trans_matrix_on_vector(float *c,float *vr,float *vi,int N);
float find_norm(float *phi_r,float *phi_i,int N);
void re_normalize(float *phi_r,float *phi_i,int N,float norm);
void log_item(char* msgFormat, ...);
time_t set_time(time_t t0);
time_t log_time(time_t t0,FILE *log);
int string_in_array(char* string_to_compare, char* string_array[], int array_size);
int determine_samples (t_non *non);
char* time_diff(time_t t0, time_t t1);
char* MPI_time(double t0);
int Eindex(int a,int b,int N);
int read_He(t_non *non,float *He,FILE *FH,int pos);
int read_Dia(t_non *non,float *He,FILE *FE,int pos);
int read_A(t_non *non,float *Anh,FILE *FH,int pos);
int read_mue(t_non *non,float *mue,FILE *FH,int pos,int x);
int read_alpha(t_non *non,float *alpha,FILE *FH,int pos,int x);
int read_over(t_non *non,float *over,FILE *FH,int pos,int x);
void muread(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj);
void mureadE(t_non *non,float *leftnr,int ti,int x,FILE *mu_traj,float *mu,float *pol);
int read_cluster(t_non *non,int pos,int *cl,FILE *FH);
void propagate_vec_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int sign);
void generateCS(float *X,float *Y,float *Z);
//void projection(float *phi,t_non *non);
int control(t_non *non);
int autodetect_singles(t_non* non);
void dipole_double(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
void dipole_double_ground(t_non *non,float *dipole,float *fr,float *fi,float *over);
void dipole_double_ES(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi);
void dipole_double_last(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi,float *over);
void dipole_double_last_ES(t_non *non,float *dipole,float *cr,float *ci,float *fr,float *fi);
float distance(float *rf,float *ri,int a,int b,int N,float box);
float distance_x(float *rf,float *ri,int a,int b,int N,float box,int x);
float distance3(float *rf,float *ri,int a,int b,int N,float *box);
float distance3_x(float *rf,float *ri,int a,int b,int N,float *box,int x);
float pbc1(float r, int x, float *box);
void diagonalize_real_nonsym(float* K, float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N);


// Index triangular matrix
// Put in the .h file to allow external referencing
inline int Sindex(int a, int b, int N) { // inline to make it quicker
    int ind;
    if (a > b) {
        //ind=a+N*b-(b*(b+1)/2);
        ind = a + b * ((N << 1) - b - 1) / 2;
    }
    else {
        //ind=b+N*a-(a*(a+1)/2);
        ind = b + a * ((N << 1) - a - 1) / 2;
    }
    return ind;
}

#endif // _NONSUBS_
