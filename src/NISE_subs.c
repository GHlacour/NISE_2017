#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "read_trajectory.h"
#include "randomlib.h"
#include "util/asprintf.h"
#include <ctype.h>
#include <cblas.h>

// Subroutines for nonadiabatic code

// Allocate 2D memory blocks
void** calloc2D(size_t nRows, size_t nCols, size_t size, size_t sizeP) {
    void** result = malloc(nRows * sizeP);
    void* data = calloc(nRows * nCols, size);
    for(int i = 0; i < nRows; i++) {
        result[i] = (unsigned char*) data + size * i * nCols;
    }

    return result;
}

// Free 2D memory blocks
void free2D(void** arr) {
    free(arr[0]);
    free(arr);
}

/* Inform user that they are trying to run code in parallel that is not parallelized */
void not_parallel(){
    printf(RED "This part of the code is not parallel!!!\n");
    printf("You may waste valuable computational resources.\n");
    printf("If running interactively you are fine, if you run\n");
    printf("trough a queue on a computer cluster, please,\n");
    printf("consider running in serial by asking for only one CPU.\n\n" RESET);
}

// Copy a vector
void copyvec(float* a, float* b, int N) {
    int i;
    for (i = 0; i < N; i++) b[i] = a[i];
}

// Set all elements of a vector to zero
void clearvec(float* a, int N) {
    int i;
    for (i = 0; i < N; i++) a[i] = 0;
}
// Set all elements of a vector to zero
void clearvec_double(double* a, int N) {
    int i;
    for (i = 0; i < N; i++) a[i] = 0;
}

// Set all elements of a vector to zero, here for integers
void clearvec_int(int *a, int N) {
    int i;
    for (i = 0; i < N; i++) a[i] = 0;
}


// Construct unit matrix
void unitmat(float *a,int N){
    int i,j;
    for (i = 0; i < N; i++){
	for (j = 0; j < N ; j++){
            a[i+N*j]=0;
    	    if (i==j) a[i+N*j]=1;
	}
    }
}

// Multiply a complex diagonal matrix on a complex vector
void vector_on_vector(float *rr,float *ir,float *vr,float *vi,int N){
    int a;
    float re,im;
    for (a=0;a<N;a++){
	re=rr[a]*vr[a]-ir[a]*vi[a];
	im=ir[a]*vr[a]+rr[a]*vi[a];
	vr[a]=re;
	vi[a]=im;
    }
}

// Multiply a real matrix on a complex vector (vr,vi)
void matrix_on_vector(float *c,float *vr,float *vi,int N){
    float *xr;
    float *xi;
    int a,b;
    xr = (float *)calloc(N * N, sizeof(float));
    xi = (float *)calloc(N * N, sizeof(float));
    // Multiply
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            xr[a]+=c[a+b*N]*vr[b];
	    xi[a]+=c[a+b*N]*vi[b];
	}
    }
    // Copy back
    copyvec(xr,vr,N);
    copyvec(xi,vi,N);
    free(xr);
    free(xi);
}

// Multiply transpose of a real matrix on a complex vector (vr,vi)
void trans_matrix_on_vector(float *c,float *vr,float *vi,int N){
    float *xr;
    float *xi;
    int a,b;
    xr = (float *)calloc(N * N, sizeof(float));
    xi = (float *)calloc(N * N, sizeof(float));
    // Multiply
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            xr[a]+=c[b+a*N]*vr[b];
            xi[a]+=c[b+a*N]*vi[b];
        }
    }
    // Copy back
    copyvec(xr,vr,N);
    copyvec(xi,vi,N);
    free(xr);
    free(xi);
}
/* Find the trace of a matrix product */
float matrix_mul_traced_DA(float *A, float *B, int N_i, int N_i3){
    /* Calculate the trace of the product of two differently shaped matrices of dimension N_i * Ni_3 */
    /* Directly computing the trace greatly reduces the number of operations needed */
    /* from N^3 to N^2 */
    // N_i3 is the shared dimension of the matrices A & B
    // N_i is the dimension of the resulting square matrix
    int i,i3;
    float the_trace;

    the_trace = 0;

    for (i=0;i<N_i;i++){
   	    for (i3=0;i3<N_i3;i3++){
            the_trace += A[i*N_i3+i3] * B[i3*N_i+i];
	    }
    }
    return the_trace;
}

// void complex_matrix_product(float *A_re, float *A_im, float *B_re, float *B_im, float *C_re, float *C_im,int N_1,int N_2,int N_3){

//     int i1, i2, i3;
//     float Aim_i1i3, Are_i1i3;
//     clearvec(C_re,N_1*N_2);
//     clearvec(C_im,N_1*N_2);

// // can make paralllel but check for loop order: warnings recieved
// #pragma omp parallel for private(i2,i3,Aim_i1i3,Are_i1i3)
// for (i1=0;i1<N_1;i1++){
//     for (i3=0;i3<N_3;i3++){
// 	    Aim_i1i3 = A_im[i1*N_3+i3];
// 	    Are_i1i3 = A_re[i1*N_3+i3];
//             for (i2=0;i2<N_2;i2++){
//                     C_re[i1*N_2+i2] += Are_i1i3 * B_re[i3*N_2+i2] - Aim_i1i3 * B_im[i3*N_2+i2];
//                     C_im[i1*N_2+i2] += Are_i1i3 * B_im[i3*N_2+i2] + Aim_i1i3 * B_re[i3*N_2+i2];
//                 }
//             }
//         }
// }

/* Computes the complex matrix product of two matrices, */
/* After converting to complex arrays */
/* Only recommended for large matrices (i.e. N >1000) */
void complex_matrix_product(float *A_re, float *A_im,
                            float *B_re, float *B_im,
                            float *C_re, float *C_im,
                            int N_1, int N_2, int N_3)
{
    // Row-major matrices:
    // A: N_1 x N_3
    // B: N_3 x N_2
    // C: N_1 x N_2

    const int sizeA = N_1 * N_3;
    const int sizeB = N_3 * N_2;
    const int sizeC = N_1 * N_2;

    // Allocate interleaved complex buffers
    float *A = (float*) malloc(sizeof(float) * 2 * sizeA);
    float *B = (float*) malloc(sizeof(float) * 2 * sizeB);
    float *C = (float*) malloc(sizeof(float) * 2 * sizeC);

    if (!A || !B || !C) {
        free(A);
        free(B);
        free(C);
        return; // allocation failed
    }

    // Convert A to interleaved complex
    for (int i = 0; i < sizeA; i++) {
        A[2*i]     = A_re[i];
        A[2*i + 1] = A_im[i];
    }

    // Convert B
    for (int i = 0; i < sizeB; i++) {
        B[2*i]     = B_re[i];
        B[2*i + 1] = B_im[i];
    }

    // Zero C
    memset(C, 0, sizeof(float) * 2 * sizeC);

    const float alpha[2] = {1.0f, 0.0f};
    const float beta[2]  = {0.0f, 0.0f};

    cblas_cgemm(CblasRowMajor,
                CblasNoTrans, CblasNoTrans,
                N_1, N_2, N_3,
                alpha,
                A, N_3,
                B, N_2,
                beta,
                C, N_2);

    // Convert result back to split format
    for (int i = 0; i < sizeC; i++) {
        C_re[i] = C[2*i];
        C_im[i] = C[2*i + 1];
    }

    free(A);
    free(B);
    free(C);
}

/* Find the norm for a complex vector */
float find_norm(float *phi_r,float *phi_i,int N){
    float norm;
    int a;
    norm=0;
    for (a=0;a<N;a++){
        norm+=phi_r[a]*phi_r[a]+phi_i[a]*phi_i[a];
    }
    return norm;
}

/* Find new norm of a complex vector and renormalize to preserve old norm */
void re_normalize(float *phi_r,float *phi_i,int N,float norm){
    int a;
    float new_norm,factor;
    new_norm=find_norm(phi_r,phi_i,N);
    if (norm!=new_norm){
        factor=sqrt(norm/new_norm);
        for (a=0;a<N;a++){
            phi_r[a]=phi_r[a]*factor;
            phi_i[a]=phi_i[a]*factor;
    
        }
    }
}

/* For testing purposes */
/* Find the sum of all matrix elements */
float matrix_sum(float *matrix,int N){
    int i,j;
    float sum;
    sum=0;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            sum=sum+matrix[N*i+j];
            // printf("element %f\n",matrix[N*i+j]);
        }
    }
    return sum;
}
/* For testing purposes */


/**
 * Method that logs a message, in which the message can be formatted like printf accepts.
 */
void log_item(char* msgFormat, ...) {
    // Parse parameters
    va_list args;
    va_start(args, msgFormat);

    // Write to log
    FILE* log = fopen("NISE.log", "a");
    if (log == NULL) {
        printf("Could not open log file!");
        exit(1);
    }
    vfprintf(log, msgFormat, args);
    fclose(log);

    va_end(args);
}

// Set time and write to screen
time_t set_time(time_t t0) {
    return log_time(t0, stdout);
}

// Set time and write to log file
time_t log_time(time_t t0, FILE* log) {
    time_t t1;
    time(&t1);

    char* text = time_diff(t0, t1);
    // fprintf(log, text);
    free(text);
    return t1;
}

/* Compare a string to an array of options */
/* The comparizon is not case sensitive */
int string_in_array(char* string_to_compare, char* string_array[], int array_size) {
    for (int i = 0; i < array_size; i++) {
        if (!strcmp_nocase(string_to_compare, string_array[i])) {
            return i+1; // return the index of the matched string
	    }
    }
    return 0; // if no match is found, return 0
}

/* Determine number of samples to use and write to log file */
int determine_samples (t_non *non){
  FILE *log;
  int N_samples;
  N_samples=(non->length-non->tmax1-1)/non->sample+1;
  if (N_samples>0) {
    printf("Making %d samples!\n",N_samples);
  } else {
    printf(RED "Insufficient data to calculate spectrum.\n" RESET);
    printf(RED "Please, lower max times or provide longer\n" RESET);
    printf(RED "trajectory.\n" RESET);
    exit(1);
  }

  if (non->end==0) non->end=N_samples;
  if (non->end>N_samples){
    printf(RED "Endpoint larger than number of samples was specified.\n" RESET);
    printf(RED "Endpoint was %d but cannot be larger than %d.\n" RESET,non->end,N_samples);
    exit(0);
  }
  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);

  return N_samples;
}

// Forms a string with the time difference between the given times
char* time_diff(time_t t0, time_t t1) {
    int control;
    int s = difftime(t1, t0);
    int h = s / 3600;
    s = s % 3600;
    int m = s / 60;
    s = s % 60;

    char* text;
    control=asprintf(&text, "Time spent: %dh %dmin %ds\n", h, m, s);
    return text;
}

// Forms a string with the times for MPI_Wtime
char* MPI_time(double t0) {
    int control;
    int ms =t0*1000; // Convert to milliseconds
    int h = ms / 3600000;
    ms = ms % 3600000;
    int m = ms / 60000;
    ms = ms % 60000;
    int s = ms / 1000;
    ms = ms % 1000;

    char* text;
    control=asprintf(&text, " %dh %dmin %ds %dms\n", h, m, s, ms);
    return text;
}

/* INDEXING FOR ELECTRONIC STATES */
int Eindex(int a, int b, int N) {
    int ind;
    if (a > b) {
        ind = a  + (b - 1) * (N + N - b - 2) / 2;
    }
    else {
        ind = b  + (a - 1) * (N + N - a - 2) / 2;
    }
//    printf("%d %d %d %d\n",a,b,ind,N);
    return ind;
}



// This subroutine generates vectors with a coordinate transformation
// corresponding to a randomly selected isotropic orientation
void generateCS(float* X, float* Y, float* Z) {
    int no;
    // Generate X vector
    X[0] = RandomGaussian(0, 1);
    X[1] = RandomGaussian(0, 1);
    X[2] = RandomGaussian(0, 1);
    // Normalize
    no = 1 / sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
    X[0] = X[0] / no;
    X[1] = X[1] / no;
    X[2] = X[2] / no;
    // Generate Y vector
    Y[0] = RandomGaussian(0, 1);
    Y[1] = RandomGaussian(0, 1);
    Y[2] = RandomGaussian(0, 1);
    // Make it orthogonal to X
    no = X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2];
    Y[0] = Y[0] - no * X[0];
    Y[1] = Y[1] - no * X[1];
    Y[2] = Y[2] - no * X[2];
    // Normalize Y
    no = 1 / sqrt(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2]);
    Y[0] = Y[0] / no;
    Y[1] = Y[1] / no;
    Y[2] = Y[2] / no;
    // Generate Z by taking cross product
    Z[0] = X[1] * Y[2] - X[2] * Y[1];
    Z[1] = -X[0] * Y[2] + X[2] * Y[0];
    Z[2] = X[0] * Y[1] - X[1] * Y[0];
    // Normalize Z
    no = 1 / sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
    Z[0] = Z[0] / no;
    Z[1] = Z[1] / no;
    Z[2] = Z[2] / no;
    return;
}

/* Test at the start if the Hamiltonian and dipole files are sensible */
int control(t_non* non) {
    float *mu_eg, *Hamil_i_e;
    FILE *H_traj, *mu_traj;
    FILE *x_traj;
    int itime, N_samples;
    int samples;
    int nn2;

    nn2 = non->singles * (non->singles + 1) / 2;
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    mu_eg = (float *)calloc(non->singles, sizeof(float));
    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        return 1;
    }
    mu_traj = fopen(non->dipoleFName, "rb");
    if (mu_traj == NULL) {
        printf("Dipole file %s not found!\n", non->dipoleFName);
        return 1;
    }
    N_samples = (non->length - non->tmax1 - 1) / non->sample + 1;
    if (N_samples < 0) {
        printf("Insufficient data to calculate spectrum.\n");
        printf("Please, lower max times or provide longer\n");
        printf("trajectory.\n");
        return 1;
    }

    // Check first element
    // Read Hamiltonian
    if (read_He(non, Hamil_i_e, H_traj, 0) != 1) {
        printf("Failed initial control\n");
        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
        return 1;
    }
    if (read_mue(non, mu_eg, mu_traj, 0, 0) != 1) {
        printf("Failed initial control\n");
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n", 0, 0);
        return 1;
    }
    if (!strcmp(non->hamiltonian, "Coupling")) { }
    else if (!strcmp(non->hamiltonian, "TransitionDipole") || !strcmp(non->hamiltonian, "ExtendedDipole")) {
        x_traj = fopen(non->positionFName, "rb");
        if (x_traj == NULL) {
          printf("Position file not found!\n");
          return 1;
        }
        fclose(x_traj);
    }

    else {
        // Check last element
        if (read_mue(non, mu_eg, mu_traj, non->length - 1, 2) != 1) {
            printf("Dipole trajectory file to short, could not fill buffer!!!\n");
            printf("ITIME %d %d\n", non->length - 1, 2);
            return 1;
        }
        // Read Hamiltonian
        if (read_He(non, Hamil_i_e, H_traj, non->length - 1) != 1) {
            printf(RED "Failed initial control\n");
            printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
            printf("Real file length shorter than specified with Length keyword!\n" RESET);
            return 1;
        }
        // Check Hamiltonian elements
        if (Hamil_i_e[0] + non->shifte > non->max1 || Hamil_i_e[0] + non->shifte < non->min1) {
            printf(RED "Warning: Hamiltonian value %f outside expected range.\n", Hamil_i_e[0] + non->shifte);
            printf("Expected frequency range: %f to %f.\n", non->min1, non->max1);
            printf("Computation will continue, but check is the number above is realistic\n");
            printf("You may have specified wrong number of sites!\n" RESET);
            printf("---------------------------------------------------------------------\n");
        }
    }
    free(mu_eg);
    free(Hamil_i_e);
    fclose(mu_traj), fclose(H_traj);
    return 0;
}

/* Routine for autodetection of the number of singles */
int autodetect_singles(t_non* non){
    float *Hamil_i_e;
    FILE *H_traj;
    int i;
    float f;
    int samples;
    int n,nn2;
    int identified;
    int identified2;

    identified=0;
    identified2=0;
    nn2=non->singles*(non->singles+1)/2;
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        return 1;
    }

    /* Check if user provided Singles is correct */
    n=non->singles;
    if (!strcmp(non->hamiltonian, "Coupling")) {
       identified2=-2;
    }
    if (!strcmp(non->hamiltonian, "TransitionDipole") || !strcmp(non->hamiltonian, "ExtendedDipole")){
       identified2=-2;
    }
    if (!strcmp(non->hamiltonian, "Full")) {
         fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
         fread(&i,sizeof(int),1,H_traj);
         fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
         fread(&f,sizeof(float),1,H_traj);
         //printf("%d %d %f\n",n,i,f);
         if (abs(i)<1000 || f==floorf(f)){
            printf("Autodetected potential singles at %d\n",n);
            if (n==non->singles){
               identified2=-1;
            }
            if (n!=non->singles){
               identified2=n;
            }
         }
      }

    /* Check alternative settings */
    if (identified2==-1){
      identified=-1;
    } else {
      for (n=1;n<non->singles*10;n++){
        if (!strcmp(non->hamiltonian, "Coupling")) {
           identified=-2;
        }
        if (!strcmp(non->hamiltonian, "TransitionDipole") || !strcmp(non->hamiltonian, "ExtendedDipole")){
           identified=-2;
        }
        if (!strcmp(non->hamiltonian, "Full")) {
	  fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
	  fread(&i,sizeof(int),1,H_traj);
          fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
          fread(&f,sizeof(float),1,H_traj);
	  //printf("%d %d %f\n",n,i,f);
          if (abs(i)<1000 || f==floorf(f)){
            printf("Autodetected potential singles at %d\n",n);
	    if (n==non->singles){
	       identified=-1;
	       break;
	    }
	    if (n!=non->singles){
	       identified=n;
	       break;
	    }
	  }
        }
      }
    }

    if (identified==0){
       printf(RED "Warning: Autodetection of sites failed. Verify that your Singles setting is correct.\n" RESET);
    }
    if (identified==-1){
       printf("Singles confirmed by auto detection.\n");
    }
    if (identified>0){
      printf(RED "Warning: Singles keyword may be specified incorrectly!\n");
      printf("Autodetection suggested %d singles.\n" RESET,identified);
    }
    fclose(H_traj);
    free(Hamil_i_e);
    return 0;
}

/* Multiply with double exciton dipole mu_ef on single states */
void dipole_double(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }

    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[index] += over[i] * cr[i];
        fi[index] += over[i] * ci[i];
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[index] += dipole[i] * cr[j];
            fi[index] += dipole[i] * ci[j];
            fr[index] += dipole[j] * cr[i];
            fi[index] += dipole[j] * ci[i];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on ground states */
void dipole_double_ground(t_non* non, float* dipole, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }

    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[index] = over[i] ;
        //! no fi since double excitation from ground state,
        //no i,j loop due to double excitation on one state only
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on single states */
void dipole_double_ES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;

    for (i = 0; i < non->singles; i++) {
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[index] += dipole[i] * cr[j];
            fi[index] += dipole[i] * ci[j];
            fr[index] += dipole[j] * cr[i];
            fi[index] += dipole[j] * ci[i];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on double states */
void dipole_double_last(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < non->singles; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }
    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[i] += over[i] * cr[index];
        fi[i] += over[i] * ci[index];
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[j] += dipole[i] * cr[index];
            fi[j] += dipole[i] * ci[index];
            fr[i] += dipole[j] * cr[index];
            fi[i] += dipole[j] * ci[index];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on double states */
void dipole_double_last_ES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < non->singles; i++) fr[i] = 0, fi[i] = 0;
    for (i = 0; i < non->singles; i++) {
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[j] += dipole[i] * cr[index];
            fi[j] += dipole[i] * ci[index];
            fr[i] += dipole[j] * cr[index];
            fi[i] += dipole[j] * ci[index];
        }
    }
    return;
}

/* Return the distance between two locations squared */
float distance(float *rf,float *ri,int a,int b,int N,float box){
  float d,r;
  int x;
  d=0;
  for (x=0;x<3;x++){
    r=rf[3*a+x]-ri[3*b+x];
    if (r>box/2) r=r-box;
    if (r<-box/2) r=r+box;
    d+=r*r;
  }
  return d;
}

/* Return the distance between two locations along a direction x */
float distance_x(float *rf,float *ri,int a,int b,int N,float box,int x){
  float r;
  r=rf[3*a+x]-ri[3*b+x];
  if (r>box/2) r=r-box;
  if (r<-box/2) r=r+box;
  return r;
}

/* Return the distance between two locations squared */
float distance3(float *rf,float *ri,int a,int b,int N,float *box){
  float d,r;
  int x;
  d=0;
  for (x=0;x<3;x++){
    r=rf[3*a+x]-ri[3*b+x];
    if (r>box[x]/2) r=r-box[x];
    if (r<-box[x]/2) r=r+box[x];
    d+=r*r;
  }
  return d;
}

/* Return the distance between two locations along a direction x */
float distance3_x(float *rf,float *ri,int a,int b,int N,float *box,int x){
  float r;
  r=rf[3*a+x]-ri[3*b+x];
  if (r>box[x]/2) r=r-box[x];
  if (r<-box[x]/2) r=r+box[x];
  return r;
}

float pbc1(float r, int x, float *box){
  // Correct for pbc if active
  if (box[0]>0.0){
     if (r>box[x]/2) r=r-box[x];
     if (r<-box[x]/2) r=r+box[x];
  }
  return r;
}

// Diagonalize real nonsymmetric matrix. Output complex eigenvalues, left and right eigenvectors.
void diagonalize_real_nonsym(float* K, float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N) {
    int INFO, lwork;
    float *work, *Kcopy;
    int i, j;
    int *pivot;
    int M;
    /* Diagonalization*/
    /* Find lwork for diagonalization */
    lwork = -1;
    work = (float *)calloc(1, sizeof(float));
    sgeev_("V", "V", &N, Kcopy, &N, eig_re, eig_im, evecL, &N, evecR, &N, work, &lwork, &INFO);
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    Kcopy = (float *)calloc(N * N, sizeof(float));
    /* Copy matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Kcopy[i * N + j] = K[i * N + j];
        }
    }

    /* Do diagonalization*/
    sgeev_("V", "V", &N, Kcopy, &N, eig_re, eig_im, evecL, &N, evecR, &N, work, &lwork, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to diagonalize a matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    free(work);

    /* Copy matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            ivecL[i * N + j] = evecL[i * N + j];
            ivecR[i * N + j] = evecR[i * N + j];
        }
    }
    return;
}


/* Integrate the rate response */
void integrate_rate_response(float *rate_response,int T,float *is13,float *isimple){
    int i;
    float simple; /* Variable for trapezium integral */
    float simp13; /* Variable for Simpsons 1/3 rule integral */
    simple=0;
    simp13=0;
    for (i=0;i<T;i++){
        if (i==0){
	    simple+=rate_response[i]/2;
	    simp13+=rate_response[i]/3;
	} else if (i%2==0){
	    simple+=rate_response[i];
            simp13+=2*rate_response[i]/3;
        } else {
	     simple+=rate_response[i];
            simp13+=4*rate_response[i]/3;
        }
    }

    /* Check for difference between integration methods */
    if (fabs(simple-simp13)/fabs(simple)>0.05){
        printf("\n");
        printf(YELLOW "Warning the timesteps may be to large for integration!\n" RESET);
        printf(YELLOW "Simple integral value: %f\n Simpson 1/3: %f\n",simple,simp13);
        printf(YELLOW "This difference is larger than 5%%.\n" RESET);
	printf(YELLOW "The trapezium rule value is used.\n\n" RESET);
    }

    /* Check for difference between initial and final value */
    if (fabs(rate_response[T-1])*50>rate_response[0]){
	    printf("\n");
            printf(YELLOW "Final value of rate response is %f %%\n",fabs(rate_response[T-1])*100/fabs(rate_response[0]));
	    printf("of the initial value. You may avearge over too\n");
	    printf("few samples (decrease the value of Samplerate) or\n");
	    printf("your chosen coherence time of %d steps, may\n",T);
	    printf("be too short for the coherence to decay.\n." RESET);
	    printf("\n");
    }

    /* Store results in variables for return */
    *isimple=simple;
    *is13=simp13;
}

/* Write a square matrix to a text file */
void write_matrix_to_file(char fname[],float *matrix,int N){
    FILE *file_handle;
    int i,j;
        file_handle=fopen(fname,"w");
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            fprintf(file_handle,"%10.14e ",matrix[i*N+j]);
        }
        fprintf(file_handle,"\n");
    }
    fclose(file_handle);
}

/* Read a square matrix from a text file */
void read_matrix_from_file(char fname[],float *matrix,int N){
    FILE *file_handle;
    int i,j;
    file_handle=fopen(fname,"r");
    if (file_handle == NULL) {
        printf("Error opening the file %s.\n",fname);
        exit(0);
    }
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            fscanf(file_handle,"%f",&matrix[i*N+j]);
        }
    }
    fclose(file_handle);
}

/* Read a vector from a text file */
void read_vector_from_file(char fname[],float *vector,int N){
    FILE *file_handle;
    int i;
    file_handle=fopen(fname,"r");
    if (file_handle == NULL) {
        printf("Error opening the file %s.\n",fname);
        exit(0);
    }
    for (i=0;i<N;i++){
        fscanf(file_handle,"%f",&vector[i]);
    }
    fclose(file_handle);
}

// Replace the extension of a filename if it matches the old extension, otherwise return NULL
char *replace_ext(const char *filename,
                  const char *old_ext,
                  const char *new_ext)
{
    size_t len = strlen(filename);
    size_t old_len = strlen(old_ext);
    size_t new_len = strlen(new_ext);

    if (len < old_len) return NULL;

    // check if filename ends with old_ext
    if (strcmp(filename + len - old_len, old_ext) != 0)
        return NULL;

    // allocate new string
    size_t new_size = len - old_len + new_len + 1;
    char *result = malloc(new_size);
    if (!result) return NULL;

    // copy base part
    memcpy(result, filename, len - old_len);

    // append new extension
    memcpy(result + (len - old_len), new_ext, new_len);

    // null terminate
    result[new_size - 1] = '\0';

    return result;
}

// Compare two strings ignoring case
int strcmp_nocase(const char *s1, const char *s2)
{
    while (*s1 && *s2) {
        unsigned char c1 = (unsigned char)*s1;
        unsigned char c2 = (unsigned char)*s2;

        c1 = (unsigned char)tolower(c1);
        c2 = (unsigned char)tolower(c2);

        if (c1 != c2)
            return c1 - c2;

        s1++;
        s2++;
    }

    return (unsigned char)tolower((unsigned char)*s1)
         - (unsigned char)tolower((unsigned char)*s2);
}

/* Save linear response function to file */
int save_time_domain_response(t_non *non,const char *filename,float *re_S_1,float *im_S_1,int pro_dim,int samples){
    int ip, t1;
    float time,re,im;
    FILE *outone;
    if (strcmp_nocase(non->outputformat, "Normal")==0) {
        outone=fopen(filename,"w");
        for (t1=0;t1<non->tmax1;t1+=non->dt1){
            fprintf(outone,"%f ",t1*non->deltat);
            for (ip=0;ip<pro_dim;ip++){
                fprintf(outone,"%e %e ",re_S_1[t1+ip*non->tmax]/samples,im_S_1[t1+ip*non->tmax]/samples);
            }
            fprintf(outone,"\n");
        }
        fclose(outone);
        return 0;
    } else if (strcmp_nocase(non->outputformat, "Binary")==0) {
        char *binary_fname = replace_ext(filename, ".dat", ".bin");
        outone=fopen(binary_fname,"wb");
        for (t1=0;t1<non->tmax1;t1+=non->dt1){
            time=t1*non->deltat;
            fwrite(&time,sizeof(float),1,outone);
            for (ip=0;ip<pro_dim;ip++){
                re=re_S_1[t1+ip*non->tmax]/samples;
                im=im_S_1[t1+ip*non->tmax]/samples;
                fwrite(&re,sizeof(float),1,outone);
                fwrite(&im,sizeof(float),1,outone);
            }
        }
        fclose(outone);
        return 0;
    } else {
        printf("\n");
        printf("To store response function use OutputFormat Normal or Binary.\n\n");
    }
    return 0;
}