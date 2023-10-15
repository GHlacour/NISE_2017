#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "read_trajectory.h"
#include "randomlib.h"
#include "util/asprintf.h"

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
int compare_string(char* string_to_compare, char* string_array[], int array_size) {
    for (int i = 0; i < array_size; i++) {
        if (!strcmp(string_to_compare, string_array[i])) {
            return i+1; // return the index of the matched string
        }
    }
    return 0; // if no match is found, return -1
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
    int s = difftime(t1, t0);
    int h = s / 3600;
    s = s % 3600;
    int m = s / 60;
    s = s % 60;

    char* text;
    asprintf(&text, "Time spent: %dh %dmin %ds\n", h, m, s);
    return text;
}

// Forms a string with the times for MPI_Wtime
char* MPI_time(double t0) {
    int ms =t0*1000; // Convert to milliseconds
    int h = ms / 3600000;
    ms = ms % 3600000;
    int m = ms / 60000;
    ms = ms % 60000;
    int s = ms / 1000;
    ms = ms % 1000;

    char* text;
    asprintf(&text, " %dh %dmin %ds %dms\n", h, m, s, ms);
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

    identified=0;
    nn2=non->singles*(non->singles+1)/2;
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        return 1;
    }
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
    if (identified==0){
       printf(RED "Warning: Autodetection of sites failed. You may need to increase Singles\n" RESET);
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
