#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif
//#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "calc_CG_2DES.h"
#include "calc_FD_CG_2DES.h"
#include <stdarg.h>
#include "project.h"
#include "propagate.h"
#include "read_trajectory.h"
#include "eq_den.h"
//#include "mpi.h"
#include "MPI_subs.h"
#include <complex.h>

void call_final_FD_CG_2DES(
  t_non *non,float *P_DA,int pro_dim,float *re_doorway,float *im_doorway,
  float *re_window_SE, float *im_window_SE,float *re_window_GB, float *im_window_GB,
  float *re_window_EA, float *im_window_EA,float *re_2DES , float *im_2DES){
}

void read_in_QF_for_FD_CG_2DES(t_non *non,double *Q1, double *Q2){
    /* Define variables for QFactors */
    int segments;
    int f_states;
    int e_states;
    int g_states;

    segments = project_dim(non);
    f_states = segments*(segments+1)/2;
    e_states = 2*segments;
    g_states = 3;

    /* Open .dat files */
    FILE *Q1sFile = fopen("Q1s.dat","r");
    if (!Q1sFile) {
        fprintf(stderr, "File not found: Q1s.dat\n");
        return;
    }
    FILE *Q2sFile = fopen("Q2s.dat","r");
    if (!Q2sFile) {
        fprintf(stderr, "File not found: Q2s.dat\n");
        return;
    }

    /* Read-in and store values for Q1s and Q2s */
    for (int i = 0; i < segments; ++i) {
        fscanf(Q1sFile, "%lf", &Q1[i]);
    }
    fclose(Q1sFile);

    for (int i = 0; i < f_states; ++i) {
        fscanf(Q2sFile, "%lf", &Q2[i]);
    }
    fclose(Q2sFile);
    
    /* Display the values for demonstration purposes.
    Can modify printLevel in input file : PrintLevel 2 */
    if(non->printLevel>1){
        printf("Q1 values:\n");
        for (int i = 0; i < segments; ++i) {
            printf("Q1_%d = %f\n", i + 1, Q1[i]);
        }
            
        printf("Q2 values:\n");
        for (int i = 0; i < f_states; ++i) {
            printf("Q2_%d = %f\n", i + 1, Q2[i]);
        }
    }
}

void calc_FD_CG_2DES(t_non *non){
    int segments;
    int f_states;
    int e_states;
    int g_states;

    segments = project_dim(non);
    f_states = segments*(segments+1)/2;
    e_states = 2*segments;

    double Q1[segments];
    double Q2[f_states];
    read_in_QF_for_FD_CG_2DES(non,Q1,Q2);
}

void FD_CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
    float *re_window_SE,float *im_window_SE,float *re_window_GB, float *im_window_GB,
    float *re_window_EA,float *im_window_EA,float *P_DA,int N, char *waittime,int wfile){
}