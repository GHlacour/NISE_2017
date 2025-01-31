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

void calc_FD_CG_2DES(t_non *non){
    /* Define variables for QFactors */
    FILE *file = fopen("QFactors.txt","r");
    if (!file) {
        fprintf(stderr, "Dile not found: QFactors.txt\n");
        return;
    }

    double *Q1 = NULL, *Q2 = NULL;
    int Q1_size = 0, Q2_size = 0;
    char line[256];

    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Q1:") != NULL) {
            char *token = strtok(line + (strchr(line, '[') - line) + 1, " ");
            while (token != NULL) {
                Q1 = realloc(Q1, (Q1_size + 1) * sizeof(double));
                sscanf(token, "%lf", &Q1[Q1_size++]);
                token = strtok(NULL, " ");
            }
        } else if (strstr(line, "Q2:") != NULL) {
            char *token = strtok(line + (strchr(line, '[') - line) + 1, " ");
            while (token != NULL) {
                Q2 = realloc(Q2, (Q2_size + 1) * sizeof(double));
                sscanf(token, "%lf", &Q2[Q2_size++]);
                token = strtok(NULL, " ");
            }
        }
    }

    fclose(file);

    // Display the values for demonstration purposes
    printf("Q1 values:\n");
    for (int i = 0; i < Q1_size; ++i) {
        printf("Q1_%d = %f\n", i + 1, Q1[i]);
    }

    printf("Q2 values:\n");
    for (int i = 0; i < Q2_size; ++i) {
        printf("Q2_%d = %f\n", i + 1, Q2[i]);
    }

    free(Q1);
    free(Q2);
	printf("Hello World!\n");
}

void FD_CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
    float *re_window_SE,float *im_window_SE,float *re_window_GB, float *im_window_GB,
    float *re_window_EA,float *im_window_EA,float *P_DA,int N, char *waittime,int wfile){
}


