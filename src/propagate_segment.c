#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <cblas.h>
#include "types.h"
#include "NISE_subs.h"
#include "randomlib.h"
#include "util/asprintf.h"
#include "propagate_segment.h"

/* Standard propagation of a single vector */
/* display is t1*x for displaying info at first step, that is when t1 and x are both zero */
/* and we have the first sample */
void propagate_vector_segments(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display, int N_i){
   int elements;
   if (non->propagation==1) propagate_vec_coupling_S_segments(non,Hamil_i_e,vecr,veci,non->ts,sign, N_i);
//    if (non->propagation==3) propagate_vec_RK4(non,Hamil_i_e,vecr,veci,non->ts,sign);
//    if (non->propagation==0){
//       if (non->thres==0 || non->thres>1){
//          propagate_vec_DIA(non,Hamil_i_e,vecr,veci,sign);
//       } else {
//          elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr,veci,sign);
//          if (samples==non->begin){
//              if (display==0){
//                  printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
//                  printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
//                  printf("Suggested truncation %f.\n",0.001);
//              }
//          }
//       }
//    }
   return;
}

/* Standard propagation of a collection of N vectors */
void propagate_matrix_segments(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display, int N_i){
   int elements;
   int N,j;
   N = N_i;
   // N=non->singles;
#pragma omp parallel for
   for (j=0;j<N_i;j++){
        propagate_vector_segments(non,Hamil_i_e,vecr+j*N,veci+j*N,sign,samples,display*j, N_i);
   }
   return;
}

/* Propagate using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S_segments(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, int sign, int N_i) {
    float f;
    int index, N;
    float *H1, *H0, *re_U, *im_U;
    int *col, *row;
    float *ocr, *oci;
    int a, b, c;
    float J;
    float cr1, cr2, ci1, ci2;
    float co, si;
    int i, k, kmax;
    int N2;

    // N = non->singles;
    N = N_i;
    N2=(N*(N-1))/2;
    f = non->deltat * icm2ifs * twoPi * sign / m;
    H0 = (float *)malloc(N*sizeof(float));
    H1 = (float *)malloc(N2*sizeof(float));
    col = (int *)malloc(N2*sizeof(int));
    row = (int *)malloc(N2*sizeof(int));
    re_U = (float *)malloc(N*sizeof(float));
    im_U = (float *)malloc(N*sizeof(float));
    ocr = (float *)malloc(N*sizeof(float));
    oci = (float *)malloc(N*sizeof(float));


    /* Build Hamiltonians H0 (diagonal) and H1 (coupling) */
    k = 0;
    for (a = 0; a < N; a++) {
        H0[a] = Hamiltonian_i[Sindex(a, a, N)]; /* Diagonal */
        for (b = a + 1; b < N; b++) {
            index = Sindex(a, b, N);
            if (fabs(Hamiltonian_i[index]) > non->couplingcut) {
                H1[k] = Hamiltonian_i[index];
                col[k] = a, row[k] = b;
                k++;
            }
        }
    }
    kmax = k;

    /* Exponentiate diagonal [U=exp(-i/2h H0 dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(0.5 * H0[a] * f);
        im_U[a] = -sin(0.5 * H0[a] * f);
    }

    for (i = 0; i < m; i++) {
        /* Multiply on vector first time */
        for (a = 0; a < N; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = ci[a] * re_U[a] + cr[a] * im_U[a];
        }

        /* Account for couplings */
        for (k = 0; k < kmax; k++) {
            a = col[k];
            b = row[k];
            J = H1[k];
            J = J * f;
            si = -sin(J);
            co = sqrt(1 - si * si);
            cr1 = co * ocr[a] - si * oci[b];
            ci1 = co * oci[a] + si * ocr[b];
            cr2 = co * ocr[b] - si * oci[a];
            ci2 = co * oci[b] + si * ocr[a];
            ocr[a] = cr1, oci[a] = ci1, ocr[b] = cr2, oci[b] = ci2;
        }

        /* Multiply on vector second time */
        for (a = 0; a < N; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = oci[a] * re_U[a] + ocr[a] * im_U[a];
        }
    }

    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}

/* Propagate using a precalculated snapshot */
void propagate_snapshot(float *U_snap_re, float *U_snap_im, float **U_comp_re, float **U_comp_im, float **temp_re, float **temp_im, int N){
    // multiply the propagator with the next snapshot
    // requires precalculated snapshots
    // overwrites the original propagator
    int i,j;
    complex_matrix_product(U_snap_re, U_snap_im, *U_comp_re, *U_comp_im, *temp_re, *temp_im, N,N,N);

    // swap pointers of the temp and comp matrices
    // this ensures proper pointing to the updated propagator after each iteration
    float *tmp;
    tmp = *U_comp_re; *U_comp_re = *temp_re; *temp_re = tmp;
    tmp = *U_comp_im; *U_comp_im = *temp_im; *temp_im = tmp;
}