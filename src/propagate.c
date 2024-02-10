#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "randomlib.h"
#include "util/asprintf.h"
#include "propagate.h"

/* Easy call standard propagation schemes with control info */



/* Standard propagation of a single vector */
/* display is t1*x for displaying info at first step, that is when t1 and x are both zero */
/* and we have the first sample */
void propagate_vector(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display){
   int elements;
   if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr,veci,non->ts,sign);
   if (non->propagation==3) propagate_vec_RK4(non,Hamil_i_e,vecr,veci,non->ts,sign);
   if (non->propagation==0){
      if (non->thres==0 || non->thres>1){
         propagate_vec_DIA(non,Hamil_i_e,vecr,veci,sign);
      } else {
         elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr,veci,sign);
         if (samples==non->begin){
             if (display==0){
                 printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
                 printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
                 printf("Suggested truncation %f.\n",0.001);
             }
         }
      }
   }
   return;
}


/* Standard propagation of a collection of N vectors */
void propagate_matrix(t_non *non,float * Hamil_i_e,float *vecr,float *veci,int sign,int samples,int display){
   int elements;
   int N,j;
   N=non->singles;
#pragma omp parallel for
   for (j=0;j<non->singles;j++){
        propagate_vector(non,Hamil_i_e,vecr+j*N,veci+j*N,sign,samples,display*j);
   }
   return;
}
        
      


/* Core functions for different propagation schemes are given below */

/* Propagate using standard matrix exponential */
void propagate_vec_DIA(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int sign) {
    float f;
    int index, N;
    float *H, *re_U, *im_U, *e;
    float *cnr, *cni;
    float *crr, *cri;
    float re, im;
    int a, b, c;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N * N, sizeof(float));
    cni = (float *)calloc(N * N, sizeof(float));
    crr = (float *)calloc(N * N, sizeof(float));
    cri = (float *)calloc(N * N, sizeof(float));
    /* Build Hamiltonian */
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; /* Diagonal */
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
    /* Exponentiate [U=exp(-i/h H dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    /* Transform to site basis */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] += H[b + a * N] * re_U[b], cni[b + a * N] += H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                crr[a + c * N] += H[b + a * N] * cnr[b + c * N], cri[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    /* The one exciton propagator has been calculated */

    for (a = 0; a < N; a++) {
        cnr[a] = 0, cni[a] = 0;
        for (b = 0; b < N; b++) {
            cnr[a] += crr[a + b * N] * cr[b] - cri[a + b * N] * ci[b];
            cni[a] += crr[a + b * N] * ci[b] + cri[a + b * N] * cr[b];
        }
    }

    for (a = 0; a < N; a++) {
        cr[a] = cnr[a], ci[a] = cni[a];
    }


    free(cnr), free(cni), free(re_U), free(im_U), free(H), free(e);
    free(crr), free(cri);
    return;
}


/* Do the propagation for t2 using the matrix exponent and using a single diagonalization */
/* for all vectors */
void propagate_t2_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,float **vr,float **vi,int sign){
    float f;
    int index, N;
    float *H, *re_U, *im_U, *e;
    float re, im;
    int a, b, c;
    int t1;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));

    /* Build Hamiltonian */
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; /* Diagonal*/
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
    /* Exponentiate [U=exp(-i/h H dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    /* Do the single t1 independent vector */
    /* Transfer to eigen basis */
    matrix_on_vector(H,cr,ci,N);
    /* Multiply with matrix exponent */
    vector_on_vector(re_U,im_U,cr,ci,N);
    /* Transfer back to site basis */
    trans_matrix_on_vector(H,cr,ci,N);

    /* Do all the t1 dependent vectors */
#pragma omp parallel for shared(non,re_U,im_U,H,vr,vi) schedule(static,1)
    for (t1=0;t1<non->tmax1;t1++){
        /* Transfer to eigen basis */
        matrix_on_vector(H,vr[t1],vi[t1],N);
        /* Multiply with matrix exponent */
        vector_on_vector(re_U,im_U,vr[t1],vi[t1],N);
        /* Transfer back to site basis */
        trans_matrix_on_vector(H,vr[t1],vi[t1],N);
    }
    free(re_U), free(im_U), free(H), free(e);
    return;
}



/* Propagate using matrix exponential sparce */
int propagate_vec_DIA_S(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int sign) {
    int elements;
    float f;
    int index, N, N2;
    float *H, *re_U, *im_U, *e;
    float *cnr, *cni;
    float *crr, *cri;
    float re, im;
    int a, b, c;

    N = non->singles;
    N2 = N * N;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N2, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N2, sizeof(float));
    cni = (float *)calloc(N2, sizeof(float));
    crr = (float *)calloc(N2, sizeof(float));
    cri = (float *)calloc(N2, sizeof(float));

    /* Build Hamiltonian */
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; /* Diagonal*/
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }
    diagonalizeLPD(H, e, N);
    /* Exponentiate [U=exp(-i/h H dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    /* Transform to site basis */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] = H[b + a * N] * re_U[b], cni[b + a * N] = H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (c = 0; c < N; c++) {
            crr[a + c * N] = 0, cri[a + c * N] = 0;
            for (b = 0; b < N; b++) {
                crr[a + c * N] += H[b + a * N] * cnr[b + c * N], cri[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    /* The one exciton propagator has been calculated */

    elements = 0;
    for (a = 0; a < N; a++) {
        cnr[a] = 0, cni[a] = 0;
        for (b = 0; b < N; b++) {
            if ((crr[a + b * N] * crr[a + b * N] + cri[a + b * N] * cri[a + b * N]) > non->thres) {
                elements++;
                cnr[a] += crr[a + b * N] * cr[b] - cri[a + b * N] * ci[b];
                cni[a] += crr[a + b * N] * ci[b] + cri[a + b * N] * cr[b];
            }
        }
    }

    for (a = 0; a < N; a++) {
        cr[a] = cnr[a], ci[a] = cni[a];
    }

    free(crr), free(cri);
    free(cnr), free(cni), free(re_U), free(im_U), free(H), free(e);

    return elements;
}

/* Propagate singles using the Runge Kutta 4 algorithm */
void propagate_vec_RK4(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign){
    float f;
    int index, N;
    float *H0;
    float *k1r,*k2r,*k3r,*k4r;
    float *k1i,*k2i,*k3i,*k4i;
    int *col, *row;
    float *ocr, *oci;
    int a, b, c;
    float J;
    float cr1, cr2, ci1, ci2;
    int i, k, kmax;
    int N2;

    /* printf("Entered the RK4 routine.\n"); */
    
    N = non->singles;
    N2=(N*(N+1))/2;
    f = non->deltat * icm2ifs * twoPi * sign / m;
    H0 = (float *)malloc(N2*sizeof(float));
    col = (int *)malloc(N2*sizeof(int));
    row = (int *)malloc(N2*sizeof(int));
    k1r = (float *)calloc(N,sizeof(float));
    k1i = (float *)calloc(N,sizeof(float));
    k2r = (float *)calloc(N,sizeof(float));
    k2i = (float *)calloc(N,sizeof(float));
    k3r = (float *)calloc(N,sizeof(float));
    k3i = (float *)calloc(N,sizeof(float));
    k4r = (float *)calloc(N,sizeof(float));
    k4i = (float *)calloc(N,sizeof(float));

    /* Build sparse Hamiltonians H0 */
    k = 0;
    for (a = 0; a < N; a++) {
        for (b = a; b < N; b++) {
            index = Sindex(a, b, N);
            if (fabs(Hamiltonian_i[index]) > non->couplingcut || a==b) {
                index = Sindex(a, b, N);
                H0[k] = f* Hamiltonian_i[index];
                col[k] = a, row[k] = b;
                k++;
            }
        }
    }
    kmax = k;

    /* We assume the Hamiltonian to be time independent! */

    /* Multi-step loop */
    for (i=0;i<m;i++){
	if (i>0) {
            clearvec(k1r,N);
            clearvec(k1i,N);
            clearvec(k2r,N);
            clearvec(k2i,N);
            clearvec(k3r,N);
            clearvec(k3i,N); 
            clearvec(k4r,N);
            clearvec(k4i,N);	    
	}
        /* Find k1 */
        for (k=0;k<kmax;k++){
            k1r[col[k]]+=H0[k]*ci[row[k]];
            k1i[col[k]]-=H0[k]*cr[row[k]];
            if (row[k]!=col[k]){
                k1r[row[k]]+=H0[k]*ci[col[k]];
                k1i[row[k]]-=H0[k]*cr[col[k]];
            }
        }
        /* Find k2 */
        for (k=0;k<kmax;k++){
            k2r[col[k]]+=H0[k]*(ci[row[k]]+k1i[row[k]]*0.5);
            k2i[col[k]]-=H0[k]*(cr[row[k]]+k1r[row[k]]*0.5);
            if (row[k]!=col[k]){
                k2r[row[k]]+=H0[k]*(ci[col[k]]+k1i[col[k]]*0.5);
                k2i[row[k]]-=H0[k]*(cr[col[k]]+k1r[col[k]]*0.5);
            }
        }
        /* Find k3 */
        for (k=0;k<kmax;k++){
            k3r[col[k]]+=H0[k]*(ci[row[k]]+k2i[row[k]]*0.5);
            k3i[col[k]]-=H0[k]*(cr[row[k]]+k2r[row[k]]*0.5);
            if (row[k]!=col[k]){
                k3r[row[k]]+=H0[k]*(ci[col[k]]+k2i[col[k]]*0.5);
                k3i[row[k]]-=H0[k]*(cr[col[k]]+k2r[col[k]]*0.5);
            }
        }
        /* Find k4 */
        for (k=0;k<kmax;k++){
            k4r[col[k]]+=H0[k]*(ci[row[k]]+k3i[row[k]]);
            k4i[col[k]]-=H0[k]*(cr[row[k]]+k3r[row[k]]);
            if (row[k]!=col[k]){
                k4r[row[k]]+=H0[k]*(ci[col[k]]+k3i[col[k]]);
                k4i[row[k]]-=H0[k]*(cr[col[k]]+k3r[col[k]]);
            }
        }

        /* Update wavefunction */
        for (k=0;k<N;k++){
            cr[k]=cr[k]+(k1r[k]+2*k2r[k]+2*k3r[k]+k4r[k])/6.0;
            ci[k]=ci[k]+(k1i[k]+2*k2i[k]+2*k3i[k]+k4i[k])/6.0;
        }
    }

    free(H0),free(col),free(row);
    free(k1r),free(k2r),free(k3r),free(k4r);
    free(k1i),free(k2i),free(k3i),free(k4i);
    return;

}

/* Propagate singles using the Runge Kutta 4 algorithm */
void propagate_vec_RK4_doubles(t_non *non,float *Hamiltonian_i,float *cr,float *ci,int m,int sign,float *Anh){
    float f;
    int index, N;
    float *H0;
    float *k1r,*k2r,*k3r,*k4r;
    float *k1i,*k2i,*k3i,*k4i;
    int *col, *row;
    float *ocr, *oci;
    int a, b, c;
    float J;
    float cr1, cr2, ci1, ci2;
    int i, k, kmax;
    int N2;

    /* printf("Entered the RK4 routine.\n"); */

    N = non->singles;
    N2=(N*(N+1))/2;
    f = non->deltat * icm2ifs * twoPi * sign / m;
    H0 = (float *)malloc(N2*sizeof(float));
    col = (int *)malloc(N2*sizeof(int));
    row = (int *)malloc(N2*sizeof(int));
    k1r = (float *)calloc(N,sizeof(float));
    k1i = (float *)calloc(N,sizeof(float));
    k2r = (float *)calloc(N,sizeof(float));
    k2i = (float *)calloc(N,sizeof(float));
    k3r = (float *)calloc(N,sizeof(float));
    k3i = (float *)calloc(N,sizeof(float));
    k4r = (float *)calloc(N,sizeof(float));
    k4i = (float *)calloc(N,sizeof(float));

    /* Build sparse Hamiltonians H0 */
    k = 0;
    for (a = 0; a < N; a++) {
        for (b = a; b < N; b++) {
            index = Sindex(a, b, N);
            if (fabs(Hamiltonian_i[index]) > non->couplingcut || a==b) {
                index = Sindex(a, b, N);
                H0[k] = f* Hamiltonian_i[index];
                col[k] = a, row[k] = b;
                k++;
            }
        }
    }
    kmax = k;


    printf("Not implemented yet!\n");
    exit(0);
}


/* Propagate using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, int sign) {
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

    N = non->singles;
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

/* Propagate threel level doubles using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S_doubles(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, float* Anh) {
    int N = non->singles;
    int N2 = N * (N + 1) / 2;
    const float f = non->deltat * icm2ifs * twoPi / m;
    float si, co;
    float si2, co2;
    float* H0 = malloc(N2* sizeof(float));
    float* H1 = malloc(N * N / 2* sizeof(float));
    int* col = malloc(N * N / 2* sizeof(int));
    int* row = malloc(N * N / 2* sizeof(int));
    float* re_U = malloc(N2* sizeof(float));
    float* im_U = malloc(N2* sizeof(float));
    float* ocr = malloc(N2* sizeof(float));
    float* oci = malloc(N2* sizeof(float));

    /* Build Hamiltonians H0 (diagonal) and H1 (coupling) */
    for (int a = 0; a < N; a++) {
        const int indexa = Sindex(a, a, N);
        for (int b = a; b < N; b++) {
            int index = Sindex(a, b, N);
            /* Diagonal */
            H0[index] = Hamiltonian_i[indexa] + Hamiltonian_i[Sindex(b, b, N)];
            if (a == b) {
                if (non->anharmonicity == 0) {
                    H0[index] -= Anh[a];
                }
                else {
                    H0[index] -= non->anharmonicity;
                }
            }
        }
    }

    /* Build Hamiltonian H1 (coupling) */
    int kmax = 0;
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            int index = b + a * ((N << 1) - a - 1) / 2;
            /* Part of Sindex, but b > a is always true here */

            if (fabsf(Hamiltonian_i[index]) > non->couplingcut) {
                H1[kmax] = Hamiltonian_i[index];
                col[kmax] = a, row[kmax] = b;
                kmax++;
            }
        }
    }

    /* Exponentiate diagonal [U=exp(-i/2h H0 dt)] */
    for (int a = 0; a < N2; a++) {
        re_U[a] = cosf(0.5f * H0[a] * f);
        im_U[a] = -sinf(0.5f * H0[a] * f);
    }

    for (int i = 0; i < m; i++) {

        /* Multiply on vector first time */
        for (int a = 0; a < N2; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = cr[a] * im_U[a] + ci[a] * re_U[a];
        }

        /* Account for couplings */
        /* Loop over couplings */
        for (int k = 0; k < kmax; k++) {
            int a = col[k];
            int b = row[k];
            int c = 0;
            float J = H1[k] * f;
            int aNsum = a * ((N << 1) - a - 1) / 2;  // copied part of Sindex
            int bNsum = b * ((N << 1) - b - 1) / 2;  // copied part of Sindex

            /* Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> */

            /* c < a,b */
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < a; c++) {
                int index1 = a + c * ((N << 1) - c - 1) / 2;
                /* part of Sindex, but always a > c */
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                /* part of Sindex, but always b > c */
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* c == a */
            si2 = -sinf(J * sqrt2);
            co2 = sqrtf(1 - si2 * si2);
            for (; c == a; c++) {
                /* yeah, one iteration. but loop for consistency across all 5 */
                int index1 = a + c * ((N << 1) - c - 1) / 2;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co2 * ocr[index1] - si2 * oci[index2];
                float ci1 = co2 * oci[index1] + si2 * ocr[index2];
                float cr2 = co2 * ocr[index2] - si2 * oci[index1];
                float ci2 = co2 * oci[index2] + si2 * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* a < c < b */
            for (; c < b; c++) { // could be 0 iterations if (b == a + 1)
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* c == b */
            for (; c == b; c++) {
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co2 * ocr[index1] - si2 * oci[index2];
                float ci1 = co2 * oci[index1] + si2 * ocr[index2];
                float cr2 = co2 * ocr[index2] - si2 * oci[index1];
                float ci2 = co2 * oci[index2] + si2 * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* c > a,b */
            for (; c < N; c++) {
                int index1 = c + aNsum;
                int index2 = c + bNsum;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

        }

        /* Multiply on vector second time */
        for (int a = 0; a < N2; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = ocr[a] * im_U[a] + oci[a] * re_U[a];
        }
    }
    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}

/* Propagate two level doubles using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S_doubles_ES(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m) {
    int N = non->singles;
    int N2 = N * (N + 1) / 2;
    int N2old= N * ( N + 1 ) / 2;
    const float f = non->deltat * icm2ifs * twoPi / m;
    float* H0 = calloc(N2, sizeof(float));
    float* H1 = calloc(N * N / 2, sizeof(float));
    int* col = calloc(N * N / 2, sizeof(int));
    int* row = calloc(N * N / 2, sizeof(int));
    float* re_U = calloc(N2, sizeof(float));
    float* im_U = calloc(N2, sizeof(float));
    float* ocr = calloc(N2, sizeof(float));
    float* oci = calloc(N2, sizeof(float));
    float co,si;

    /* Build Hamiltonians H0 (diagonal) and H1 (coupling) */
    for (int a = 0; a < N; a++) {
        const int indexa = Sindex(a, a, N);
        for (int b = a; b < N; b++) {
            int index = Sindex(a, b, N);
            /* Diagonal */
            H0[index] = Hamiltonian_i[indexa] + Hamiltonian_i[Sindex(b, b, N)];
        }
    }

    /* Build Hamiltonian H1 (coupling) */
    int kmax = 0;
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            int index = b + a * ((N << 1) - a - 1) / 2;
            /* Part of Sindex, but b > a is always true here */

            if (fabsf(Hamiltonian_i[index]) > non->couplingcut) {
                H1[kmax] = Hamiltonian_i[index];
                col[kmax] = a, row[kmax] = b;
                kmax++;
            }
        }
    }

    /* Exponentiate diagonal [U=exp(-i/2h H0 dt)] */
    for (int a = 0; a < N2; a++) {
        re_U[a] = cosf(0.5f * H0[a] * f);
        im_U[a] = -sinf(0.5f * H0[a] * f);
    }

    /* Loop over Trotter steps */
    for (int i = 0; i < m; i++) {

        /* Multiply on vector first time */
        for (int a = 0; a < N2; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = cr[a] * im_U[a] + ci[a] * re_U[a];
        }

        /* Account for couplings */
        /* Loop over couplings */
        for (int k = 0; k < kmax; k++) {
            int a = col[k];
            int b = row[k];
            float J = H1[k] * f;

            /* Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> */
            /* TODO speedup */
            int c=0;
            int aNsum = a * ((N << 1) - a - 1) / 2; /* copied part of Sindex */
            int bNsum = b * ((N << 1) - b - 1) / 2; /* copied part of Sindex */

            /* Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> */

            /* c < a,b */
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < a; c++) {
                int index1 = a + c * ((N << 1) - c - 1) / 2;
                /* part of Sindex, but always a > c */
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                /* part of Sindex, but always b > c */
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* c == a */
            for (; c == a; c++) {
                /* yeah, one iteration. but loop for consistency across all 5 */
            }

            /* a < c < b */
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < b; c++) { /* could be 0 iterations if (b == a + 1) */
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            /* c == b */
            for (; c == b; c++) {
            }

            /* c > a,b */
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < N; c++) {
                int index1 = c + aNsum;
                int index2 = c + bNsum;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }
        }

        /* Multiply on vector second time */
        for (int a = 0; a < N2; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = ocr[a] * im_U[a] + oci[a] * re_U[a];
        }
    }

    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}


/* Diagonalize with LAPACK (destructive version) */
void diagonalizeLPD(float* H, float* v, int N) {
    int INFO, lwork;
    float *work, *Hcopy;
    int i, j;
    /* Find lwork for diagonalization */
    lwork = -1;
    work = (float *)calloc(1, sizeof(float));
    ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    Hcopy = (float *)calloc(N * N, sizeof(float));
    /* Copy Hamiltonian */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Hcopy[i * N + j] = H[i * N + j];
        }
    }

    /* Call LAPACK routine */
    ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to diagonalize a matrix...\n");
        exit(0);
    }

    /* Move eigenvectors */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            /* Converting from FORTRAN format */
            H[i * N + j] = Hcopy[j * N + i];
        }
    }
    /* Free space */
    free(Hcopy), free(work);
    return;
}

/* Build square Hamiltonian for triangular and diagonalize */
void build_diag_H(float* Hamiltonian_i, float* H, float* e, int N) {
    int a, b, c;
    /* Build square Hamiltonian from triagonal matrix */
    for (a = 0; a < N; a++) {
        /* Fill diagonal elements */
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; 
        /* Fill couplings */
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
}

/* Create truncated time-evolution operator */
int time_evolution_mat(t_non* non, float* Hamiltonian_i, float* Ur, float* Ui, int* R, int* C, int m) {
    float f, g;
    float *H, *re_U, *im_U, *e;
    float *cr, *ci;
    float *cnr, *cni;
    int indexA, indexB, N;
    int a, b, c, d;
    int elements;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi / m;
    H = (float *)calloc(N * N, sizeof(float));
    cr = (float *)calloc(N * N, sizeof(float));
    ci = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N * N, sizeof(float));
    cni = (float *)calloc(N * N, sizeof(float));
    /* Build Hamiltonian */
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }
    diagonalizeLPD(H, e, N);
    /* Exponentiate [U=exp(-i/h H dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }
    /* Transform to site basis */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] += H[b + a * N] * re_U[b], cni[b + a * N] += H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                cr[a + c * N] += H[b + a * N] * cnr[b + c * N], ci[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    /* The one exciton propagator has been calculated */
    elements = 0;
    /* Make sparse */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            if ((cr[a + b * N] * cr[a + b * N] + ci[a + b * N] * ci[a + b * N]) > non->thres) {
                Ur[elements] = cr[a + b * N], Ui[elements] = ci[a + b * N];
                R[elements] = a, C[elements] = b;
                elements++;
            }
        }
    }
    free(H), free(cr), free(ci), free(re_U), free(im_U), free(e), free(cnr), free(cni);
    return elements;
}

/* This function propagates the doubles for three level systems with the */
/* sparse time-evolution operator algorithm. The singles time evolution */
/* operator must be provided in Ur, Ui and fr, and fi contain the propagated *//* vectors */
void propagate_double_sparce(t_non* non, float* Ur, float* Ui, int* R, int* C, float* fr, float* fi, int elements,
                             int m, float* Anh) {
    int a, b, indexA, indexB;
    int N, Nf, s, t, u, v;
    float g, sqrt12;
    float *vr, *vi;
    float Uar, Ubr, Uai, Ubi;
    float *si, *co;
    float f;
    float norm1, norm2;
    float fm;
    int i;

    sqrt12 = 1.0 / sqrt2;
    f = non->deltat * icm2ifs * twoPi;
    fm = f * 0.5 / m;
    N = non->singles;
    Nf = non->singles * (non->singles + 1) / 2;
    vr = (float *)calloc(Nf, sizeof(float));
    vi = (float *)calloc(Nf, sizeof(float));
    co = (float *)calloc(N, sizeof(float));
    si = (float *)calloc(N, sizeof(float));

    if (non->anharmonicity != 0) {
        for (i = 0; i < N; i++) {
            co[i] = cos(fm * non->anharmonicity), si[i] = sin(fm * non->anharmonicity);
        }
    }
    else {
        for (i = 0; i < N; i++) {
            co[i] = cos(fm * Anh[i]), si[i] = sin(fm * Anh[i]);
        }
    }
    /* Repeat m times */
    for (i = 0; i < m; i++) {

        /* Anharmonicity */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = fr[indexA];
                vi[indexA] = fi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = co[a] * vr[indexA] - si[a] * vi[indexA];
            fi[indexA] = co[a] * vi[indexA] + si[a] * vr[indexA];
        }
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = 0;
                vi[indexA] = 0;
            }
        }

        for (a = 0; a < elements; a++) {
            s = R[a], t = C[a], Uar = Ur[a], Uai = Ui[a];
            for (b = 0; b <= a; b++) {
                u = R[b], v = C[b], Ubr = Ur[b], Ubi = Ui[b];
                g = 1;
                indexA = Sindex(s, u, N), indexB = Sindex(t, v, N);
                if (s != u) g = sqrt2;
                if (v != t) g = sqrt2;
                if (s != u && v != t) g = 1;
                vr[indexA] += (Uar * Ubr - Uai * Ubi) * fr[indexB] * g;
                vr[indexA] -= (Uai * Ubr + Uar * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uar * Ubr - Uai * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uai * Ubr + Uar * Ubi) * fr[indexB] * g;
            }
        }
        /* Anharmonicity */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                fr[indexA] = vr[indexA];
                fi[indexA] = vi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = co[a] * vr[indexA] - si[a] * vi[indexA];
            fi[indexA] = co[a] * vi[indexA] + si[a] * vr[indexA];
        }
    }
    free(si);
    free(co);
    free(vr);
    free(vi);
}


/* This function propagates the doubles for two level systems with the sparse */
/* time-evolution operator algorithm. The singles time evolution operator */
/* must be provided in Ur, Ui and fr, and fi contain the propagated vectors */
void propagate_double_sparce_ES(t_non* non, float* Ur, float* Ui, int* R, int* C, float* fr, float* fi, int elements,int m) {
    int a, b, indexA, indexB;
    int N, Nf, s, t, u, v;
    float g, sqrt12;
    float *vr, *vi;
    float Uar, Ubr, Uai, Ubi;
    float f;
    float norm1, norm2;
    float fm;
    int i;

    sqrt12 = 1.0 / sqrt2;
    f = non->deltat * icm2ifs * twoPi;
    fm = f * 0.5 / m;
    N = non->singles;
    Nf = non->singles * (non->singles + 1) / 2;
    vr = (float *)calloc(Nf, sizeof(float));
    vi = (float *)calloc(Nf, sizeof(float));

    /* Repeat m times */
    for (i = 0; i < m; i++) {

        /* Copy */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = fr[indexA];
                vi[indexA] = fi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = vr[indexA] ;
            fi[indexA] = vi[indexA] ;
        }
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = 0;
                vi[indexA] = 0;
            }
        }
        for (a = 0; a < elements; a++) {
            s = R[a], t = C[a], Uar = Ur[a], Uai = Ui[a];
            for (b = 0; b <= a; b++) {
                u = R[b], v = C[b], Ubr = Ur[b], Ubi = Ui[b];
                g = 1;
                indexA = Sindex(s, u, N), indexB = Sindex(t, v, N);
                if (s != u) g = 0; // No coupling to double excited states
                if (v != t) g = 0; 
                if (s != u && v != t) g = 1;
                vr[indexA] += (Uar * Ubr - Uai * Ubi) * fr[indexB] * g;
                vr[indexA] -= (Uai * Ubr + Uar * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uar * Ubr - Uai * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uai * Ubr + Uar * Ubi) * fr[indexB] * g;
            }   
                    }
        /* Copy */  
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                fr[indexA] = vr[indexA];
                fi[indexA] = vi[indexA];
            }   
        }   
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = vr[indexA] ;
            fi[indexA] = vi[indexA] ;
        }   
    }   
    free(vr);
    free(vi);
}


