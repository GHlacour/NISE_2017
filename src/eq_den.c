#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "propagate.h"
#include "absorption.h"
//#include "luminescence.h"
#include "1DFFT.h"
#include "project.h"
#include "eq_den.h"
#include <complex.h>



/*Here we defined the equilibrium density operator*/
void eq_den(float *Hamiltonian_i, float *rho_l, int N, t_non *non){
 
  int index;
  float *H,*e,*e_1;
  float *rho_r;
  //float *rho_l;
  float *diag_sum;  /*Here is store the diagonal element for one segment*/
  //N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  e_1=(float *)calloc(N,sizeof(float));
  rho_r=(float *)calloc(N*N,sizeof(float));
  int a,b,c,pro_dim;
  int site_num,site_num_1,site_num_2,seg_num;/*site num is the number of the site, which used to make sure the index number later*/
  float kBT=non->temperature*0.6950859; // Here 0.6950859 is the K_b in the unit of cm-1
  float i_u = 0; 
  float u = 0; //initialized u
  /* projection */
  pro_dim=project_dim(non);
  diag_sum = (float *)calloc(pro_dim,sizeof(float));
  /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
  if (non->Npsites==non->singles){
    zero_coupling(Hamiltonian_i,non);
  } else {
    printf("Segment number and the projection number are different");
    exit(1);
  }

  // Build Hamiltonian, convert the triangular Hamiltonian to the  square Hamiltonian
  for (a=0;a<N;a++){
    /* Fill diagonal elements */
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; 
    /* Fill couplings */
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }

/*Here we diagonalize the square Hamiltonian, the H is replaced by the eigenvector in eigen basis, e is the eigenvalue */
  diagonalizeLPD(H,e,N);
  /* Exponentiate [U=exp(-H/kBT)] */
  for (a=0;a<N;a++){
      e_1[a]=exp(-e[a]/kBT);
  }
  

/*Here we need to normalize rho, exp(-KB*H)/[Trexp(-KB*H)] within every segment*/
/*We first sum the diagonal element in one segment*/
  for (seg_num=0;seg_num<pro_dim;seg_num++){
    for (site_num=0;site_num<non->singles;site_num++){
      if (seg_num == non->psites[site_num]){
        diag_sum[seg_num]+= e_1[site_num];
      }
    }
  }

    /*Secondly, we normalize the density matrix within one segment*/
    //for (seg_num=0;seg_num<pro_dim;seg_num++){
  for (site_num_1=0;site_num_1<non->singles;site_num_1++){
    seg_num=non->psites[site_num_1];
    e_1[site_num_1]=e_1[site_num_1]/diag_sum[seg_num];
  }
  /* Transform back to site basis */ 
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          rho_r[b+a*N]+=H[b+a*N]*e_1[a];
          rho_l[b+a*N]=0;
      }
  }  
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          for (c=0;c<N;c++){
            rho_l[a+c*N]+=H[b+a*N]*rho_r[b+c*N];
          }
      }
  }
  free(H);
  free(e_1);
  free(e);
  free(rho_r);
  return ;
}

// Multiply a real matrix on a real vector (vr,vi)
void matrix_on_real_vector(float *mat,float *vr,int N){
    float *xr;
    float *xi;
    int a,b;
    xr = (float *)calloc(N * N, sizeof(float));
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
          xr[a]+=mat[a+b*N]*vr[b];
	}
    }
    // Copy back
    copyvec(xr,vr,N);
    free(xr);
}


/* Multiply with double exciton dipole mu_ef on single states */
void dipole_double_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    int seg_num;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;

    for (i = 0; i < non->singles; i++) {
       seg_num=non->psites[i];
        for (j = i + 1; j < non->singles; j++) {
            if (seg_num==non->psites[j]){
              index = Sindex(i, j, non->singles);
              fr[index] += dipole[i] * cr[j];
              fi[index] += dipole[i] * ci[j];
              fr[index] += dipole[j] * cr[i];
              fi[index] += dipole[j] * ci[i];
          }
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on double states */
void dipole_double_inverse_CG2DES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    int seg_num;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < non->singles; i++) fr[i] = 0, fi[i] = 0;

    for (i = 0; i < non->singles; i++) {
        seg_num=non->psites[i];
        for (j = i + 1; j < non->singles; j++) {
            if (seg_num==non->psites[j]){
               index = Sindex(i, j, non->singles);
               fr[j] += dipole[i] * cr[index];
               fi[j] += dipole[i] * ci[index];
               fr[i] += dipole[j] * cr[index];
               fi[i] += dipole[j] * ci[index];
            } 
        }
    }
    return;
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

// calculate the inverse real matrix for the left and right side
void inversie_real_matrix(float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N) {
    int INFO, lwork;
    float *work;
    int i, j;
    int *pivot;
    int M;

    /* Inverse right eigenvectors*/
    pivot = (int *)calloc(N,sizeof(int));
    sgetrf_(&N, &N, ivecR, &N, pivot, &INFO); //LU factorization
    if (INFO != 0) {
        printf("Something went wrong trying to factorize right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; 
    work = (float *)calloc(1, sizeof(float));
    sgetri_(&N, ivecR, &N, pivot, work, &lwork, &INFO); //Find lwork for diagonalization
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    sgetri_(&N, ivecR, &N, pivot, work, &lwork, &INFO); //Do inversion
    if (INFO != 0) {
        printf("Something went wrong trying to inverse right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    free(work), free(pivot);

    /* Inverse left eigenvectors*/
    pivot = (int *)calloc(N,sizeof(int));
    sgetrf_(&N, &N, ivecL, &N, pivot, &INFO); //LU factorization
    if (INFO != 0) {
        printf("Something went wrong trying to factorize left eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; 
    work = (float *)calloc(1, sizeof(float));
    sgetri_(&N, ivecL, &N, pivot, work, &lwork, &INFO); // Find lwork for diagonalization
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    sgetri_(&N, ivecL, &N, pivot, work, &lwork, &INFO); //Do inversion
    if (INFO != 0) {
        printf("Something went wrong trying to inverse left eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    /* Free space */
     free(work), free(pivot);    
}


// calculate the inverse real matrix for the left and right side
void inversie_complex_matrix(float* eig_re, float* eig_im, float* evecL, float* evecR,  float _Complex* ivecL_com,  float _Complex* ivecR_com, int N) {
    int INFO, lwork;
    float *work;
    int i, j;
    int *pivot;
    int M;
    float _Complex *work_com;

    /* Inverse right eigenvectors*/
    for (int i = 0; i < N; i++)
        {// If the current and next eigenvalues form a complex conjugate pair
            if (i < N - 1 && eig_im[i] != 0.0 && eig_im[i + 1] == -eig_im[i]) {
                for (int j = 0; j < N; j++) {
                    // u(j) = VL(:,j) + i*VL(:,j+1)
                    ivecR_com[i * N + j] = evecR[i * N + j]-evecR[(1 + i) * N + j]*_Complex_I;
                    ivecR_com[(1+i) * N + j] = evecR[i * N + j]+evecR[(1 + i) * N + j]*_Complex_I;         
                }
                // Skip the next eigenvector (since it's part of the complex conjugate pair)
                i++;
            } else {
                for (int j = 0; j < N; j++) {
                    // u(j) = VL(:,j)
                    ivecR_com[i * N + j] = evecR[i * N + j]+0*_Complex_I;
                }
            }
        }

    pivot = (int *)calloc(N,sizeof(int));
    cgetrf_(&N, &N, ivecR_com, &N, pivot, &INFO); //LU factorization
    if (INFO != 0) {
        printf("Something went wrong trying to factorize right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; 
    work_com = ( float _Complex *)calloc(1, sizeof( float _Complex));
    cgetri_(&N, ivecR_com, &N, pivot, work_com, &lwork, &INFO); //Find lwork for diagonalization
    lwork = work_com[0];
    free(work_com);
    work_com = ( float _Complex *)calloc(lwork, sizeof( float _Complex));
    cgetri_(&N, ivecR_com, &N, pivot, work_com, &lwork, &INFO); //Do inversion
    if (INFO != 0) {
        printf("Something went wrong trying to inverse right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    free(work_com), free(pivot);

    /* Inverse left eigenvectors*/
    for (int i = 0; i < N; i++)
        {// If the current and next eigenvalues form a complex conjugate pair
            if (i < N - 1 && eig_im[i] != 0.0 && eig_im[i + 1] == -eig_im[i]) {
                for (int j = 0; j < N; j++) {
                    // u(j) = VL(:,j) + i*VL(:,j+1)
                    ivecL_com[i * N + j] = evecL[i * N + j]-evecL[(1 + i) * N + j]*I;
                    ivecL_com[(1+i) * N + j] = evecL[i * N + j]+evecL[(1 + i) * N + j]*I;         
                }
                // Skip the next eigenvector (since it's part of the complex conjugate pair)
                i++;
            } else {
                for (int j = 0; j < N; j++) {
                    // u(j) = VL(:,j)
                    ivecL_com[i * N + j] = evecL[i * N + j]+0*I;
                }
            }
        }

    pivot = (int *)calloc(N,sizeof(int));
    cgetrf_(&N, &N, ivecL_com, &N, pivot, &INFO); //LU factorization
    if (INFO != 0) {
        printf("Something went wrong trying to factorize right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; 
    work_com = ( float _Complex *)calloc(1, sizeof( float _Complex));
    cgetri_(&N, ivecL_com, &N, pivot, work_com, &lwork, &INFO); //Find lwork for diagonalization
    lwork = work_com[0];
    free(work_com);
    work_com = ( float _Complex *)calloc(lwork, sizeof( float _Complex));
    cgetri_(&N, ivecL_com, &N, pivot, work_com, &lwork, &INFO); //Do inversion
    if (INFO != 0) {
        printf("Something went wrong trying to inverse right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }   
        /* Free space */
    free(work_com), free(pivot);
    
    return;
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

/* Write a Hamiltonian to a text file */
void write_ham_to_file(char fname[],float *Hamiltonian_i,int N){
  FILE *file_handle;
  int i,j;
  int index;
  // Build Hamiltonian, convert the triangular Hamiltonian to the  square Hamiltonian
  file_handle=fopen(fname,"w");
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      fprintf(file_handle,"%10.14e ",Hamiltonian_i[Sindex(i,j,N)]);
    }
    fprintf(file_handle,"\n");
  }
  fclose(file_handle);
}
