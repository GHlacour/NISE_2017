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
#include "luminescence.h"
#include "1DFFT.h"
#include "project.h"

void luminescence(t_non *non){
  // Initialize variables
  float *re_S_1,*im_S_1; // The first-order response function
  float *re_S_1x,*im_S_1x;
  float *re_S_1y,*im_S_1y;
  float *re_S_1z,*im_S_1z;
  float *mu_eg,*Hamil_i_e;
  float *mu_xyz;
  // Aid arrays
  float *vecr,*veci;
  //,*vecr_old,*veci_old;

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *outone,*log;
  FILE *C_traj;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i;
  int t1,fft;
  int elements;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  fft=0; 
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;
  printf("Temperature %f.\n",non->temperature);

  // Allocate memory
  re_S_1=(float *)calloc(non->tmax,sizeof(float));
  im_S_1=(float *)calloc(non->tmax,sizeof(float));
  re_S_1x=(float *)calloc(non->tmax,sizeof(float));
  im_S_1x=(float *)calloc(non->tmax,sizeof(float));
  re_S_1y=(float *)calloc(non->tmax,sizeof(float));
  im_S_1y=(float *)calloc(non->tmax,sizeof(float));
  re_S_1z=(float *)calloc(non->tmax,sizeof(float));
  im_S_1z=(float *)calloc(non->tmax,sizeof(float));
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));

  /* Open Trajectory files */
  H_traj=fopen(non->energyFName,"rb");
  if (H_traj==NULL){
    printf("Hamiltonian file not found!\n");
    exit(1);
  }

  mu_traj=fopen(non->dipoleFName,"rb");
  if (mu_traj==NULL){
    printf("Dipole file %s not found!\n",non->dipoleFName);
    exit(1);
  }

  // Here we want to call the routine for checking the trajectory files

  itime=0;
  // Do calculation
  N_samples=(non->length-non->tmax1-1)/non->sample+1;
  if (N_samples>0) {
    printf("Making %d samples!\n",N_samples);
  } else {
    printf("Insufficient data to calculate spectrum.\n");
    printf("Please, lower max times or provide longer\n");
    printf("trajectory.\n");
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

  vecr=(float *)calloc(non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles,sizeof(float));
//  vecr_old=(float *)calloc(non->singles,sizeof(float));
//  veci_old=(float *)calloc(non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_xyz=(float *)calloc(non->singles*3,sizeof(float));

/* Read coupling, this is done if the coupling and transition-dipoles are *
 *    * time-independent and only one snapshot is stored */
  if (!strcmp(non->hamiltonian,"Coupling")){
    C_traj=fopen(non->couplingFName,"rb");
    if (C_traj==NULL){
      printf("Coupling file not found!\n");
      exit(1);
    }
    if (read_He(non,Hamil_i_e,C_traj,-1)!=1){
      printf("Coupling trajectory file to short, could not fill buffer!!!\n");
      exit(1);
    }
    fclose(C_traj);
    /* Reading in single fixed transition dipole vector matrix */
    for (x=0;x<3;x++){
      if (read_mue(non,mu_xyz+non->singles*x,mu_traj,0,x)!=1){
         printf("Dipole trajectory file to short, could not fill buffer!!!\n");
         printf("ITIME %d %d\n",0,x);
         exit(1);
      }
    }
  }


  // Loop over samples
  for (samples=non->begin;samples<non->end;samples++){

    // Calculate linear response    
    ti=samples*non->sample;
    for (x=0;x<3;x++){
      /* Read mu(ti) */
      if (!strcmp(non->hamiltonian,"Coupling")){
          copyvec(mu_xyz+non->singles*x,vecr,non->singles);
      } else {
          if (read_mue(non,vecr,mu_traj,ti,x)!=1){
	      printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	      printf("ITIME %d %d\n",ti,x);
	      exit(1);
          }
      }
      clearvec(veci,non->singles);
//      copyvec(vecr,vecr_old,non->singles);
      copyvec(vecr,mu_eg,non->singles);

      // Add Boltzman weight
      bltz_weight(mu_eg,Hamil_i_e,non);

      // Loop over delay
      for (t1=0;t1<non->tmax;t1++){
	tj=ti+t1;
	/* Read Hamiltonian */
        if (!strcmp(non->hamiltonian,"Coupling")){
            if (read_Dia(non,Hamil_i_e,H_traj,tj)!=1){
              printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
              exit(1);
            }
        } else {
	    if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	        exit(1);
           }
	}
	
        /* Read mu(tj) */
        if (!strcmp(non->hamiltonian,"Coupling")){
          copyvec(mu_xyz+non->singles*x,mu_eg,non->singles);
        } else {
	    if (read_mue(non,mu_eg,mu_traj,tj,x)!=1){
	      printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	      printf("JTIME %d %d\n",tj,x);
	      exit(1);
	    }
        }

	// Do projection on selected sites if asked
	if (non->Npsites>0){
	  projection(mu_eg,non);
	}

	/* Add Boltzman weight */
	if (non->is==0){
	    bltz_weight(mu_eg,Hamil_i_e,non);
	} else {
            bltz_weight_itime(mu_eg,Hamil_i_e,non);
        }	    
	// Find response
	calc_S1(re_S_1,im_S_1,t1,non,vecr,veci,mu_eg);
        if (x==0) calc_S1(re_S_1x,im_S_1x,t1,non,vecr,veci,mu_eg);
        if (x==1) calc_S1(re_S_1y,im_S_1y,t1,non,vecr,veci,mu_eg);
        if (x==2) calc_S1(re_S_1z,im_S_1z,t1,non,vecr,veci,mu_eg);

	/* Probagate vector */
        propagate_vector(non,Hamil_i_e,vecr,veci,1,samples,t1*x);
	
      }
    }
  
    // Log time
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);
    time_now=log_time(time_now,log);
    fclose(log);
  }

  free(vecr);
  free(veci);
//  free(vecr_old);
//  free(veci_old);
  free(mu_eg);
  free(mu_xyz);
  free(Hamil_i_e);

  // The calculation is finished, lets write output
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  fclose(mu_traj),fclose(H_traj);

  outone=fopen("TD_Lum.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f %e %e\n",t1*non->deltat,re_S_1[t1]/samples,im_S_1[t1]/samples);
  }
  fclose(outone);

  /* Do Forier transform and save */
  do_1DFFT(non,"Luminescence.dat",re_S_1,im_S_1,samples);
  do_1DFFT(non,"Luminescence_x.dat",re_S_1x,im_S_1x,samples);
  do_1DFFT(non,"Luminescence_y.dat",re_S_1y,im_S_1y,samples);
  do_1DFFT(non,"Luminescence_z.dat",re_S_1z,im_S_1z,samples);

  free(re_S_1),free(im_S_1);
  free(re_S_1x),free(im_S_1x);
  free(re_S_1y),free(im_S_1y);
  free(re_S_1z),free(im_S_1z);

  printf("----------------------------------------------\n");
  printf(" Luminescence calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

void calc_LUM(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=mu[i]*cr[i];
    im_S_1[t1]+=mu[i]*ci[i];
  }
  return;
}

void bltz_weight(float *mu_eg,float *Hamiltonian_i,t_non *non){
  int index,N;
  float *H,*e,*c2,*c1;
  float *cnr;
  float *crr;
  N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  c2=(float *)calloc(N,sizeof(float));
  c1=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N*N,sizeof(float));
  crr=(float *)calloc(N*N,sizeof(float));
  int a,b,c;
  float kBT=non->temperature*0.695; // Kelvin to cm-1
  float Q,iQ;

  // Build Hamiltonian
  for (a=0;a<N;a++){
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; // Diagonal
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }
  diagonalizeLPD(H,e,N);
 
  // Exponentiate [U=exp(-H/kBT)]
  for (a=0;a<N;a++){
    c2[a]=exp(-e[a]/kBT);
    Q=Q+c2[a];
  }
  iQ=1.0/Q;

  // Transform to site basis
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      cnr[b+a*N]+=H[b+a*N]*c2[b]*iQ;
    }
  }  
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
        crr[a+c*N]+=H[b+a*N]*cnr[b+c*N];
      }
    }
  }
  // The weights in the site basis were calculated

  // Multiply weights to dipole operator
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      c1[a]+=crr[a+b*N]*mu_eg[b];
    }
  }

  // Update dipole operator
  for (a=0;a<N;a++){
    mu_eg[a]=c1[a];
  }
  free(H);
  free(c2);
  free(c1);
  free(e);
  free(crr);
  free(cnr);
  return;
}

/* This routine multiplies with the equilibrium density matrix through */
/* imaginary time propagation */
void bltz_weight_itime(float *cr,float *Hamiltonian_i,t_non *non){
  float f;
  int index, N;
  float *H1, *H0, *re_U;
  int *col, *row;
  float *ocr;
  float J;
  float eJ,emJ;
  float cr1, cr2, ci1, ci2;
  float coh, sih;
  int i, k, kmax;
  int a,b,c;
  float kBT=non->temperature*0.695; // Kelvin to cm-1
  float Q,iQ,norm;
  int m;
  m=non->is;

  N = non->singles;
  f = 1.0 /(m*kBT);
  H0 = (float *)calloc(N, sizeof(float));
  H1 = (float *)calloc(N * N, sizeof(float));
  col = (int *)calloc(N * N / 2, sizeof(int));
  row = (int *)calloc(N * N / 2, sizeof(int));
  re_U = (float *)calloc(N, sizeof(float));
  ocr = (float *)calloc(N, sizeof(float));
  Q=0;
  norm=0;

  /* Find initial norm */
  for (a=0;a<N;a++){
    norm+=cr[a]*cr[a];
  } 

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

    /* Exponentiate diagonal [U=exp(-H0/2kBT )] */
    for (a = 0; a < N; a++) {
        re_U[a] = exp(-0.5 * H0[a] * f);
    }

    for (i = 0; i < m; i++) {
      /* Multiply on vector first time */
      for (a = 0; a < N; a++) {
         ocr[a] = cr[a] * re_U[a];
      }

      /* Account for couplings */
      for (k = 0; k < kmax; k++) {
            a = col[k];
            b = row[k];
            J = H1[k];
            J = -J * f;
	    eJ=exp(J)*0.5;
	    emJ=exp(-J)*0.5;
            sih = eJ-emJ;
            coh = eJ+emJ;
            cr1 = coh * ocr[a]+sih * ocr[b];
            cr2 = coh * ocr[b]+sih * ocr[a];
            ocr[a] = cr1, ocr[b] = cr2;
        }
              /* Multiply on vector second time */
        for (a = 0; a < N; a++) {
            cr[a] = ocr[a] * re_U[a] ;
        }
    }

  /* Find final norm */
  for (a=0;a<N;a++){
    Q+=cr[a]*cr[a];
  }
  iQ=sqrt(norm/Q);
  /* Renormalize (divide by partition function) */
  for (a=0;a<N;a++){
    cr[a]=cr[a]*iQ;
  }

    free(ocr), free(re_U), free(H1), free(H0);
    free(col), free(row);
         
}
