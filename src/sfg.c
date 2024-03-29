#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "read_trajectory.h"
#include "propagate.h"
#include "sfg.h"
#include "1DFFT.h"
#include "project.h"

/* This subroutine is for calculating the Linear SFG 
 * We follow the theory of for example Roy, Groenbaum, and Skinner
 * J. Chem. Phys. 141 18C502 (2024) */

/* The surface normal is assumed to be along the z-axis
 * Then the SSP signal corresponds to 1/2 (\chi_{xxz)}+\chi_{yyz})
 * and the PPP signal is chi_{zzz} */

void sfg(t_non *non){
  /* Initialize variables */
  float *SSP_re_S_1, *SSP_im_S_1, *PPP_re_S_1, *PPP_im_S_1; /* Response functions */
  float *alpha,*Hamil_i_e; /* Variables for Raman tensor values (alpha) and Hamiltonian */
  float *mu_eg,*mu_xyz; 
  float *mu_p;

  /* Auxillary arrays for propagation */
  float *vecr,*veci; /* vector real and imaginary */

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*alpha_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i; /* x runs over the Raman indices (xx xy xz yy yz zz) */
  int y; /* y run over dipole indices (x y z) */
  int t1,fft;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  // Check for projection
  pro_dim=project_dim(non);

  // Allocate memory
  SSP_re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  SSP_im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  PPP_re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  PPP_im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));

  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_p=(float *)calloc(non->singles,sizeof(float)); 
  mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
  vecr=(float *)calloc(non->singles,sizeof(float));
  veci=(float *)calloc(non->singles,sizeof(float));
  alpha=(float *)calloc(non->singles,sizeof(float));

  /* Open Trajectory files */
  open_files(non,&H_traj,&mu_traj,&Cfile);

  alpha_traj=fopen(non->alphaFName,"rb");
  if (alpha_traj==NULL){
    printf("Polarizability file %s not found!\n",non->alphaFName);
    exit(1);
  }

  /* Initialize sample numbers */
  N_samples=determine_samples(non);
  Ncl=0;

  // Here we want to call the routine for checking the trajectory files
  control(non);


  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

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

  /* Read coupling, this is done if the coupling and transition-dipoles are *
   * time-independent and only one snapshot is stored */
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
  }

  /* Loop over samples */
  for (samples=non->begin;samples<non->end;samples++){

    /* Calculate linear response */
    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
	      printf("Cluster trajectory file to short, could not fill buffer!!!\n");
	      printf("ITIME %d\n",ti);
	      exit(1);
      }

      /* Configuration belong to cluster */
      if (non->cluster==cl){
	      Ncl++;
      }
    }

    /* Include snapshot if it is in the cluster or if no clusters are defined */
    if (non->cluster==-1 || non->cluster==cl){

      for (y=0;y<3;y++){ /* Loop over polarizations */
        /* Read mu(ti) into auxillary array vecr and clear the imaginary part */
        read_dipole(non,mu_traj,vecr,mu_xyz,y,ti);
        clearvec(veci,non->singles);

       /* Loop over coherence time */
       for (t1=0;t1<non->tmax;t1++){
	   tj=ti+t1;
	   /* Read Hamiltonian */
	   read_Hamiltonian(non,Hamil_i_e,H_traj,tj);

	   for (x=0;x<6;x++){
	       /* Read alpha(tj) */
	       if (read_alpha(non,alpha,alpha_traj,tj,x)!=1){
	           printf("Polarizability trajectory file to short, could not fill buffer!!!\n");
	           printf("JTIME %d %d\n",tj,x);
	           exit(1);
	       }

	       /* Do projection on selected sites if asked */
	       if (non->Npsites==0){
                   /* Find response without projection */
                   if(x==5 && y==2) calc_SFG_PPP(PPP_re_S_1,PPP_im_S_1,t1,non,vecr,veci,alpha);
                   if((x==0 || x==3) && y==2) calc_SFG_SSP(SSP_re_S_1,SSP_im_S_1,t1,non,vecr,veci,alpha);
               } else if (non->Npsites<non->singles){
                   projection(alpha,non);
                   /* Find response with projection on single segment */
                   if(x==5 && y==2) calc_SFG_PPP(PPP_re_S_1,PPP_im_S_1,t1,non,vecr,veci,alpha);
                   if((x==0 || x==3) && y==2) calc_SFG_SSP(SSP_re_S_1,SSP_im_S_1,t1,non,vecr,veci,alpha);
               } else {
                   /* Find response with projection on multiple segments */
                   for (ip=0;ip<pro_dim;ip++){
                       multi_projection(alpha,mu_p,non,ip);
		       if(x==5 && y==2) calc_SFG_PPP(PPP_re_S_1+non->tmax*ip,PPP_im_S_1+non->tmax*ip,t1,non,vecr,veci,alpha);
                       if((x==0 || x==3) && y==2) calc_SFG_SSP(SSP_re_S_1+non->tmax*ip,SSP_im_S_1+non->tmax*ip,t1,non,vecr,veci,alpha);
                   }
	       }
	    }
	    /* Propagate vector */
	    propagate_vector(non,Hamil_i_e,vecr,veci,1,samples,t1*y);
         }
       }
    } /* Cluster loop */

    /* Log time */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);

    time_now=log_time(time_now,log);
    fclose(log);
  }

  free(vecr);
  free(veci);
  free(alpha);
  free(mu_p);
  free(mu_xyz);
  free(mu_eg);
  free(Hamil_i_e);

  // The calculation is finished, lets write output
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");
  fclose(log);

  if (non->cluster!=-1){
    printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
    //    printf("Changing normalization to number of samples used.");
    //    samples=Ncl;
    if (samples==0){ // Avoid dividing by zero
      samples=1;
    }
  }

  fclose(alpha_traj),fclose(H_traj);
  if (non->cluster!=-1){
    fclose(Cfile);
  }

  /* Save time domain response */
  outone=fopen("TD_SFG_PPP.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    //fprintf(outone,"%f %e %e\n",t1*non->deltat,VV_re_S_1[t1]/samples,VV_im_S_1[t1]/samples);
     fprintf(outone,"%f ",t1*non->deltat);
     for (ip=0;ip<pro_dim;ip++){
        fprintf(outone,"%e %e ",PPP_re_S_1[t1+ip*non->tmax]/samples,PPP_im_S_1[t1+ip*non->tmax]/samples);
     }
     fprintf(outone,"\n");
  }
  fclose(outone);

  outone=fopen("TD_SFG_SSP.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    //fprintf(outone,"%f %e %e\n",t1*non->deltat,2*VH_re_S_1[t1]/samples,2*VH_im_S_1[t1]/samples);
     fprintf(outone,"%f ",t1*non->deltat);
     for (ip=0;ip<pro_dim;ip++){
        fprintf(outone,"%e %e ",SSP_re_S_1[t1+ip*non->tmax]/samples,SSP_im_S_1[t1+ip*non->tmax]/samples);
     }
     fprintf(outone,"\n");
  }
  fclose(outone);

  /* Do Forier transform and save */
  do_1DFFT(non,"SFG_PPP.dat",PPP_re_S_1,PPP_im_S_1,samples);
  do_1DFFT(non,"SFG_SSP.dat",SSP_re_S_1,SSP_im_S_1,samples);

  free(PPP_re_S_1),free(PPP_im_S_1),free(SSP_re_S_1),free(SSP_im_S_1);

  printf("----------------------------------------------\n");
  printf(" SFG calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}

/* This subroutine updates the response function for PPP */
void calc_SFG_PPP(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=alpha[i]*cr[i];
    im_S_1[t1]+=alpha[i]*ci[i];
  }
  return;
}

/* This subroutine updates the response function for SSP */
void calc_SFG_SSP(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=alpha[i]*cr[i]/2;
    im_S_1[t1]+=alpha[i]*ci[i]/2;
  }
  return;
}
