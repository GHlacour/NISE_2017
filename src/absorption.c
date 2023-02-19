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
#include "absorption.h"
#include "1DFFT.h"
#include "project.h"

/* This subroutine is for calculating the Linear Absorption */
/* both for IR and UVvis */
void absorption(t_non *non){
  /* Initialize variables */
  float *re_S_1,*im_S_1; /* The first-order response function */
  float *mu_eg,*Hamil_i_e;
  float *mu_p;
  float *mu_xyz;
  /* Auxillary arrays */
  float *vecr,*veci;

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i;
  int t1,fft;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);

  /* Find average frequency for shifting the frequencies and remove */ 
  /* high-frequency oscillations in the response function */
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  /* Check if projection keyword is defined and select number of subsets in */
  /* projection */
  pro_dim=project_dim(non);

  /* Allocate memory */
  re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  vecr=(float *)calloc(non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_p=(float *)calloc(non->singles,sizeof(float)); 
  mu_xyz=(float *)calloc(non->singles*3,sizeof(float));

  /* Open Trajectory files */
  open_files(non,&H_traj,&mu_traj,&Cfile);
  
  /* Here we want to call the routine for checking the trajectory files */
  /* before we start the calculation */   
  control(non);

  itime=0;
  /* Initialize sample numbers */
  N_samples=determine_samples(non);
  Ncl=0;

  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

  /* Loop over samples */
  for (samples=non->begin;samples<non->end;samples++){

    /* Calculate linear response */   
    ti=samples*non->sample;

    /* If clusters are selectred open cluster information file */
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
      /* Loop over the three polarization directions */  
      for (x=0;x<3;x++){
        /* Read mu(ti) into auxillary array vecr and clear the imaginary part */
        read_dipole(non,mu_traj,vecr,mu_xyz,x,ti);        
        clearvec(veci,non->singles);

        /* Loop over coherence time */
        for (t1=0;t1<non->tmax;t1++){
	        tj=ti+t1;
	        /* Read Hamiltonian */
          read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
	
  	      /* Read mu(tj) */
          read_dipole(non,mu_traj,mu_eg,mu_xyz,x,tj);

	        /* Do projection on selected sites if asked and calculate response function */
	        if (non->Npsites==0){
             /* Find response without projection */
             calc_S1(re_S_1,im_S_1,t1,non,vecr,veci,mu_eg);
          } else if (non->Npsites<non->singles){
            projection(mu_eg,non);
            /* Find response with projection on single segment */
            calc_S1(re_S_1,im_S_1,t1,non,vecr,veci,mu_eg);
          } else {
            /* Find response with projection on multiple segments */
            for (ip=0;ip<pro_dim;ip++){
              multi_projection(mu_eg,mu_p,non,ip);
              calc_S1(re_S_1+non->tmax*ip,im_S_1+non->tmax*ip,t1,non,vecr,veci,mu_p);
            }
	        }

	        /* Probagate vector */
          propagate_vector(non,Hamil_i_e,vecr,veci,1,samples,t1*x);
        }
      }
    } /* Loop over possible Clusters */
  
    /* Update Log file with time and sample numner */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);        
    time_now=log_time(time_now,log);
    fclose(log);
  }

  /* The calculation is finished we can close all auxillary arrays before writing */
  /* output to file. */
  free(vecr);
  free(veci);
  free(mu_eg);
  free(mu_xyz);
  free(Hamil_i_e);

  /* The calculation is finished, lets write output */
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  /* Print information on number of realizations included belonging to the selected */
  /* cluster and close the cluster file. (Only to be done if cluster option is active.) */
  if (non->cluster!=-1){
      printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
      fclose(Cfile);
  }

  /* Close Trajectory Files */
  fclose(mu_traj),fclose(H_traj);

  /* Save time domain response */
  outone=fopen("TD_Absorption.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
     fprintf(outone,"%f ",t1*non->deltat);
     for (ip=0;ip<pro_dim;ip++){
        fprintf(outone,"%e %e ",re_S_1[t1+ip*non->tmax]/samples,im_S_1[t1+ip*non->tmax]/samples);
     }
     fprintf(outone,"\n");
  }
  fclose(outone);

  /* Do Forier transform and save */
  do_1DFFT(non,"Absorption.dat",re_S_1,im_S_1,samples);

  /* Free response function arrary */
  free(re_S_1),free(im_S_1);

  printf("----------------------------------------------\n");
  printf(" Absorption calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

/* This subroutine updates the response function */
void calc_S1(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=mu[i]*cr[i];
    im_S_1[t1]+=mu[i]*ci[i];
  }
  return;
}
