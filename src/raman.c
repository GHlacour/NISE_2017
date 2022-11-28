#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "absorption.h"
#include "1DFFT.h"
#include "project.h"

// This subroutine is for calculating the Linear Raman
void raman(t_non *non){
  // Initialize variables
  float *VV_re_S_1, *VV_im_S_1, *VH_re_S_1, *VH_im_S_1; // first-order res. func. for parallel(VV) and perpendicular(VH)
  float *alpha,*Hamil_i_e; //variables for raman tensor values (alpha) and hamiltonian
  float *mu_p;
  // Aid arrays
  float *vecr,*veci; //vector real and imaginary

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*alpha_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i;//x runs over the raman indices (xx xy xz yy yz zz)
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
  VV_re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  VV_im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  VH_re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  VH_im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));

  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));

  /* Open Trajectory files */
  H_traj=fopen(non->energyFName,"rb");
  if (H_traj==NULL){
    printf("Hamiltonian file not found!\n");
    exit(1);
  }

  alpha_traj=fopen(non->alphaFName,"rb");
  if (alpha_traj==NULL){
    printf("Raman file %s not found!\n",non->alphaFName);
    exit(1);
  }

  /* Open file with cluster information if appicable */
  if (non->cluster!=-1){
    Cfile=fopen("Cluster.bin","rb");
    if (Cfile==NULL){
      printf("Cluster option was activated but no Cluster.bin file provided.\n");
      printf("Please, provide cluster file or remove Cluster keyword from\n");
      printf("input file.\n");
      exit(0);
    }
    Ncl=0; // Counter for snapshots calculated
  }

  // Here we want to call the routine for checking the trajectory files
  control(non);

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
  mu_p=(float *)calloc(non->singles,sizeof(float));
  alpha=(float *)calloc(non->singles,sizeof(float));

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

      // Configuration belong to cluster
      if (non->cluster==cl){
	      Ncl++;
      }
    }

    // Include snapshot if it is in the cluster or if no clusters are defined
    if (non->cluster==-1 || non->cluster==cl){

      for (x=0;x<6;x++){ //loop over polarizations
        /* Read alpha(ti) */

        if (read_alpha(non,vecr,alpha_traj,ti,x)!=1){
	        printf("Raman trajectory file to short, could not fill buffer!!!\n");
	        printf("ITIME %d %d\n",ti,x);
	        exit(1);
        }

        clearvec(veci,non->singles);
        copyvec(vecr,alpha,non->singles);

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

	        /* Read alpha(tj) */
	        if (read_alpha(non,alpha,alpha_traj,tj,x)!=1){
	          printf("Raman trajectory file to short, could not fill buffer!!!\n");
	         printf("JTIME %d %d\n",tj,x);
	          exit(1);
	        }

	  // Do projection on selected sites if asked
	  if (non->Npsites==0){
             // Find response without projection
            if(x==0 || x==3 || x==5) calc_Ram_VV(VV_re_S_1,VV_im_S_1,t1,non,vecr,veci,alpha);
            if(x==1 || x==2 || x==4) calc_Ram_VH(VH_re_S_1,VH_im_S_1,t1,non,vecr,veci,alpha);
          } else if (non->Npsites<non->singles){
             projection(alpha,non);
             // Find response with projection on single segment
            if(x==0 || x==3 || x==5) calc_Ram_VV(VV_re_S_1,VV_im_S_1,t1,non,vecr,veci,alpha);
            if(x==1 || x==2 || x==4) calc_Ram_VH(VH_re_S_1,VH_im_S_1,t1,non,vecr,veci,alpha);
          } else {
             // Find response with projection on multiple segments
             for (ip=0;ip<pro_dim;ip++){
                multi_projection(alpha,mu_p,non,ip);
                if(x==0 || x==3 || x==5) calc_Ram_VV(VV_re_S_1+non->tmax*ip,VV_im_S_1+non->tmax*ip,t1,non,vecr,veci,mu_p);
                if(x==1 || x==2 || x==4) calc_Ram_VH(VH_re_S_1+non->tmax*ip,VH_im_S_1+non->tmax*ip,t1,non,vecr,veci,mu_p);
             }
	  }

	        // Find response of Raman
          //if(x==0 || x==3 || x==5) calc_Ram(VV_re_S_1,VV_im_S_1,t1,non,vecr,veci,alpha);
          //if(x==1 || x==2 || x==4) calc_Ram_VH(VH_re_S_1,VH_im_S_1,t1,non,vecr,veci,alpha);

	        // Probagate vector
	        if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr,veci,non->ts,1);
	        if (non->propagation==0){
	          if (non->thres==0 || non->thres>1){
	            propagate_vec_DIA(non,Hamil_i_e,vecr,veci,1);
	          } else {
	            elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr,veci,1);
	            if (samples==non->begin){
	              if (t1==0){
		              if (x==0){
		                printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
		                printf("Pressent tuncation %f.\n",non->thres/((non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts)));
		                printf("Suggested truncation %f.\n",0.001);
		              }
	              }
	            }
	          }
          }
        }
      }
    } // Cluster loop

    // Log time
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);

    time_now=log_time(time_now,log);
    fclose(log);
  }

  free(vecr);
  free(veci);
  free(alpha);
  free(mu_p);
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
  outone=fopen("TD_Raman_VV.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    //fprintf(outone,"%f %e %e\n",t1*non->deltat,VV_re_S_1[t1]/samples,VV_im_S_1[t1]/samples);
     fprintf(outone,"%f ",t1*non->deltat);
     for (ip=0;ip<pro_dim;ip++){
        fprintf(outone,"%e %e ",VV_re_S_1[t1+ip*non->tmax]/samples,VV_im_S_1[t1+ip*non->tmax]/samples);
     }
     fprintf(outone,"\n");
  }
  fclose(outone);

    outone=fopen("TD_Raman_VH.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    //fprintf(outone,"%f %e %e\n",t1*non->deltat,2*VH_re_S_1[t1]/samples,2*VH_im_S_1[t1]/samples);
     fprintf(outone,"%f ",t1*non->deltat);
     for (ip=0;ip<pro_dim;ip++){
        fprintf(outone,"%e %e ",VH_re_S_1[t1+ip*non->tmax]/samples,VH_im_S_1[t1+ip*non->tmax]/samples);
     }
     fprintf(outone,"\n");
  }
  fclose(outone);

  /* Do Forier transform and save */
  do_1DFFT(non,"Raman_VV.dat",VV_re_S_1,VV_im_S_1,samples);
  do_1DFFT(non,"Raman_VH.dat",VH_re_S_1,VH_im_S_1,samples);

  free(VV_re_S_1),free(VV_im_S_1),free(VH_re_S_1),free(VH_im_S_1);

  printf("----------------------------------------------\n");
  printf(" Raman calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}

// This subroutine updates the response function for VV and VH
void calc_Ram_VV(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=alpha[i]*cr[i];
    im_S_1[t1]+=alpha[i]*ci[i];
  }
  return;
}

void calc_Ram_VH(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=2*alpha[i]*cr[i];
    im_S_1[t1]+=2*alpha[i]*ci[i];
  }
  return;
}