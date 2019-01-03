#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types3.1.h"
#include "NISE3.1subs.h"
#include "absorption.h"

void absorption(t_non *non){
  // Initialize variables
  float *re_S_1,*im_S_1; // The first-order response function
  float *mu_eg,*Hamil_i_e;

  // Aid arrays
  float *vecr,*veci,*vecr_old,*veci_old;

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*mu_traj;
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

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  // Allocate memory
  re_S_1=(float *)calloc(non->tmax,sizeof(float));
  im_S_1=(float *)calloc(non->tmax,sizeof(float));
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

  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);

  vecr=(float *)calloc(non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles,sizeof(float));
  vecr_old=(float *)calloc(non->singles,sizeof(float));
  veci_old=(float *)calloc(non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));

  // Loop over samples
  for (samples=non->begin;samples<non->end;samples++){

    // Calculate linear response    
    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
	printf("Cluster trajectory file to short, could not fill buffer!!!\n");
	printf("ITIME %d\n",ti);
	exit(1);
      }
      //      printf("%d\n",cl);
      // Configuration belong to cluster
      if (non->cluster==cl){
	Ncl++;
      }
    }
    if (non->cluster==-1 || non->cluster==cl){
      
    for (x=0;x<3;x++){
      // Read mu(ti)
      if (read_mue(non,vecr,mu_traj,ti,x)!=1){
	printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	printf("ITIME %d %d\n",ti,x);
	exit(1);
      }
      clearvec(veci,non->singles);
      copyvec(vecr,vecr_old,non->singles);
      copyvec(vecr,mu_eg,non->singles);
      // Loop over delay
      for (t1=0;t1<non->tmax;t1++){
	tj=ti+t1;
	// Read Hamiltonian
	if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	  exit(1);
	}
	
	// Read mu(tj)
	if (read_mue(non,mu_eg,mu_traj,tj,x)!=1){
	  printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	  printf("JTIME %d %d\n",tj,x);
	  exit(1);
	}

	// Do projection on selected sites if asked
	if (non->Npsites>0){
	  projection(mu_eg,non);
	}
	
	// Find response
	calc_S1(re_S_1,im_S_1,t1,non,vecr,veci,mu_eg);

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
  free(vecr_old);
  free(veci_old);
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

  fclose(mu_traj),fclose(H_traj);
  if (non->cluster!=-1){
    fclose(Cfile);
  }

  outone=fopen("TD_Absorption.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f %e %e\n",t1*non->deltat,re_S_1[t1]/samples,im_S_1[t1]/samples);
  }
  fclose(outone);

  fft=0;
  if (fft<non->tmax1*2) fft=2*non->tmax1;
 
  // Fourier transform 1D spectrum
  fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
  fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
  fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);
    
  for (i=0;i<=fft;i++){
    fftIn[i][0]=0;
    fftIn[i][1]=0;
  }
  for (i=0;i<non->tmax1;i++){
    fftIn[i][0]=im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    fftIn[i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    fftIn[fft-i][0]=-im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    fftIn[fft-i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
  }

  fftw_execute(fftPlan);
  outone=fopen("Absorption.dat","w");
  for (i=fft/2;i<=fft-1;i++){
    if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
      fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
    }
  }
  for (i=0;i<=fft/2-1;i++){
    if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
      fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]*0);
    }
  }
    
  fclose(outone);
  free(re_S_1),free(im_S_1);

  printf("----------------------------------------------\n");
  printf(" Absorption calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

void calc_S1(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu){
  int i;
  for (i=0;i<non->singles;i++){
    re_S_1[t1]+=mu[i]*cr[i];
    im_S_1[t1]+=mu[i]*ci[i];
  }
  return;
}
