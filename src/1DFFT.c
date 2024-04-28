#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "types.h"
#include "nrutil.h"
#include "project.h"
#include "1DFFT.h"
#include "NISE_subs.h"
#include <stdio.h>
#include <stdlib.h>
#include "propagate.h"
#include "read_trajectory.h"
#include "omp.h"


void Lineshape_FFT(t_non *non){
  FILE *outone,*outone1;
  int samples,i;
  int pro_dim,a,ip;
  int loop_num = 0;
  /* Floats */
  float shift1;
  float *re_S_1,*im_S_1; /* The first-order response function */
  float *re_S_2,*im_S_2; /* The first-order response function */
  float *re_S_fin,*im_S_fin; /* The first-order response function */
  float dummyf;
    /* Check if projection keyword is defined and select number of subsets in */
  /* projection */
  samples =1;
  /* Allocate memory */
  pro_dim=project_dim(non);
  re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  re_S_2=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_2=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  re_S_fin=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_fin=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;
  /* Open Trajectory files */
  outone=fopen("TD_Absorption.dat","r");
  if (outone==NULL){
    printf("TD_Absorption.dat file not found!\n");
    exit(1);
  }
  for (i=0;i<non->tmax;i++){
    fscanf(outone,"%f",&dummyf);
      ip=0;
      fscanf(outone,"%f %f",&re_S_1[i+non->tmax*ip],&im_S_1[i+non->tmax*ip]);
  }    
  printf("\n\nCompleted reading the TD_Absorption.dat file.\n");
  
  outone1=fopen("lineshape.dat","r");
  if (outone1==NULL){
    printf("lineshape.dat file not found!\n");
    exit(1);
  }
      // Read data from the file
  for (i=0;i<non->tmax;i++){
    fscanf(outone1,"%f",&dummyf);
    fscanf(outone1,"%f %f",&re_S_2[i],&im_S_2[i]);
  }  
  printf("\n\nCompleted reading the lineshape.dat file.\n");

  for (a = 0; a < non->tmax; a++) {
    re_S_fin[a]=re_S_1[a]*re_S_2[a]-im_S_1[a]*im_S_2[a];
    im_S_fin[a]=re_S_1[a]*im_S_2[a]+im_S_1[a]*re_S_2[a];
  }
  printf("Performing the 1 D Forier transform and save.\n"); 

  do_1DFFT(non,"Absorption_vibronic.dat",re_S_fin,im_S_fin,samples);
    // Close the file
    fclose(outone);
    fclose(outone1);

  /* Free response function arrary */
  free(re_S_1),free(im_S_1);
  free(re_S_2),free(im_S_2);
  free(re_S_fin),free(im_S_fin); 

  printf("----------------------------------------------\n");
  printf(" Absorption calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");
  return;
}

void ONE_DFFT(t_non *non){
  FILE *outone;
  /* Floats */
  float shift1;
  int samples;
  int pro_dim,a,i,ip;
  int loop_num = 0;
  float *re_S_1,*im_S_1; /* The first-order response function */
  float dummyf;
    /* Check if projection keyword is defined and select number of subsets in */
  /* projection */
  samples=1;
  pro_dim=project_dim(non);
  /* Find average frequency for shifting the frequencies and remove */ 
  /* high-frequency oscillations in the response function */
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  /* Allocate memory */
  re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));

  /* Open Trajectory files */
  outone=fopen("TD_Absorption.dat","r");
  if (outone==NULL){
    printf("TD_Absorption.dat file not found!\n");
    exit(1);
  }
  for (i=0;i<non->tmax;i++){
    fscanf(outone,"%f",&dummyf);
      ip=0;
      fscanf(outone,"%f %f",&re_S_1[i+non->tmax*ip],&im_S_1[i+non->tmax*ip]);
  }
  printf("\n\nCompleted reading the TD_Absorption.dat file.\n");
  printf("Performing the 1 D Forier transform and save.\n"); 
  do_1DFFT(non,"Absorption.dat",re_S_1,im_S_1,samples);

  // Close the file
  fclose(outone);
  /* Free response function arrary */
  free(re_S_1),free(im_S_1);

  printf("----------------------------------------------\n");
  printf(" Absorption calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");
  return;
}


/* Do 1D Fourier transform */
void do_1DFFT(t_non *non,char fname[],float *re_S_1,float *im_S_1,int samples){

  /* Floats */
  float shift1;
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;
  float *spec_r,*spec_i;
  /* Integers */
  int i,fft;
  int pro_dim,ip;
  /* Files */
  FILE *outone;

  shift1=non->shifte;
  pro_dim=project_dim(non);
  fft=non->fft;
  if (fft<non->tmax1) fft=non->tmax1;

    /* Fourier transform 1D spectrum */
  fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*2));
  fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*2));
  fftPlan = fftw_plan_dft_1d(fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);
  
  spec_r=(float *)calloc(fft*2*pro_dim,sizeof(float));  
  spec_i=(float *)calloc(fft*2*pro_dim,sizeof(float));
  for (ip=0;ip<pro_dim;ip++){  
     for (i=0;i<=fft;i++){
        fftIn[i][0]=0;
        fftIn[i][1]=0;
     }
     for (i=0;i<non->tmax1;i++){
        fftIn[i][0]=im_S_1[i+non->tmax*ip]/samples;
        fftIn[i][1]=re_S_1[i+non->tmax*ip]/samples;
        if (non->lifetime > 0.0){
           fftIn[i][0]*=exp(-i*non->deltat/(2*non->lifetime));
           fftIn[i][1]*=exp(-i*non->deltat/(2*non->lifetime));
        }
        if (non->homogen > 0.0){
           fftIn[i][0]*=exp(-i*non->deltat/(2*non->homogen));
           fftIn[i][1]*=exp(-i*non->deltat/(2*non->homogen));
        }
        if (non->inhomogen > 0.0){
           fftIn[i][0]*=exp(-i*non->deltat*i*non->deltat/(2*non->inhomogen*non->inhomogen));
           fftIn[i][1]*=exp(-i*non->deltat*i*non->deltat/(2*non->inhomogen*non->inhomogen));
        }
/*    fftIn[fft-i][0]=-im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    fftIn[fft-i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));*/
     }
  /* Scale first point */
     fftIn[0][0]=fftIn[0][0]*0.5;
     fftIn[0][1]=fftIn[0][1]*0.5;

     fftw_execute(fftPlan);
     for (i=0;i<2*fft;i++){
          spec_r[i+fft*2*ip]=fftOut[i][1];
          spec_i[i+fft*2*ip]=fftOut[i][0];
     }
  }
  outone=fopen(fname,"w");
  for (i=fft/2;i<=fft-1;i++){
    if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){
//      fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
      fprintf(outone,"%f ",-((fft-i)/non->deltat/c_v/fft-shift1));
      for (ip=0;ip<pro_dim;ip++){
         fprintf(outone,"%e %e ",spec_r[i+fft*2*ip],spec_i[i+fft*2*ip]);
      }
      fprintf(outone,"\n");
    }
  }
  for (i=0;i<=fft/2-1;i++){
    if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){
//      fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
      fprintf(outone,"%f ",-((-i)/non->deltat/c_v/fft-shift1));
      for (ip=0;ip<pro_dim;ip++){
          fprintf(outone,"%e %e ",spec_r[i+fft*2*ip],spec_i[i+fft*2*ip]);
      }
      fprintf(outone,"\n");
    }
  }
    
  fclose(outone);
}

/* Do 1D Fourier transform */
void do_1DFFTold(t_non *non,char fname[256],float *re_S_1,float *im_S_1,int samples){

  /* Floats */
  float shift1;
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;
  /* Integers */
  int i,fft;
  /* Files */
  FILE *outone;

  shift1=non->shifte;

  fft=0;
  if (fft<non->tmax1*2) fft=2*non->tmax1;
 
  /* Fourier transform 1D spectrum */
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
  outone=fopen(fname,"w");
  for (i=fft/2;i<=fft-1;i++){
    if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
      fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
/*        fprintf(outone,"%f ",-((fft-i)/non->deltat/c_v/fft-shift1));
        for (ip=0;ip<pro_dim;ip++){
	  fprintf(outone,"%e %e ",spec_r[i+fft*2*ip],spec_i[i+fft*2*ip]);
        }
        fprintf(outone,"\n");*/
    }
  }
  for (i=0;i<=fft/2-1;i++){
    if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
      fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
/*      fprintf(outone,"%f ",-((-i)/non->deltat/c_v/fft-shift1));
      for (ip=0;ip<pro_dim;ip++){
          fprintf(outone,"%e %e ",spec_r[i+fft*2*ip],spec_i[i+fft*2*ip]);
      }
      fprintf(outone,"\n");*/
    }
  }
    
  fclose(outone);
}
