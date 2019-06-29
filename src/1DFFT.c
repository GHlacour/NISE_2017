#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "types.h"
#include "nrutil.h"

/* Do 1D Fourier transform */
void do_1DFFT(t_non *non,char fname[256],float *re_S_1,float *im_S_1,int samples){

  /* Floats */
  float shift1;
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;
  /* Integers */
  int i,fft;
  /* Files */
  FILE *outone;

  shift1=non->shifte;

  fft=non->fft;
  if (fft<non->tmax1) fft=non->tmax1;

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
/*    fftIn[fft-i][0]=-im_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));
    fftIn[fft-i][1]=re_S_1[i]/samples*exp(-i*non->deltat/(2*non->lifetime));*/
  }
  /* Scale first point */
  fftIn[0][0]=fftIn[0][0]*0.5;
  fftIn[0][1]=fftIn[0][1]*0.5;

  fftw_execute(fftPlan);
  outone=fopen(fname,"w");
  for (i=fft/2;i<=fft-1;i++){
    if (-((fft-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((fft-i)/non->deltat/c_v/fft-shift1)<non->max1){
      fprintf(outone,"%f %e %e\n",-((fft-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
    }
  }
  for (i=0;i<=fft/2-1;i++){
    if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){
      fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
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
    }
  }
  for (i=0;i<=fft/2-1;i++){
    if (-((-i)/non->deltat/c_v/fft-shift1)>non->min1 && -((-i)/non->deltat/c_v/fft-shift1)<non->max1){ 
      fprintf(outone,"%f %e %e\n",-((-i)/non->deltat/c_v/fft-shift1),fftOut[i][1],fftOut[i][0]);
    }
  }
    
  fclose(outone);
}
