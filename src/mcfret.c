#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"

/* Distribute the MCFRET calculations */
void mcfret(t_non *non){
    /* Define variables and arrays */
   
    /* Response functions for emission and absorption: real and imaginary part*/
    float *re_S_Abs,*im_S_Abs;
    float *re_S_Emi,*im_S_Emi;
    /*Vjesno: Question: Compared to DOS.c script we do not need veci and vecr because we do not want to project on the states?!*/
    /*Hamiltonian of the whole system - all donors and acceptors included*/
    float *Hamil_i_e;

/* Floats */
    float shift1;
  /* 1D Fourier transform */
    fftw_complex *fftIn,*fftOut;
    /*was is das?*/
    fftw_plan fftPlan;

/*File handles: Do we need the same ones like in DOS.c*/

    FILE *H_traj,*mu_traj;
    FILE *C_traj;
    FILE *outone,*log;
    FILE *Cfile;
  /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int x,ti,tj,i,j;
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
    /*Allocate memory for all the variables*/
    /*Allocate memory for the response functions*/
    nn2=non->singles*(non->singles+1)/2;
    re_S_Abs=(float *)calloc(nn2*nn2*non->tmax,sizeof(float));
    im_S_Abs=(float *)calloc(nn2*nn2*non->tmax,sizeof(float));
    re_S_Emi=(float *)calloc(nn2*nn2*non->tmax,sizeof(float));
    im_S_Emi=(float *)calloc(nn2*nn2*non->tmax,sizeof(float));
    /*Allocate memory for the Hamiltonian matrix*/
    Hamil_i_e=(float *)calloc(nn2,sizeof(float));

     /* Open Trajectory files */
    H_traj=fopen(non->energyFName,"rb");
    if (H_traj==NULL){
        printf("Hamiltonian file not found!\n");
        exit(1);
    }
      /*Here we want to call the routine for checking the trajectory files*/ 
    control(non);

    itime=0;
  /* Do calculation*/
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

    /* Read coupling */
    /*Question: Vjesno - Can we already do something here to nullify certain couplings??*/
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


  
    /* Call the MCFRET Routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Autodetect")) || (!strcmp(non->technique, "MCFRET-Absorption"))
        || (!strcmp(non->technique, "MCFRET-Emission")) || (!strcmp(non->technique, "MCFRET-Coupling")) || (!strcmp(non->technique,
        "MCFRET-Rate")) || (!strcmp(non->technique, "MCFRET-Analyse")) ) {
        printf("Perfroming the MCFRET calculation.\n");
    }

    /* Call the MCFRET Response for Absorption matrix */
    if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Absorption")  ){
        printf("Calculation Absortion Matrix for MCFRET.\n");
        mcfret_response_function(non,re_S_Abs,im_S_Abs,0);
        
    }

    /* Call the MCFRET Response for Absorption matrix */
    if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Emission")  ){
        printf("Calculation Emssion Matrix for MCFRET.\n");
        mcfret_response_function(non,re_S_Emi,im_S_Emi,1);
    }

    free(re_S_Abs),free(im_S_Abs);
    free(re_S_Emi),free(im_S_Emi);
    return;
}

void mcfret_autodetect(t_non *non, float treshold);

void mcfret_response_function(t_non *non, float *re_S_1, float *im_S_1,int emission){

    return;
}

void mcfret_coupling(t_non *non);

void mcfret_rate(t_non *non);

void mcfret_validate(t_non *non);
void mcfret_analyse(t_non *non);
