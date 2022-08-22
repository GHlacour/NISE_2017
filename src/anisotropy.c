#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "anisotropy.h"

void anisotropy(t_non *non){
  // Initialize variables
  float avall,flucall;
  float *Hamil_i_e,*H,*e,*mu_eg;

  // Aid arrays
  float *vecr,*veci,*vecr_old,*veci_old;
  float *pos_i,*pos_f;
  float *Anis,*Ori;

  /* Floats */
  float shift1;
  float pr,pi;
  float norm,sum,sum2;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *outone,*log;

  /* Integers */
  int nn2,N;
  int itime,N_samples;
  int samples;
  int ti,tj,i,j;
  int t1,fft;
  int elements;
  int Nsam;
  int counts;
  int a,b,c,d;
  int x;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  // Allocate memory
  N=non->singles;
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  mu_eg=(float *)calloc(N,sizeof(float));
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  Anis=(float *)calloc(non->tmax,sizeof(float));
  Ori=(float *)calloc(non->tmax,sizeof(float));
  pos_i=(float *)calloc(non->singles*3,sizeof(float));
  pos_f=(float *)calloc(non->singles*3,sizeof(float));

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
  Nsam=non->end-non->begin;

  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);

  if (non->propagation==1){
    printf("==============================================================\n");
    printf(RED "You cannot use the Coupling propagation scheme for Anisotropy\n");
    printf("calculations! Please, use the 'Propagation Sparse' option.\n");
    printf("Use 'Threshold 0.0' for highest accuracy.\n");
    printf("Aborting run now.\n" RESET);
    printf("==============================================================\n");
    exit(0);
  }

  /* Loop over samples */
  for (samples=non->begin;samples<non->end;samples++){
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
    veci=(float *)calloc(non->singles*non->singles,sizeof(float));
    /* Initialize */
    for (a=0;a<non->singles;a++) vecr[a+a*non->singles]=1.0;

    ti=samples*non->sample;      
    for (x=0;x<3;x++){
      // Read mu(ti)
      if (read_mue(non,pos_i+x*non->singles,mu_traj,ti,x)!=1){
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n",ti,x);
        exit(1);
      }
    }
    for (a=0;a<non->singles;a++){
      norm=0;
      for (x=0;x<3;x++) norm+=pos_i[a+x*non->singles]*pos_i[a+x*non->singles];
      norm=1./sqrt(norm);
      for (x=0;x<3;x++) pos_i[a+x*non->singles]*=norm;
    }
    for (t1=0;t1<non->tmax;t1++){
      tj=ti+t1;
      /* Read Hamiltonian */
      if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
        exit(1);
      }
      // Read mu(ti)
      for (x=0;x<3;x++){
         if (read_mue(non,pos_f+x*non->singles,mu_traj,tj,x)!=1){
            printf("Dipole trajectory file to short, could not fill buffer!!!\n"
);
            printf("ITIME %d %d\n",ti,x);
            exit(1);
         }
      }

      for (a=0;a<non->singles;a++){
         norm=0;
         for (x=0;x<3;x++) norm+=pos_f[a+x*non->singles]*pos_f[a+x*non->singles
];
         norm=1./sqrt(norm);
         for (x=0;x<3;x++) pos_f[a+x*non->singles]*=norm;
      }

      // Calculate anisotropy evolution
      for (a=0;a<non->singles;a++){
        sum=0;
        sum2=0;
        for (x=0;x<3;x++){
          norm=0;
          for (b=0;b<non->singles;b++){
            norm+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*pos_i[a+x*non->singles]*pos_f[b+x*non->singles];
            norm+=veci[a+b*non->singles]*veci[a+b*non->singles]*pos_i[a+x*non->singles]*pos_f[b+x*non->singles];
          }
          sum+=norm;
          sum2+=pos_i[a+x*non->singles]*pos_f[a+x*non->singles];
        }
        Anis[t1]+=sum*sum;
        Ori[t1]+=sum2*sum2;
      }

      build_diag_H(Hamil_i_e,H,e,non->singles);

      
      for (a=0;a<non->singles;a++){
        /* Probagate vector */
        if (non->thres==0 || non->thres>1){
          propagate_vec_DIA(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
        } else {
          elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
          if (samples==non->begin && a==0){
            if (t1==0){
              printf("Sparce matrix efficiency: %f pct.\n",(1-(1.0*elements/(non->singles*non->singles)))*100);
              printf("Pressent tuncation %f.\n",non->thres/(non->deltat*icm2ifs*twoPi/non->ts)*(non->deltat*icm2ifs*twoPi/non->ts));
              printf("Suggested truncation %f.\n",0.001);
            }
          }
        }
      }
    }

  }
  /* Correct for when not starting at sample zero */
  samples=samples-non->begin;
  /* Write anisotropy */
  outone=fopen("Ani.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
     Anis[t1]=0.2*(3*Anis[t1]/samples/non->singles-1);
     Ori[t1]=0.2*(3*Ori[t1]/samples/non->singles-1);
     fprintf(outone,"%f %e %e\n",t1*non->deltat,Anis[t1],Ori[t1]);
  }
  fclose(outone);

  free(Hamil_i_e);
  free(H);
  free(e);
  free(Anis);
  free(Ori);
  free(pos_i);
  free(pos_f);
  // The calculation is finished, lets write output
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Population Transfer!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  fclose(H_traj);
  fclose(mu_traj); 

  printf("----------------------------------------------\n");
  printf(" Anisotropy calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

