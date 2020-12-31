#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "calc_Diffusion.h"
#include "1DFFT.h"

void calc_Diffusion(t_non *non){
  // Initialize variables
  float *re_S_1,*im_S_1; // The first-order response function
  float *Hamil_i_e;
  // Aid arrays
  float *vecr,*veci;
  float *Pop,*pos_i,*pos_f,*Ori,*Anis;
  float dist[3];

  /* Floats */
  float shift1;
  float box_size;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile,*P_traj;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i,j;
  int t1,fft;
  int elements;
  int cl,Ncl;
  int a,b;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  // Allocate memory
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  Pop=(float *)calloc(non->tmax,sizeof(float));
  pos_i=(float *)calloc(non->singles*3,sizeof(float));
  pos_f=(float *)calloc(non->singles*3,sizeof(float));
  Ori=(float *)calloc(non->tmax,sizeof(float));
  Anis=(float *)calloc(non->tmax,sizeof(float));

  /* Open Trajectory files */
  H_traj=fopen(non->energyFName,"rb");
  if (H_traj==NULL){
    printf("Hamiltonian file not found!\n");
    exit(1);
  }
  P_traj=fopen("Position.bin","rb");
  if (P_traj==NULL){
    printf("Position file %s not found!\n","Position.bin");
    exit(1);
  }
  // Read box size
  fread(&box_size,sizeof(float),1,P_traj);

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

  vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles*non->singles,sizeof(float));

  /* Read coupling */
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

    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
        printf("Cluster trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d\n",ti);
        exit(1);
      }
      if (non->cluster==cl){
        Ncl++;
      }
    }
    if (non->cluster==-1 || non->cluster==cl){
      /* Initialize */
      for (a=0;a<non->singles;a++) vecr[a+a*non->singles]=1.0;
      
      /* Read initial coordinates */
      fseek(P_traj,sizeof(float)*(1+ti*non->singles*3),SEEK_SET);
      fread(pos_i,sizeof(float),non->singles*3,P_traj);
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

        /* Read final coordinates */
        fseek(P_traj,sizeof(float)*(1+tj*non->singles*3),SEEK_SET);
        fread(pos_f,sizeof(float),non->singles*3,P_traj);
        /* Calculate distance evolution of wave */
        for (a=0;a<non->singles;a++){
          for (b=0;b<non->singles;b++){
            Pop[t1]+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*distance(pos_f,pos_i,b,a,non->singles,box_size);
            Pop[t1]+=veci[a+b*non->singles]*veci[a+b*non->singles]*distance(pos_f,pos_i,b,a,non->singles,box_size);
          }
        }

        /* Calculate distance evolution of center */
        for (a=0;a<non->singles;a++){
          dist[0]=0,dist[1]=0,dist[2]=0;
          for (x=0;x<3;x++){
            for (b=0;b<non->singles;b++){
              dist[x]+=vecr[a+b*non->singles]*vecr[a+b*non->singles]*distance_x(pos_f,pos_i,b,a,non->singles,box_size,x);
              dist[x]+=veci[a+b*non->singles]*veci[a+b*non->singles]*distance_x(pos_f,pos_i,b,a,non->singles,box_size,x);
            }
          }
          Ori[t1]+=dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
        }
        for (a=0;a<non->singles;a++){
          /* Propagate vector */
          if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,non->ts,1);
          if (non->propagation==0){
            if (non->thres==0 || non->thres>1){
              propagate_vec_DIA(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
            } else {
              elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+a*non->singles,veci+a*non->singles,1);
            if (samples==0 && a==0){
              if (t1==0){
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

    // Log time
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);
          
    time_now=log_time(time_now,log);
    fclose(log);
  }

  free(vecr);
  free(veci);
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

  fclose(H_traj);
  if (non->cluster!=-1){
    fclose(Cfile);
  }

  /* Save time domain response */
  outone=fopen("Dif.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f %e %e\n",t1*non->deltat,Pop[t1]/samples/non->singles,Ori[t1]/samples/non->singles);
  }
  fclose(outone);

  printf("----------------------------------------------\n");
  printf(" Diffusion calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

