#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "analyse.h"

void analyse(t_non *non){
  // Initialize variables
  float *average_frequency;
  float *average_coupling;
  float *average_H;
  float avall,flucall;
  float *fluctuation;
  float *Jfluctuation;
  float *mu_eg,*Hamil_i_e,*H,*e;
  float participation_ratio;
  float *cEig,*dip2,*cDOS;

  // Aid arrays
  float *vecr,*veci,*vecr_old,*veci_old;

  /* Floats */
  float shift1;
  float x;

  /* File handles */
  FILE *H_traj,*mu_traj,*PDF_traj;
  FILE *outone,*log;
  FILE *Cfile;

  /* Integers */
  int nn2,N;
  int itime,N_samples;
  int samples;
  int ti,tj,i,j;
  int t1,fft;
  int elements;
  int Nsam;
  int counts;
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
  N=non->singles;
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  average_H=(float *)calloc(nn2,sizeof(float));
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  average_frequency=(float *)calloc(N,sizeof(float));
  average_coupling=(float *)calloc(N,sizeof(float)); // Coupling strength = sum of all couplings for one molecule
  fluctuation=(float *)calloc(N,sizeof(float));
  Jfluctuation=(float *)calloc(N,sizeof(float));
  cEig=(float *)calloc(N,sizeof(float));
  cDOS=(float *)calloc(N,sizeof(float));
  dip2=(float *)calloc(N,sizeof(float));

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


  itime=0;
  // Do calculation
  N_samples=(non->length-1)/non->sample+1;
  if (N_samples>0) {
    printf("Making %d samples!\n",N_samples);
  } else {
    printf("Insufficient data to calculate spectrum.\n");
    printf("Please, lower max times or provide longer\n");
    printf("trajectory.\n");
    exit(1);
  }

  if (non->end==0) non->end=N_samples;
  Nsam=non->end-non->begin;

  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);


  participation_ratio=0;
  avall=0;
  flucall=0;
  counts=0;
  Ncl=0;

  // Loop over samples first time
  for (samples=non->begin;samples<non->end;samples++){

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
    // Read Hamiltonian
    if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
      exit(1);
    }
    build_diag_H(Hamil_i_e,H,e,N);
    participation_ratio+=calc_participation_ratio(N,H);
    find_dipole_mag(non,dip2,samples,mu_traj,H);
    counts=find_cEig(cEig,cDOS,dip2,H,e,N,non->min1,non->max1,counts,non->shifte);
    // Find Averages
    for (i=0;i<non->singles;i++){
      average_frequency[i]+=Hamil_i_e[Sindex(i,i,N)];
      avall+=Hamil_i_e[Sindex(i,i,N)];
      for (j=0;j<non->singles;j++){
        if (j>=i){
          average_H[Sindex(i,j,N)]+=Hamil_i_e[Sindex(i,j,N)];
        }
        if (j!=i){
          average_coupling[i]+=Hamil_i_e[Sindex(i,j,N)];
        }
      }
    }     
  }
  }
  if (Ncl>0) Nsam=Ncl;
  // Normalize average_frequencies
  for (i=0;i<non->singles;i++){
    average_frequency[i]/=Nsam;
    average_coupling[i]/=Nsam;
    for(j=i;j<non->singles;j++){
      average_H[Sindex(i,j,non->singles)]/=Nsam;
    }
  }
  avall/=(Nsam*non->singles);   

  // Loop over samples second time
  for (samples=non->begin;samples<non->end;samples++){
    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
        printf("Cluster trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d\n",ti);
        exit(1);
      }
      //      printf("%d\n",cl);
      // Configuration belong to cluster
    }
    if (non->cluster==-1 || non->cluster==cl){

    // Read Hamiltonian
    if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
      exit(1);
    }
    // Find standard deviation for frequencies
    for (i=0;i<non->singles;i++){
      x=(Hamil_i_e[Sindex(i,i,N)]-average_frequency[i]);
      fluctuation[i]+=x*x;
      x=(Hamil_i_e[Sindex(i,i,N)]-avall);
      flucall+=x*x;
      x=0;
      for (j=0;j<non->singles;j++){
        if (i!=j) x+=Hamil_i_e[Sindex(i,j,N)];
      }
      x=x-average_coupling[i];
      Jfluctuation[i]+=x*x;
    }
  }
  }
  // Normalize fluctuations and take square root
  for (i=0;i<non->singles;i++){
    fluctuation[i]/=Nsam;
    fluctuation[i]=sqrt(fluctuation[i]);
    Jfluctuation[i]/=Nsam;
    Jfluctuation[i]=sqrt(Jfluctuation[i]);
  }
  flucall/=(Nsam*non->singles);   
  flucall=sqrt(flucall);

  // Write Average Hamiltonian in GROASC format
  outone=fopen("Av_Hamiltonian.txt","w");
    if (outone==NULL){
    printf("Problem encountered opening Analyse.dat for writing.\n");
    printf("Disk full or write protected?\n");
    exit(1);
  }
  fprintf(outone,"0 ");
  for (i=0;i<non->singles;i++){
    for (j=i;j<non->singles;j++){
      if (i==j) average_H[Sindex(i,j,non->singles)]+=non->shifte;
      fprintf(outone,"%f ",average_H[Sindex(i,j,non->singles)]);
    }
  }
  fprintf(outone,"\n");
  fclose(outone);

  outone=fopen("Analyse.dat","w");
  if (outone==NULL){
    printf("Problem encountered opening Analyse.dat for writing.\n");
    printf("Disk full or write protected?\n");
    exit(1);
  }

  participation_ratio/=(non->singles*Nsam);
  printf("===================================\n");
  printf("Result of Hamiltonian analysis:\n");
  printf("Delocalization size according to\n");
  printf("Thouless, Phys. Rep. 13:93 (1974)\n");
  printf("R=%f\n",participation_ratio);
  printf("Average site frequency %f cm-1.\n",avall+non->shifte);
  printf("Overall standard deviation of site\n");
  printf("frequencies from averall average:\n");
  printf("%f cm-1.\n",flucall);
  printf("===================================\n");

  // Print output to file
  fprintf(outone,"# Using frequency range %f to %f cm-1\n",non->min1,non->max1);
  fprintf(outone,"# Site AvFreq. SDFreq AvJ SDJ cEig cDOS\n");
  for (i=0;i<non->singles;i++){
    fprintf(outone,"%d %f %f %f %F %f %f\n",i,average_frequency[i]+non->shifte,fluctuation[i],average_coupling[i],Jfluctuation[i],cEig[i]/counts,cDOS[i]/counts);
  }

  free(average_frequency);
  free(fluctuation);
  free(average_coupling);
  free(Jfluctuation);
  // free(mu_eg);
  free(Hamil_i_e);
  free(average_H);
  free(cEig);
  free(cDOS);
  free(H);
  free(e);
  free(dip2);
  
  // The calculation is finished, lets write output
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  fclose(mu_traj),fclose(H_traj);
  fclose(outone);
  if (non->cluster!=-1){
    fclose(Cfile);
  } 
 

  printf("----------------------------------------------\n");
  printf(" Analyse calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

float calc_participation_ratio(int N,float *H){
  int i,j,a,b;
  float inter,parti;

  parti=0;
  for (i=0;i<N;i++){
    inter=0;
    // Loop over sites
    for (j=0;j<N;j++){
      inter+=H[i+N*j]*H[i+N*j]*H[i+N*j]*H[i+N*j];
    }
    parti+=1.0/inter;
  }

  return parti;
}

int find_cEig(float *cEig,float *cDOS,float *dip2,float *H,float *e,int N,float min,float max,int counts,float shift){
  int i,j;
  
  // Loop over eigenstates
  for (i=0;i<N;i++){
    // Only for eigenstates in given range
    if (e[i]>min-shift && e[i]<max-shift){
      counts++;
      // Loop over sites
      for (j=0;j<N;j++){
	cEig[j]+=dip2[i]*H[i+N*j]*H[i+N*j];
	cDOS[j]+=H[i+N*j]*H[i+N*j];
      }
    }
  }
  return counts;
}  

// Find dipole magnitude
void find_dipole_mag(t_non *non,float *dip2,int step,FILE *mu_traj,float *H){
  float *dip,*dipeb;
  int i,j,x,N;

  N=non->singles;
  dip=(float *)calloc(N,sizeof(float));
  dipeb=(float *)calloc(N,sizeof(float));

  for (i=0;i<N;i++){
    dip2[i]=0;
  }
  for (x=0;x<3;x++){
    // Read mu(ti)
    if (read_mue(non,dip,mu_traj,step,x)!=1){
      printf("Dipole trajectory file to short, could not fill buffer!!!\n");
      printf("ITIME %d %d\n",step,x);
      exit(1);
    }
    // Transform to eigen basis
    for (i=0;i<N;i++){
      dipeb[i]=0;
      for (j=0;j<N;j++){
	dipeb[i]+=H[i+j*N]*dip[j]; // i is eigen state, j site
      }
    }
    for (i=0;i<N;i++){
      dip2[i]+=dipeb[i]*dipeb[i];
    }
  }
  free(dip);
  free(dipeb);
  return;
}
