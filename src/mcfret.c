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




   
void mcfret(t_non *non){
  /* Call the MCFRET Routine */
   if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Autodetect")) || (!strcmp(non->technique, "MCFRET-Absorption"))
        || (!strcmp(non->technique, "MCFRET-Emission")) || (!strcmp(non->technique, "MCFRET-Coupling")) || (!strcmp(non->technique,
        "MCFRET-Rate")) || (!strcmp(non->technique, "MCFRET-Analyse")) ) {
        printf("Perfroming the MCFRET calculation.\n");
    }

/* Define variables and arrays */
   /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int x,ti,tj,i,j;
    int t1,fft;
    int elements;
    int cl,Ncl;

    /* Response functions for emission and absorption: real and imaginary part*/
    float *re_S_Abs,*im_S_Abs;
    float *re_S_Emi,*im_S_Emi;
    /*Hamiltonian of the whole system - all donors and acceptors included*/
    float *Hamil_i_e;
    /*Vectors representing time dependent states: real and imaginary part*/
    float *vecr, *veci;
    /*Density matrix*/
    float *dens_mat;

    float shift1;
  /* 1D Fourier transform */
    fftw_complex *fftIn,*fftOut;
    /*was is das?*/
    fftw_plan fftPlan;  

  /* Time parameters */
    time_t time_now,time_old,time_0;
  /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;


  /*File handles*/

    FILE *H_traj,*mu_traj;
    FILE *C_traj;
    FILE *absorption_matrix, *emission_matrix,*log;
    FILE *Cfile;

    /*Allocate memory for all the variables*/
    /*Allocate memory for the response functions*/
    nn2=non->singles*(non->singles+1)/2;
    re_S_Abs=(float *)calloc(nn2*non->tmax,sizeof(float));
    im_S_Abs=(float *)calloc(nn2*non->tmax,sizeof(float));
    re_S_Emi=(float *)calloc(nn2*non->tmax,sizeof(float));
    im_S_Emi=(float *)calloc(nn2*non->tmax,sizeof(float));

    /*Allocating memory for the real and imaginary part of the wave function that we need to propagate*/
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
    veci=(float *)calloc(non->singles*non->singles,sizeof(float));

    /*Allocating memory for the density matrix*/
    dens_mat=(float *)calloc(non->singles*non->singles,sizeof(float));

 /* Open Trajectory files */
    H_traj=fopen(non->energyFName,"rb");
    if (H_traj==NULL){
        printf("Hamiltonian file not found!\n");
        exit(1);
    }
      /*Here we want to call the routine for checking the trajectory files*/ 
    control(non);

    itime=0;

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

    
    /* Read coupling: Since in this case couplings are static, we want to read them in only once*/
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

    /*Looping over samples: Each sample represents a different starting point on the Hamiltonian trajectory*/
  
    for (samples=non->begin;samples<non->end;samples++){ 
        ti=samples*non->sample;
        if (non->cluster!=-1){
            if (read_cluster(non,ti,&cl,Cfile)!=1){
	              printf("Cluster trajectory file to short, could not fill buffer!!!\n");
	              printf("ITIME %d\n",ti);
	              exit(1);
            }
            /*Configuration belong to cluster*/ 
        if (non->cluster==cl){
	          Ncl++;
        }
        }
        if (non->cluster==-1 || non->cluster==cl){
      /*Initialize time-evolution operator*/ 
            unitmat(vecr,non->singles);
      /*Do projection on selected sites if asked*/
            if (non->Npsites>0){
                for (i=0;i<non->singles;i++){
                    j=i*non->singles;
                    projection(vecr+j,non);
                }
            }
            clearvec(veci,non->singles*non->singles);

      /*Loop over delay*/ 
            for (t1=0;t1<non->tmax;t1++){
	              tj=ti+t1;
	    /* Read Hamiltonian */
	              if (!strcmp(non->hamiltonian,"Coupling")){
                    if (read_Dia(non,Hamil_i_e,H_traj,tj)!=1){
                        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                        exit(1);
                    }
                }
                else {
	                  if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	                      printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	                      exit(1);
  	                }
                }
             /* Call the MCFRET Response for Absorption matrix */
                if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Absorption")  ){
                    printf("Calculation Absortion Matrix for MCFRET.\n");
                    mcfret_response_function(Hamil_i_e, re_S_Abs, im_S_Abs,tj,non,vecr,veci, dens_mat,0);
        
                }

            /* Call the MCFRET Response for Emission matrix */
                if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Emission")  ){
                    printf("Calculation Emssion Matrix for MCFRET.\n");
                    mcfret_response_function(Hamil_i_e, re_S_Emi, im_S_Emi,tj,non,vecr,veci, dens_mat,1);
                }
                #pragma omp parallel for
	              for (i=0; i<non->singles; i++){
	                  j=i*non->singles;
	              /* Propagate vector*/
	                  if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr+j,veci+j,non->ts,1);
	                      if (non->propagation==0){
	                          if (non->thres==0 || non->thres>1){
	                              propagate_vec_DIA(non,Hamil_i_e,vecr+j,veci+j,1);
	                          } 
                        else {
	                          elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+j,veci+j,1);
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
            
              }/*We are closing the loop over time delays -t1 times*/
         } /*We are closing the cluster loop*/

         /*Create NISE log file: Why not immediately in the begining?*/ 
         log=fopen("NISE.log","a");
         fprintf(log,"Finished sample %d\n",samples);
          
         time_now=log_time(time_now,log);
         fclose(log);
    }/*Closing the loop over samples*/
    
    /*The calculation is finished, lets write output*/
    log=fopen("NISE.log","a");
    fprintf(log,"Finished Calculating Response!\n");
    fprintf(log,"Writing to file!\n");  
    fclose(log);

    if (non->cluster!=-1){
        printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
        if (samples==0){ /*Avoid dividing by zero*/ 
            samples=1;
        }
    }

    fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }

  /* Save time domain response */
    absorption_matrix=fopen("TD_absorption_matrix.dat","w");
    emission_matrix=fopen("TD_emission_matrix.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
        fprintf(absorption_matrix,"%f %e %e\n",t1*non->deltat,re_S_Abs[t1]/samples,im_S_Abs[t1]/samples);
        fprintf(emission_matrix, "%f %e %e\n",t1*non->deltat,re_S_Emi[t1]/samples,im_S_Emi[t1]/samples);
    }
    fclose(absorption_matrix);
    fclose(emission_matrix);
    
    
    /*Free the memory*/
    
    free(re_S_Abs);
    free(im_S_Abs);
    free(re_S_Emi);
    free(im_S_Emi);

    
    free(vecr);	
    free(veci);

    free(dens_mat);

    
}
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non){
  /*This function will create a density matrix where every term is weighted with a Boltzmann weight*/
  int index,N;
  float *H,*e,*c2,*c1;
  float *cnr;
  /*float *crr;*/
  N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  c2=(float *)calloc(N,sizeof(float));
  /*c1=(float *)calloc(N,sizeof(float));*/
  cnr=(float *)calloc(N*N,sizeof(float));
  /*crr=(float *)calloc(N*N,sizeof(float));*/
  int a,b,c;
  float kBT=non->temperature*0.695; // Kelvin to cm-1
  float Q,iQ;

  /*Build Hamiltonian*/
  for (a=0;a<N;a++){
      H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; /*Diagonal*/
      for (b=a+1;b<N;b++){
          H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
          H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      }
  }
  diagonalizeLPD(H,e,N);
 
  /*Exponentiate [U=exp(-H/kBT)]*/
  for (a=0;a<N;a++){
      c2[a]=exp(-e[a]/kBT);
      Q=Q+c2[a];
  }
  iQ=1.0/Q;

  /*Transform to site basis*/ 
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          cnr[b+a*N]+=H[b+a*N]*c2[b]*iQ;
      }
  }  
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          for (c=0;c<N;c++){
              density_matrix[a+c*N]+=H[b+a*N]*cnr[b+c*N];
          }
      }
  }
  
  free(H);
  free(c2);
  free(c1);
  free(e);
  /*free(crr);*/
  free(cnr);
  return;
}

void mcfret_autodetect(t_non *non, float treshold);

void mcfret_response_function(float* Hamiltonian, float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci, float *density_matrix, int emission){
  int i,j, k;
  int N,nn2;
  N=non->singles;
  nn2=N*(N-1)/2;
  if (emission==1){
      printf("Calculating emission matrix");
      /*density_matrix(density_matrix,Hamiltonian,non);*/
      /*Calculate emission matrix*/
      for (i=0;i<non->singles;i++){
          j=i+i*non->singles;
          for (k=i+1; k<non->singles; k++){
            re_S_1[t1*nn2+(N*i+i*k)]=cr[j]*density_matrix[N*i+i*k];
            /*Minus sign to take care of complex conjugate*/
            im_S_1[t1*nn2+(N*i+i*k)]=-ci[j]*density_matrix[N*i+i*k];
          }
      }
  }
      else {
        unitmat(cr,non->singles);
        printf("Calculating absorption matrix");
        for (i=0;i<non->singles;i++){
            j=i+i*non->singles;
            for (k=i+1; k<non->singles; k++){
              /*We do not have to weight this with density matrix*/
                re_S_1[t1*nn2+(N*i+i*k)]=cr[i+N*j]*density_matrix[N*i+i*k];
                im_S_1[t1*nn2+(N*i+i*k)]=ci[i+N*j]*density_matrix[N*i+i*k];
            }
        
        }
      }

return;
}

void mcfret_coupling(t_non *non){


return;
}
void mcfret_rate(t_non *non);

void mcfret_validate(t_non *non);
void mcfret_analyse(t_non *non);
