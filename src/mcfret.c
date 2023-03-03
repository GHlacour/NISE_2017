#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"

/* Main MCFRET routine only calling and combining the other subroutines */ 
void mcfret(t_non *non){
    int nn2;
    int segments;
  /* Response functions for emission and absorption: real and imaginary part*/
    float *re_Abs,*im_Abs;
    float *re_Emi,*im_Emi;
    float *J;
    float *rate_matrix;

  /*Allocate memory for the response functions*/
    nn2=non->singles*non->singles;
    re_Abs=(float *)calloc(nn2*non->tmax,sizeof(float));
    im_Abs=(float *)calloc(nn2*non->tmax,sizeof(float));
    re_Emi=(float *)calloc(nn2*non->tmax,sizeof(float));
    im_Emi=(float *)calloc(nn2*non->tmax,sizeof(float));
    J=(float *)calloc(nn2,sizeof(float));

  /* The rate matrix is determined by the integral over t1 for */
  /* Tr [ J * Abs(t1) * J * Emi(t1) ] */

    segments=project_dim(non);
    if (segments<2){
      printf(RED "Too few segments defined for MCFRET calculation!" RESET);
      exit(0);
    }
    rate_matrix=(float *)calloc(segments*segments,sizeof(float));

  /* Tell the user that we are in the MCFRET Routine */
     if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Autodetect")) || (!strcmp(non->technique, "MCFRET-Absorption"))
      || (!strcmp(non->technique, "MCFRET-Emission")) || (!strcmp(non->technique, "MCFRET-Coupling")) || (!strcmp(non->technique,
      "MCFRET-Rate")) || (!strcmp(non->technique, "MCFRET-Analyse")) ) {
          printf("Perfroming the MCFRET calculation.\n");
    }

/* Call the absorption routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Absorption"))){
        mcfret_response_function(re_Abs,im_Abs,non,0);
    }

/* Call the emission routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Emission"))){
        mcfret_response_function(re_Emi,im_Emi,non,1);
    }

/* Call the coupling routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Coupling"))){
        mcfret_coupling(J,non);
    }

/* Call the rate routine routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Rate"))){
        if ((!strcmp(non->technique, "MCFRET-Rate"))){
            /* Read in absorption, emission and coupling from file if needed */
        }
        mcfret_rate(rate_matrix,segments,re_Abs,im_Abs,re_Emi,im_Emi,J,non);
    }

/* NOTE!!! We still need to write the calculated stuff to file */
    write_matrix_to_file("RateMatrix.dat",rate_matrix,segments);

    free(re_Abs);
    free(im_Abs);
    free(re_Emi);
    free(im_Emi);
    free(J);
    free(rate_matrix);
    return;
}

/* Calculate Absorption/Emission matrix (depending on emission variable 0/1) */
void mcfret_response_function(float *re_S_1,float *im_S_1,t_non *non,int emission){
/* Define variables and arrays */
   /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int x,ti,tj,i,j;
    int t1;
    int elements;
    int cl,Ncl;

    /*Hamiltonian of the whole system - all donors and acceptors included*/
    float *Hamil_i_e;
    /*Vectors representing time dependent states: real and imaginary part*/
    float *vecr, *veci;
    float shift1; 

  /* Time parameters */
    time_t time_now,time_old,time_0;
  /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;

  /*File handles*/
    FILE *H_traj;
    FILE *C_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;

    /*Allocate memory for all the variables */
    /*Allocating memory for the real and imaginary part of the wave function that we need to propagate*/
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
    veci=(float *)calloc(non->singles*non->singles,sizeof(float));

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
      /* Initialize time-evolution operator */
        if (emission==0){
          	    /* Read Hamiltonian */
	          if (!strcmp(non->hamiltonian,"Coupling")){
                if (read_Dia(non,Hamil_i_e,H_traj,ti)!=1){
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }
            } else {
	              if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
	                  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	                  exit(1);
  	            }
            }
            /* Remove couplings between segments */
            multi_projection_Hamiltonian(Hamil_i_e,non);
            /* Use the thermal equilibrium as initial state */
            density_matrix(vecr,Hamil_i_e,non);
        } else { 
            unitmat(vecr,non->singles);
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
            } else {
	              if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	                  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
  	            }
            }
                /* Remove couplings between segments */
            multi_projection_Hamiltonian(Hamil_i_e,non);
                
                /* Update the MCFRET Response  */
            mcfret_response_function_sub(re_S_1, im_S_1,t1,non,vecr,veci);        
        
            #pragma omp parallel for
	          for (i=0; i<non->singles; i++){
	              j=i*non->singles;
	              /* Propagate vector*/
	              if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr+j,veci+j,non->ts,1);
	                  if (non->propagation==0){
	                      if (non->thres==0 || non->thres>1){
	                          propagate_vec_DIA(non,Hamil_i_e,vecr+j,veci+j,1);
	                      } else {
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

    /* Normalize response */
    for (t1=0;t1<non->tmax1*nn2;t1+=1){
        re_S_1[t1]=re_S_1[t1]/samples;
        im_S_1[t1]=im_S_1[t1]/samples;
    }
  /* Save time domain response */
  /*
    absorption_matrix=fopen("TD_absorption_matrix.dat","w");
    emission_matrix=fopen("TD_emission_matrix.dat","w");
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
        fprintf(absorption_matrix,"%f %e %e\n",t1*non->deltat,re_S_Abs[t1]/samples,im_S_Abs[t1]/samples);
        fprintf(emission_matrix, "%f %e %e\n",t1*non->deltat,re_S_Emi[t1]/samples,im_S_Emi[t1]/samples);
    }
    fclose(absorption_matrix);
    fclose(emission_matrix);
    */
    
    /*Free the memory*/
    free(vecr);	
    free(veci);  
}

void mcfret_response_function_sub(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci){
  int i,k;
  int N,nn2;
  N=non->singles;
  nn2=N*N;

  /* Update response matrix */
  for (i=0;i<non->singles;i++){
      for (k=0; k<non->singles; k++){
        /* We store response function so we can do matrix multiplication */
          re_S_1[t1*nn2+(N*i+k)]+=cr[N*i+k]; 
          im_S_1[t1*nn2+(N*i+k)]+=ci[N*i+k];
      }
  }

  return;
}

/* Find the average couplings but only between different segments */
void mcfret_coupling(float *J,t_non *non){
  /* Define variables and arrays */
   /* Integers */
    int nn2;
    int N_samples;
    int samples;
    int ti,i,j;
    int cl,Ncl;

    /*Hamiltonian of the whole system - all donors and acceptors included*/
    float *Hamil_i_e;
    float shift1; 

  /* Time parameters */
    time_t time_now,time_old,time_0;
  /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;

  /*File handles*/
    FILE *H_traj;
    FILE *C_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;

 /* Open Trajectory files */
    H_traj=fopen(non->energyFName,"rb");
    if (H_traj==NULL){
        printf("Hamiltonian file not found!\n");
        exit(1);
    }
      /*Here we want to call the routine for checking the trajectory files*/ 
    control(non);

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
      /* Initialize time-evolution operator */
          	    /* Read Hamiltonian */
	          if (!strcmp(non->hamiltonian,"Coupling")){
                if (read_Dia(non,Hamil_i_e,H_traj,ti)!=1){
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }
            } else {
	              if (read_He(non,Hamil_i_e,H_traj,ti)!=1){
	                  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
	                  exit(1);
  	            }
            }
            /* Extract couplings between segments */
            multi_projection_Coupling(Hamil_i_e,non);
            for (i=0;i<non->singles;i++){
              for (j=0;j<non->singles;j++){
                J[non->singles*i+j]=+Hamil_i_e[Sindex(i,j,non->singles)];
                J[non->singles*j+i]=+Hamil_i_e[Sindex(i,j,non->singles)];
              }
            }
        } /*We are closing the cluster loop*/

         /*Create NISE log file: Why not immediately in the begining?*/ 
        log=fopen("NISE.log","a");
        fprintf(log,"Finished sample %d\n",samples);
          
        time_now=log_time(time_now,log);
        fclose(log);
    }/*Closing the loop over samples*/
    
    /* Divide with total number of samples */
    for (i=0;i<non->singles;i++){
        for (j=0;j<non->singles;j++){
            J[non->singles*i+j]=J[non->singles*i+j]/samples;
        }
    }

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
    return;
}


/* Find MCFRET segments using an automatic scheme */
void mcfret_autodetect(t_non *non, float treshold);

/* Calculate actual rate matrix */
void mcfret_rate(float *rate_matrix,int segments,float *re_Abs,float *im_Abs,
    float *re_Emi,float *im_Emi,float *J,t_non *non){
    int nn2,N;
    int si,sj;
    int i,j,k;
    int *ns; /* Segment dimensions */
    int t1;
    float *rate_response;
    float *re_aux_mat,*im_aux_mat;
    float *re_aux_mat2,*im_aux_mat2;
    float *Zeros;
    N=non->singles;
    nn2=non->singles*non->singles;

    rate_response=(float *)calloc(non->tmax,sizeof(float));
    re_aux_mat=(float *)calloc(nn2,sizeof(float));
    im_aux_mat=(float *)calloc(nn2,sizeof(float));
    re_aux_mat2=(float *)calloc(nn2,sizeof(float));
    im_aux_mat2=(float *)calloc(nn2,sizeof(float));
    Zeros=(float *)calloc(nn2,sizeof(float));

    /* Do one rate at a time - so first we loop over segments */
    /* Tr [ J * Abs(t1) * J * Emi(t1) ] */
    for (si=0;si<segments;si++){
      for (sj=0;sj<segments;sj++){
        /* Exclude rate between same segments */
        if (sj!=si){
        /* Loop over time delay */
            for (t1=0;t1<non->tmax;t1++){
            /* Matrix multiplication - J Emi */
                segment_matrix_mul(J,Zeros,re_Emi+nn2*t1,im_Emi+nn2*t1,
                    re_aux_mat,im_aux_mat,non->psites,segments,si,sj,sj,N);
            /* Matrix multiplication - Abs (J Emi) */
                segment_matrix_mul(re_Abs+nn2*t1,im_Abs+nn2*t1,re_aux_mat,im_aux_mat,
                    re_aux_mat2,im_aux_mat2,non->psites,segments,si,si,sj,N);
            /* Matrix multiplication - J (Abs J Emi) */
                segment_matrix_mul(J,Zeros,re_aux_mat2,im_aux_mat2,
                    re_aux_mat,im_aux_mat,non->psites,segments,sj,si,sj,N);
            /* Here trace should be */
                rate_response[t1]=trace_rate(re_aux_mat,N);
            }
            /* Update rate matrix */
            rate_matrix[si+segments+sj]=integrate_rate_response(rate_response,non->tmax);
        }

      }
    }

    free(rate_response);
    free(ns);
    return;
}

/* Check if mcfret rates are in the incoherent limit */
void mcfret_validate(t_non *non);

/* Analyse rate matrix */
void mcfret_analyse(t_non *non);

/*This function will create a density matrix where every term is weighted with a Boltzmann weight*/
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non){
  int index,N;
  float *H,*e,*c2;
  float *cnr;
  /*float *crr;*/
  N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  c2=(float *)calloc(N,sizeof(float));
  cnr=(float *)calloc(N*N,sizeof(float));

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
  /* Find the inverse of the partition function */
  iQ=1.0/Q;

  /*Transform back to site basis*/ 
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
  free(e);
  free(cnr);
  return;
}

/* Matrix multiplication for different segments */
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,
    float *rC,float *iC,int *psites,int segments,int si,int sj,int sk,int N){
    int i,j,k;
    /* Set initial values of results matrix to zero to be sure */
    clearvec(rC,N*N);
    clearvec(iC,N*N);
    for (i=0;i<N;i++){
        if (psites[i]==si){
            for (j=0;i<N;i++){
                if (psites[j]==sj){
                    for (k=0;i<N;i++){
                        if (psites[k]==sk){
                            rC[i*N+j]+=rA[i*N+k]*rB[k*N+j]-iA[i*N+k]*iB[k*N+j];
                            iC[i*N+j]+=rA[i*N+k]*iB[k*N+j]+iA[i*N+k]*rB[k*N+j];
                        } 
                    } 
                }
            }
        }
    }
    return;
} 

/* Find the trace of the matrix */
float trace_rate(float *matrix,int N){
  int i;
  float trace;
  for (i=0;i<N;i++){
    trace=trace+matrix[N*i+i];
  }
  return trace;
}

/* Integrate the rate response */
float integrate_rate_response(float *rate_response,int T){
    int i;
    float simple; /* Variable for naieve box integral */
    float simp13; /* Variable for Simpsons 1/3 rule integral */
    for (i=0;i<T;i++){
        simple+=rate_response[i];
        if (i%2==0){
          simp13+=rate_response[i]/3;
        } else {
          simp13+=4*rate_response[i]/3;
        }
    }
    if (abs(simple-simp13)/abs(simp13)>0.1){
      printf(YELLOW "Warning the timesteps may be to large for integration!\n" RESET);
    }
    return simp13;
}

/* Write a square matrix to a text file */
float write_matrix_to_file(char fname[],float *matrix,int N){
  FILE *file_handle;
  int i,j;
  file_handle=fopen(fname,"w");
  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      fprintf(file_handle,"%f ",matrix[i*N+j]);
    }
    fprintf(file_handle,"\n");
  }
}