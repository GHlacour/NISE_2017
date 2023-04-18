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
#include "propagate.h"
#include "read_trajectory.h"

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
          printf("Performing the MCFRET calculation.\n");
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
    int segments;

    /*Hamiltonian of the whole system - all donors and acceptors included*/
    float *Hamil_i_e;
    /*Vectors representing time dependent states: real and imaginary part*/
    float *vecr, *veci;
    /* Transition dipoles for coupling on the fly */
    float *mu_xyz;
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
    FILE *mu_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;
    FILE *absorption_matrix; 

    /*Allocate memory for all the variables */
    /*Allocating memory for the real and imaginary part of the wave function that we need to propagate*/
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
    veci=(float *)calloc(non->singles*non->singles,sizeof(float));
    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));
    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /*Here we want to call the routine for checking the trajectory files*/ 
    control(non);

    itime=0;

    /* Initialize sample numbers */
    N_samples=determine_samples(non);
    segments=project_dim(non);
    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);

    /* Read coupling, this is done if the coupling and transition-dipoles are */
    /* time-independent and only one snapshot is stored */
    read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

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
            if (emission==1){
                /* Read Hamiltonian */
                read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
	          
                /* Remove couplings between segments */
                multi_projection_Hamiltonian(Hamil_i_e,non);

                /* Use the thermal equilibrium as initial state */
                density_matrix(vecr,Hamil_i_e,non,segments);
		   // write_matrix_to_file("Density.dat",vecr,non->singles);
            } else { 
                unitmat(vecr,non->singles);
		//write_matrix_to_file("Unit.dat",vecr,non->singles);
            }
            clearvec(veci,non->singles*non->singles);
        
            /*Loop over delay*/ 
            for (t1=0;t1<non->tmax;t1++){
	              tj=ti+t1;
	              /* Read Hamiltonian */
                read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
	          
                /* Remove couplings between segments */
                multi_projection_Hamiltonian(Hamil_i_e,non);
                
                /* Update the MCFRET Response  */
                mcfret_response_function_sub(re_S_1, im_S_1,t1,non,vecr,veci);        
                if (emission==0){         
                   propagate_matrix(non,Hamil_i_e,vecr,veci,-1,samples,t1*x);
                } else {
		   propagate_matrix(non,Hamil_i_e,vecr,veci,1,samples,t1*x);
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
    for (t1=0;t1<non->tmax1*non->singles*non->singles;t1++){
        re_S_1[t1]=re_S_1[t1]/samples;
        im_S_1[t1]=im_S_1[t1]/samples;
    }
  /* Save time domain response */
   if (emission==0){ 
    absorption_matrix=fopen("TD_absorption_matrix.dat","w");
   } else {
    absorption_matrix=fopen("TD_emission_matrix.dat","w");
   }
   fprintf(absorption_matrix,"Samples %d\n",samples);
    for (t1=0;t1<non->tmax1;t1+=non->dt1){
        fprintf(absorption_matrix,"%f %e %e\n",t1*non->deltat,re_S_1[t1],im_S_1[t1]);
 //       fprintf(emission_matrix, "%f %e %e\n",t1*non->deltat,re_S_Emi[t1]/samples,im_S_Emi[t1]/samples);
    }
    fclose(absorption_matrix);
   // fclose(emission_matrix);
    
    
    /*Free the memory*/
    free(vecr);	
    free(veci);  
    free(Hamil_i_e);
    free(mu_xyz);
}

/* Sub routine for adding up the calculated response in the response function */
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
    float *mu_xyz;
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
    FILE *mu_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;

    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));

    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /*Here we want to call the routine for checking the trajectory files*/ 
    control(non);

 /* Initialize sample numbers */
  N_samples=determine_samples(non);
  Ncl=0;

  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);
  samples=1;


    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);


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
            /* Read Hamiltonian */
            read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
          	    
            /* Extract couplings between segments */
            multi_projection_Coupling(Hamil_i_e,non);
            for (i=0;i<non->singles;i++){
              for (j=i+1;j<non->singles;j++){
                J[non->singles*i+j]+=Hamil_i_e[Sindex(i,j,non->singles)];
                J[non->singles*j+i]+=Hamil_i_e[Sindex(i,j,non->singles)];
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

    fclose(mu_traj),fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }
    free(Hamil_i_e);
    free(mu_xyz);  
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
    float rate;
    float *re_aux_mat,*im_aux_mat;
    float *re_aux_mat2,*im_aux_mat2;
    float *Zeros;
    FILE *ratefile;
    N=non->singles;
    nn2=non->singles*non->singles;

    rate_response=(float *)calloc(non->tmax,sizeof(float));
    re_aux_mat=(float *)calloc(nn2,sizeof(float));
    im_aux_mat=(float *)calloc(nn2,sizeof(float));
    re_aux_mat2=(float *)calloc(nn2,sizeof(float));
    im_aux_mat2=(float *)calloc(nn2,sizeof(float));
    Zeros=(float *)calloc(nn2,sizeof(float));
  
    ratefile=fopen("RateFile.dat","w");
    /* Do one rate at a time - so first we loop over segments */
    /* Tr [ J * Abs(t1) * J * Emi(t1) ] */
    for (si=0;si<segments;si++){
      for (sj=0;sj<segments;sj++){
        /* Exclude rate between same segments */
        if (sj!=si){
        /* Loop over time delay */
//	    fprintf(ratefile,"%d %d\n",si,sj);
            for (t1=0;t1<non->tmax;t1++){
            /* Matrix multiplication - J Emi */
                segment_matrix_mul(J,Zeros,re_Emi+nn2*t1,im_Emi+nn2*t1,
                    re_aux_mat,im_aux_mat,non->psites,segments,si,sj,sj,N);
//		fprintf(ratefile,"JE %f %f %f %f\n",re_aux_mat[0],re_aux_mat[1],re_aux_mat[2],re_aux_mat[3]);
            /* Matrix multiplication - Abs (J Emi) */
                segment_matrix_mul(re_Abs+nn2*t1,im_Abs+nn2*t1,re_aux_mat,im_aux_mat,
                    re_aux_mat2,im_aux_mat2,non->psites,segments,si,si,sj,N);
//		fprintf(ratefile,"AJE %f %f %f %f\n",re_aux_mat2[0],re_aux_mat2[1],re_aux_mat2[2],re_aux_mat2[3]);
            /* Matrix multiplication - J (Abs J Emi) */
                segment_matrix_mul(J,Zeros,re_aux_mat2,im_aux_mat2,
                    re_aux_mat,im_aux_mat,non->psites,segments,sj,si,sj,N);
            /* Here trace should be */
                rate_response[t1]=trace_rate(re_aux_mat,non->singles);
		fprintf(ratefile,"%d %f\n",t1,rate_response[t1]);
//		fprintf(ratefile,"JAJE %f %f %f %f\n",re_aux_mat[0],re_aux_mat[1],re_aux_mat[2],re_aux_mat[3]);
//		fprintf(ratefile,"%f %f %f %f\n",J[0],J[1],J[2],J[3]);
//                fprintf(ratefile,"%f %f %f %f\n",re_Emi[0+nn2*t1],re_Emi[1+nn2*t1],re_Emi[2+nn2*t1],re_Emi[3+nn2*t1]);
//	        fprintf(ratefile,"%f %f %f %f\n",re_Abs[0+nn2*t1],re_Abs[1+nn2*t1],re_Abs[2+nn2*t1],re_Abs[3+nn2*t1]);	
            }
            /* Update rate matrix */
            rate=integrate_rate_response(rate_response,non->tmax)*non->deltat*icm2ifs*icm2ifs*twoPi*twoPi*1000;
            rate_matrix[si*segments+sj]=rate;
            rate_matrix[si*segments+si]=-rate;
        }

      }
    }
    fclose(ratefile);

    free(rate_response);
    free(ns);
    return;
}

/* Check if mcfret rates are in the incoherent limit */
void mcfret_validate(t_non *non);

/* Analyse rate matrix */
void mcfret_analyse(t_non *non);

/*This function will create a density matrix where every term is weighted with a Boltzmann weight*/
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int segments){
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
  float *Q,iQ;

  Q=(float *)calloc(segments,sizeof(float));  
  //printf("Seg %d \n",segments);
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
//      Q=Q+c2[a];
  }
  /* Find the inverse of the partition function */
//  iQ=1.0/Q;

  /*Transform back to site basis*/ 
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          cnr[b+a*N]+=H[b+a*N]*c2[b];
      }
  }  
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
          for (c=0;c<N;c++){
              density_matrix[a+c*N]+=H[b+a*N]*cnr[b+c*N];
          }
      }
  }
  
  /* Re-normalize */
  for (a=0;a<N;a++){
     Q[non->psites[a]]+=density_matrix[a+a*N];
  }
  for (a=0;a<N;a++){
      for (b=0;b<N;b++){
	  density_matrix[a+b*N]=density_matrix[a+b*N]/Q[non->psites[a]];
      }
  }      
    

  free(H);
  free(c2);
  free(e);
  free(cnr);
  free(Q);
  return;
}

/* Matrix multiplication for different segments */
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,
    float *rC,float *iC,int *psites,int segments,int si,int sk,int sj,int N){
    int i,j,k;
    /* Set initial values of results matrix to zero to be sure */
    clearvec(rC,N*N);
    clearvec(iC,N*N);
    for (i=0;i<N;i++){
        if (psites[i]==si){
            for (j=0;j<N;j++){
                if (psites[j]==sj){
                    for (k=0;k<N;k++){
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
  trace=0;
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
  fclose(file_handle);
}
