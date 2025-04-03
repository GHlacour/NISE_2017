#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif
//#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "CG_subs.h"
#include "calc_FD_CG_2DES.h"
#include <stdarg.h>
#include "project.h"
#include "propagate.h"
#include "read_trajectory.h"
#include "eq_den.h"
//#include "mpi.h"
#include "MPI_subs.h"
#include <complex.h>

/* ======================================================== */
/* MAIN ROUTINE INTEGRATING CG DOORWAY AND WINDOW FUNCTIONS */
/* ======================================================== */
void calc_FD_CG_2DES(t_non *non){
  /* This is the main routine which will control the flow of the FD
  coarse grained 2DES calculation. It uses the general NISE input. */
  int segments;
  int f_states;
  int e_states;
  int g_states;

  double *Q1;
  double *Q2;

  float *re_doorway, *im_doorway; 
  float *re_window_SE, * im_window_SE;
  float *re_window_GB, * im_window_GB;
  float *re_window_EA, * im_window_EA;
  int i;
  segments = project_dim(non);
  f_states = segments*(segments+1)/2;
  e_states = 2*segments;

  Q1 = (double *)calloc(segments,sizeof(double));
  Q2 = (double *)calloc(f_states,sizeof(double));

  printf("Reading in segments: %i \n", segments);
  printf("Reading in Quantum Yields from Q1s.dat Q2s.dat files.\n"); 
  read_in_QF_for_FD_CG_2DES(non,Q1,Q2,segments);

  re_doorway   = (float *)calloc(non->tmax*9*segments,sizeof(float));
  im_doorway   = (float *)calloc(non->tmax*9*segments,sizeof(float));
  re_window_SE = (float *)calloc(non->tmax*9*segments,sizeof(float));
  im_window_SE = (float *)calloc(non->tmax*9*segments,sizeof(float));   
  re_window_GB = (float *)calloc(non->tmax*9*segments,sizeof(float));
  im_window_GB = (float *)calloc(non->tmax*9*segments,sizeof(float)); 
  re_window_EA = (float *)calloc(non->tmax*9*segments,sizeof(float));
  im_window_EA = (float *)calloc(non->tmax*9*segments,sizeof(float)); 

  
  printf("Performing the FD_CG_2DES calculation.\n"); 

  if (!strcmp(non->technique, "FD_CG_2DES") ||  (!strcmp(non->technique, "FD_CG_2DES_doorway"))) {
      CG_doorway(non, re_doorway, im_doorway);
  }  
   if (!strcmp(non->technique, "FD_CG_2DES") ||  (!strcmp(non->technique, "FD_CG_2DES_window_SE")) ) {
      CG_window_SE(non, re_window_SE, im_window_SE); 
   } 
    if (!strcmp(non->technique, "FD_CG_2DES") ||  (!strcmp(non->technique, "FD_CG_2DES_window_GB"))){
      CG_window_GB(non, re_window_GB, im_window_GB); 
    }

    if (!strcmp(non->technique, "FD_CG_2DES") ||  (!strcmp(non->technique, "FD_CG_2DES_window_EA"))){
      CG_window_EA(non, re_window_EA, im_window_EA); 
    }

  /* Call the rate routine routine */
  if (!strcmp(non->technique, "FD_CG_2DES")||  (!strcmp(non->technique, "FD_CG_full_2DES_segments")) || (!strcmp(non->technique, "FD_CG_2DES_waitingtime"))){
    
      printf("Starting calculation of the 2DES with specific time delay\n");
      
      if ((!strcmp(non->technique, "FD_CG_2DES_waitingtime"))){
        /* Read in absorption, emission and coupling from file if needed */
        printf("Calculating spectroscopy from precalculated doorway function, window function\n");
        read_doorway_window_from_file(non,"CG_2DES_doorway.dat",im_doorway,re_doorway,non->tmax1);
        read_doorway_window_from_file(non,"CG_2DES_windows_EA.dat",im_window_EA,re_window_EA,non->tmax1);
        read_doorway_window_from_file(non,"CG_2DES_windows_GB.dat",im_window_GB,re_window_GB,non->tmax1);
        read_doorway_window_from_file(non,"CG_2DES_windows_SE.dat",im_window_SE,re_window_SE,non->tmax1);
        printf("Completed reading pre-calculated data.\n");
      }
      
      call_final_FD_CG_2DES(non,segments,re_doorway,im_doorway,re_window_SE,im_window_SE,re_window_GB,im_window_GB,
                            re_window_EA,im_window_EA, Q1, Q2);
  }
  
  free(re_doorway),      free(im_doorway);
  free(re_window_SE),    free(im_window_SE);
  free(re_window_GB),    free(im_window_GB);
  free(re_window_EA),    free(im_window_EA);  

  free(Q1),free(Q2);
  return;
}

void call_final_FD_CG_2DES(
  t_non *non,int segments,float *re_doorway,float *im_doorway,
  float *re_window_SE, float *im_window_SE,float *re_window_GB, float *im_window_GB,
  float *re_window_EA, float *im_window_EA, double *Q1, double *Q2){
    /* Define variables for multiple waiting time use */
    FILE *WTime;
    char waittime[16];
    int wfile;
    float *P_DA;

    P_DA = (float *)calloc(segments*segments,sizeof(float));
    /* Check if Waitingtime.dat is defined */
    sprintf(waittime,"");
    wfile=0;
    if (access("Waitingtime.dat", F_OK) !=-1){
      printf("Waitingtime file found by calc_FD_CG_2DES. \n");
      WTime=fopen("Waitingtime.dat","r");
      wfile=1;
    }
    /* While loop over waiting times */
    
    /* Change the value of non->tmax2 to the wanted waiting time */
    while (wfile>-1){
      /* Read new waiting time */
      if (wfile==1){
        if (fscanf(WTime,"%s",waittime)==1){ // Should this not be waittime instead of &waittime?
          printf("Calculating FD-CG2DES for %s fs \n",waittime);
          non->tmax2 = floor(atof(waittime)/(non->deltat));
        } else {
          fclose(WTime);
          break;
        }
      } else {
        wfile=-1;
      }
      
      CG_P_DA(non,P_DA,segments);
      
      FD_CG_full_2DES_segments(non,re_doorway,im_doorway,re_window_SE,im_window_SE,
        re_window_GB,im_window_GB,re_window_EA,im_window_EA,P_DA,segments,waittime,wfile,Q1,Q2);
    }
    free(P_DA);
}

/* -=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-. */
/* READ IN QUANTUM YIELDS FROM SEPARATE Q1 AND Q2 FILES */
/* -=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-._.-=-. */
void read_in_QF_for_FD_CG_2DES(t_non *non, double *Q1, double *Q2, int segments){
    /* This function is used to read in the quantum yields Q1s
        and Q2s and make them available for further calculations*/
    /* Define variables for QFactors */
    int f_states;
    int e_states;
    int g_states;

    f_states = segments*(segments+1)/2;
    e_states = 2*segments;
    g_states = 3;

    /* Open .dat files */
    FILE *Q1sFile = fopen("Q1s.dat","r");
    if (!Q1sFile) {
        fprintf(stderr, "File not found: Q1s.dat\n");
        return;
        exit(0);
    }
    FILE *Q2sFile = fopen("Q2s.dat","r");
    if (!Q2sFile) {
        fprintf(stderr, "File not found: Q2s.dat\n");
        return;
    }

    /* Read-in and store values for Q1s and Q2s */
    for (int i = 0; i < segments; ++i) {
        fscanf(Q1sFile, "%lf", &Q1[i]);
    }
    fclose(Q1sFile);

    for (int i = 0; i < f_states; ++i) {
        fscanf(Q2sFile, "%lf", &Q2[i]);
    }
    fclose(Q2sFile);
    
    /* Display the values for demonstration purposes.
    Can modify printLevel in input file : PrintLevel 2 */
    if(non->printLevel>1){
        printf("Q1 values:\n");
        for (int i = 0; i < segments; ++i) {
            printf("Q1_%d = %f\n", i + 1, Q1[i]);
        }
            
        printf("Q2 values:\n");
        for (int i = 0; i < f_states; ++i) {
            printf("Q2_%d = %f\n", i + 1, Q2[i]);
        }
    }
}



/* ==================================================== */
/* COMBINE DOORWAY AND WINDOW FUNCTIONS FOR CG SEGMENTS */
/* ==================================================== */
void FD_CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
    float *re_window_SE,float *im_window_SE,float *re_window_GB, float *im_window_GB,
    float *re_window_EA,float *im_window_EA,float *P_DA,int N, char *waittime,int wfile,
    double *Q1, double *Q2){
    int t1,t2,t3;
    int S,R,W; // Segment indices
    int indext1,indext2,indext3;
    int sampleCount=1;
    int pol,molPol;
    int px[4];
    int RR,RW;
    float **re_2DES_NR_sum,  **re_2DES_R_sum, **im_2DES_NR_sum, **im_2DES_R_sum;
    float polWeight;
    float factorWEA1;
    float factorWEA2;
    float factorR;
    float factorRR;

    re_2DES_NR_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    re_2DES_R_sum  = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_NR_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_R_sum  = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    
    /* We repeat everyting for each waiting time t2 */
    t2 = non->tmax2-1;
    /* We repeat everyting for each labframe polarization */
    for (pol=0;pol<3;pol++){
      /* Loop over the 21 microscopic polarizations */
	    for (molPol=0;molPol<21;molPol++){
        polWeight=polarweight(pol,molPol);
	      polar(px,molPol);	  
	      /* Loop over times */
	      for (t1=0;t1<non->tmax1;t1++){
          for (t3=0;t3<non->tmax3;t3++){
            /* Loop over segments */
            for (S=0; S<N; S++){
              for (R=0; R<N; R++){
                /* First do GB */
		            indext1=CG_index(non,S,px[0],px[1],t1);
		            indext2=R*N+S;
                factorR=Q1[R]*polWeight*P_DA[indext2];
		            indext3=CG_index(non,R,px[2],px[3],t3);
                RR = Sindex(R, R, N); // R + R * ((N * 2) - R - 1) / 2; // Sindex(R, R, N);
                
		            /* Ground state bleach */
                re_2DES_R_sum[t3][t1]-=factorR*re_doorway[indext1]*re_window_GB[indext3];
                re_2DES_R_sum[t3][t1]+=factorR*im_doorway[indext1]*im_window_GB[indext3];
                im_2DES_R_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_GB[indext3];
                im_2DES_R_sum[t3][t1]+=factorR*im_doorway[indext1]*re_window_GB[indext3];
                re_2DES_NR_sum[t3][t1]-=factorR*re_doorway[indext1]*re_window_GB[indext3];
                re_2DES_NR_sum[t3][t1]-=factorR*im_doorway[indext1]*im_window_GB[indext3];
                im_2DES_NR_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_GB[indext3];
                im_2DES_NR_sum[t3][t1]-=factorR*im_doorway[indext1]*re_window_GB[indext3];
                
		            /* Stimulated Emission */
                re_2DES_R_sum[t3][t1]-=factorR*re_doorway[indext1]*re_window_SE[indext3];
                re_2DES_R_sum[t3][t1]+=factorR*im_doorway[indext1]*im_window_SE[indext3];
                im_2DES_R_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_SE[indext3];
                im_2DES_R_sum[t3][t1]+=factorR*im_doorway[indext1]*re_window_SE[indext3];
                re_2DES_NR_sum[t3][t1]-=factorR*re_doorway[indext1]*re_window_SE[indext3];
                re_2DES_NR_sum[t3][t1]-=factorR*im_doorway[indext1]*im_window_SE[indext3];
                im_2DES_NR_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_SE[indext3];
                im_2DES_NR_sum[t3][t1]-=factorR*im_doorway[indext1]*re_window_SE[indext3];
                
                /* Excited State Absorption 1 */
                factorR=-factorR;
                re_2DES_R_sum[t3][t1]+=factorR*re_doorway[indext1]*re_window_EA[indext3];
                re_2DES_R_sum[t3][t1]+=factorR*im_doorway[indext1]*im_window_EA[indext3];
                im_2DES_R_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_EA[indext3];
                im_2DES_R_sum[t3][t1]-=factorR*im_doorway[indext1]*re_window_EA[indext3];
                re_2DES_NR_sum[t3][t1]+=factorR*re_doorway[indext1]*re_window_EA[indext3];
                re_2DES_NR_sum[t3][t1]-=factorR*im_doorway[indext1]*im_window_EA[indext3];
                im_2DES_NR_sum[t3][t1]+=factorR*re_doorway[indext1]*im_window_EA[indext3];
                im_2DES_NR_sum[t3][t1]+=factorR*im_doorway[indext1]*re_window_EA[indext3];

                /* Excited State Absorption 2 */
                /* like Ground state bleach */
                  factorRR = Q2[RR]*polWeight*P_DA[indext2];
                  re_2DES_R_sum[t3][t1]+=factorRR*re_doorway[indext1]*re_window_EA[indext3];
                  re_2DES_R_sum[t3][t1]+=factorRR*im_doorway[indext1]*im_window_EA[indext3];
                  im_2DES_R_sum[t3][t1]+=factorRR*re_doorway[indext1]*im_window_EA[indext3];
                  im_2DES_R_sum[t3][t1]-=factorRR*im_doorway[indext1]*re_window_EA[indext3];
                  re_2DES_NR_sum[t3][t1]+=factorRR*re_doorway[indext1]*re_window_EA[indext3];
                  re_2DES_NR_sum[t3][t1]-=factorRR*im_doorway[indext1]*im_window_EA[indext3];
                  im_2DES_NR_sum[t3][t1]+=factorRR*re_doorway[indext1]*im_window_EA[indext3];
                  im_2DES_NR_sum[t3][t1]+=factorRR*im_doorway[indext1]*re_window_EA[indext3];
                for (W=0; W<N; W++){
                  if (W != R) {                    
                    RW = Sindex(R,W,N); //R + W * ((N * 2) - W - 1) / 2; //Sindex(R,W,21);
                    factorWEA1=2*Q1[R]*polWeight*P_DA[indext2];
                    indext3=CG_index(non,W,px[2],px[3],t3);
                    /* Excited State Absorption 1 */
                    re_2DES_R_sum[t3][t1]-=factorWEA1*re_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_R_sum[t3][t1]+=factorWEA1*im_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_R_sum[t3][t1]+=factorWEA1*re_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_R_sum[t3][t1]+=factorWEA1*im_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_NR_sum[t3][t1]-=factorWEA1*re_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_NR_sum[t3][t1]-=factorWEA1*im_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_NR_sum[t3][t1]+=factorWEA1*re_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_NR_sum[t3][t1]-=factorWEA1*im_doorway[indext1]*re_window_GB[indext3];
                    /* Excited State Absorption 2 */
                    /* like Ground state bleach */
                    factorWEA2=(-1)*Q2[RW]*polWeight*P_DA[indext2];
                    re_2DES_R_sum[t3][t1]-=factorWEA2*re_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_R_sum[t3][t1]+=factorWEA2*im_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_R_sum[t3][t1]+=factorWEA2*re_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_R_sum[t3][t1]+=factorWEA2*im_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_NR_sum[t3][t1]-=factorWEA2*re_doorway[indext1]*re_window_GB[indext3];
                    re_2DES_NR_sum[t3][t1]-=factorWEA2*im_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_NR_sum[t3][t1]+=factorWEA2*re_doorway[indext1]*im_window_GB[indext3];
                    im_2DES_NR_sum[t3][t1]-=factorWEA2*im_doorway[indext1]*re_window_GB[indext3];
                  }
                }
              }
            }
          }
        }
      }
    
    char WFileName[256];
	  /* Write response functions to file */
    if (wfile==1){
      if (pol==0){
        sprintf(WFileName,"RparI_%sfs.dat",waittime);
        print2D(WFileName, im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        sprintf(WFileName,"RparII_%sfs.dat",waittime);
        print2D(WFileName,  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      } else if (pol==1){
        sprintf(WFileName,"RperI_%sfs.dat",waittime);
        print2D(WFileName, im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        sprintf(WFileName,"RperII_%sfs.dat",waittime);
        print2D(WFileName,  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      } else{
        sprintf(WFileName,"RcroI_%sfs.dat",waittime);
        print2D(WFileName, im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        sprintf(WFileName,"RcroII_%sfs.dat",waittime);
        print2D(WFileName,  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      }
    } else{
      if (pol==0){
        print2D("RparI.dat", im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        print2D("RparII.dat",  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      } else if (pol==1){
        print2D("RperI.dat", im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        print2D("RperII.dat",  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      } else{
        print2D("RcroI.dat", im_2DES_R_sum,  re_2DES_R_sum,  non, sampleCount);
        print2D("RcroII.dat",  im_2DES_NR_sum, re_2DES_NR_sum, non, sampleCount);
      }
    }
    

	  /* Loop over times to clean response functions */
    for (t1=0;t1<non->tmax1;t1++){
      for (t3=0;t3<non->tmax3;t3++){
        re_2DES_R_sum[t3][t1]=0;
        im_2DES_R_sum[t3][t1]=0;
        re_2DES_NR_sum[t3][t1]=0;
        im_2DES_NR_sum[t3][t1]=0;
      }
    }
    }
    
    free(re_2DES_NR_sum) ,free(re_2DES_R_sum),free(im_2DES_NR_sum) ,free(im_2DES_R_sum);
  }