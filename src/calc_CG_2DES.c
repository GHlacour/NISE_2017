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
#include "calc_CG_2DES.h"
#include <stdarg.h>
#include "project.h"
#include "propagate.h"
#include "read_trajectory.h"
#include "eq_den.h"
#include "CG_subs.h"
//#include "mpi.h"
#include "MPI_subs.h"
#include <complex.h>

    
/* This is the main routine which will control the flow of the 
coarse grained 2DES calculation. It uses the general NISE input. */
void calc_CG_2DES(t_non *non){
    float *re_doorway, *im_doorway; 
    float *re_window_SE, * im_window_SE;
    float *re_window_GB, * im_window_GB;
    float *re_window_EA, * im_window_EA;
    float *re_2DES , *im_2DES;
    float *P_DA;
    int pro_dim,i;
    pro_dim=project_dim(non);
    
    re_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    re_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));   
    re_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float)); 
    re_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float)); 
    P_DA=(float *)calloc(pro_dim*pro_dim,sizeof(float));
    
    printf("Performing the CG_2DES calculation.\n"); 

    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway"))) {
        CG_doorway(non, re_doorway, im_doorway);
      }  
     if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ) {
        CG_window_SE(non, re_window_SE, im_window_SE); 
     } 
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_GB"))){
        CG_window_GB(non, re_window_GB, im_window_GB); 
      }
  
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_EA"))){
        CG_window_EA(non, re_window_EA, im_window_EA); 
      }
    /* Call the rate routine routine */
    if (!strcmp(non->technique, "CG_2DES")||  (!strcmp(non->technique, "CG_full_2DES_segments")) || (!strcmp(non->technique, "CG_2DES_waitingtime"))){
      CG_P_DA(non,P_DA,pro_dim);
        printf("Starting calculation of the 2DES with spcifical time delay\n");
        if ((!strcmp(non->technique, "CG_2DES_waitingtime"))){
            /* Read in absorption, emission and coupling from file if needed */
	    printf("Calculating spectroscopy from precalculated doorway function, window function\n");
	    read_doorway_window_from_file(non,"CG_2DES_doorway.dat",im_doorway,re_doorway,non->tmax1);
      read_doorway_window_from_file(non,"CG_2DES_windows_EA.dat",im_window_EA,re_window_EA,non->tmax1);
      read_doorway_window_from_file(non,"CG_2DES_windows_GB.dat",im_window_GB,re_window_GB,non->tmax1);
      read_doorway_window_from_file(non,"CG_2DES_windows_SE.dat",im_window_SE,re_window_SE,non->tmax1);
      printf("Completed reading pre-calculated data.\n");
        }
        call_final_CG_2DES(non,P_DA,pro_dim,re_doorway,im_doorway,
                              re_window_SE,im_window_SE,
                              re_window_GB,im_window_GB,
                              re_window_EA,im_window_EA,
                              re_2DES, im_2DES);
        /*CG_full_2DES_segments(non,re_doorway,im_doorway,
                              re_window_SE,im_window_SE,
                              re_window_GB, im_window_GB,
                              re_window_EA,im_window_EA,
				                      P_DA,pro_dim);*/
    }
    free(re_doorway),      free(im_doorway);
    free(re_window_SE),    free(im_window_SE);
    free(re_window_GB),    free(im_window_GB);
    free(re_window_EA),    free(im_window_EA);  
    free(P_DA);
    return;
}

/* This routine control the calculation for different waiting times. 
   It is called by the main routine after the reservation of memory for
   various variables. */
void call_final_CG_2DES(
  t_non *non,float *P_DA,int pro_dim,float *re_doorway,float *im_doorway,
  float *re_window_SE, float *im_window_SE,float *re_window_GB, float *im_window_GB,
  float *re_window_EA, float *im_window_EA,float *re_2DES , float *im_2DES){

  /* Define variables for multiple waiting time use */
  FILE *WTime;
  char waittime[16];
  int wfile;
                          
  /* Check if Waitingtime.dat is defined */
  sprintf(waittime,"");
  wfile=0;
  if (access("Waitingtime.dat", F_OK) !=-1){
    printf("Waitingtime file found by calc_CG_2DES. \n");
    WTime=fopen("Waitingtime.dat","r");
    wfile=1;
  }
  /* While loop over waiting times */
  /* Change the value of non->tmax2 to the wanted waiting time */
  while (wfile>-1){
    /* Read new waiting time */
    if (wfile==1){
      if (fscanf(WTime,"%s",waittime)==1){ // Should this not be waittime instead of &waittime?
        printf("Calculating CG2DES for %s fs \n",waittime);
        non->tmax2 = floor(atof(waittime)/(non->deltat));
      } else {
        fclose(WTime);
        break;
      }
    } else {
      wfile=-1;
    }

    CG_P_DA(non,P_DA,pro_dim);

    CG_full_2DES_segments(non,re_doorway,im_doorway,re_window_SE,im_window_SE,
      re_window_GB,im_window_GB,re_window_EA,im_window_EA,P_DA,pro_dim,waittime,wfile);
  }
}


/* Combine the doorway and window functions for the segments */
void CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
  float *re_window_SE,float *im_window_SE,float *re_window_GB,float *im_window_GB,
  float *re_window_EA,float *im_window_EA,float *P_DA,int N,char *waittime, int wfile){

  int t1,t2,t3;
  int S,R; // Segment indices
  int indext1,indext2,indext3;
  int sampleCount=1;
  int pol,molPol;
  int px[4];
  float **re_2DES_NR_sum,  **re_2DES_R_sum, **im_2DES_NR_sum, **im_2DES_R_sum;
  float polWeight;

  re_2DES_NR_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  re_2DES_R_sum  = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_NR_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_R_sum  = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));

  /* We repeat everyting for each waiting time t2 */
  //for (t2=0;t2<non->tmax2;t2++){
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
	        for (S=0;S<N;S++){
		        for (R=0;R<N;R++){
	            /* First do GB */
              //indext1=S*9*non->tmax1+px[0]*3*non->tmax1+px[1]*non->tmax1+t1;
		          indext1=CG_index(non,S,px[0],px[1],t1);
		          indext2=R*N+S;
	            //indext3=R*9*non->tmax3+px[2]*3*non->tmax3+px[3]*non->tmax3+t3;
		          indext3=CG_index(non,R,px[2],px[3],t3);
		          /* Ground state bleach */
              re_2DES_R_sum[t3][t1]-=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_GB[indext3];
              re_2DES_R_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_GB[indext3];
              im_2DES_R_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_GB[indext3];
              im_2DES_R_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_GB[indext3];
              re_2DES_NR_sum[t3][t1]-=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_GB[indext3];
              re_2DES_NR_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_GB[indext3];
              im_2DES_NR_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_GB[indext3];
              im_2DES_NR_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_GB[indext3];

		          /* Stimulated Emission */
              re_2DES_R_sum[t3][t1]-=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_SE[indext3];
              re_2DES_R_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_SE[indext3];
              im_2DES_R_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_SE[indext3];
              im_2DES_R_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_SE[indext3];
              re_2DES_NR_sum[t3][t1]-=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_SE[indext3];
              re_2DES_NR_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_SE[indext3];
              im_2DES_NR_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_SE[indext3];
              im_2DES_NR_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_SE[indext3];
		          
              /* Excited State Absorption */
              re_2DES_R_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_EA[indext3];
              re_2DES_R_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_EA[indext3];
              im_2DES_R_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_EA[indext3];
              im_2DES_R_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_EA[indext3];
              re_2DES_NR_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*re_window_EA[indext3];
              re_2DES_NR_sum[t3][t1]-=polWeight*im_doorway[indext1]*P_DA[indext2]*im_window_EA[indext3];
              im_2DES_NR_sum[t3][t1]+=polWeight*re_doorway[indext1]*P_DA[indext2]*im_window_EA[indext3];
              im_2DES_NR_sum[t3][t1]+=polWeight*im_doorway[indext1]*P_DA[indext2]*re_window_EA[indext3];
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

