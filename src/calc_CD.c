#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "calc_CD.h"
#include "1DFFT.h"
#include "project.h"

void calc_CD(t_non *non){
  // Initialize variables
  float *re_S_1,*im_S_1; // The first-order response function
  float *re_S_1j,*im_S_1j; // Before avaraging
  float *mu_eg,*Hamil_i_e;
  float *mu_p;
  float *pos;
  float posj;
  // Aid arrays
  float *vecr,*veci;
  //,*vecr_old,*veci_old;

  /* Floats */
  float shift1;
  /* 1D Fourier transform */
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;

  /* File handles */
  FILE *H_traj,*mu_traj,*pos_traj;
  FILE *outone,*log;
  FILE *Cfile;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int x,ti,tj,i;
  int y,z;
  int j,N;
  int t1,fft;
  int elements;
  int cl,Ncl;
  int sign;
  int pro_dim,ip;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

  /* Check for projection */
  pro_dim=project_dim(non);

  // Allocate memory
  re_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));
  im_S_1=(float *)calloc(non->tmax*pro_dim,sizeof(float));  
  nn2=non->singles*(non->singles+1)/2;
  N=non->singles;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));
  re_S_1j=(float *)calloc(non->tmax*N*pro_dim,sizeof(float));
  im_S_1j=(float *)calloc(non->tmax*N*pro_dim,sizeof(float));

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

  pos_traj=fopen(non->positionFName,"rb");
  if (pos_traj==NULL){
    printf("Position file %s not found!\n",non->positionFName);
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
  if (non->end>N_samples){
    printf(RED "Endpoint larger than number of samples was specified.\n" RESET);
    printf(RED "Endpoint was %d but cannot be larger than %d.\n" RESET,non->end,N_samples);
    exit(0);
  }

  log=fopen("NISE.log","a");
  fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
  fclose(log);

  vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles*non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_p=(float *)calloc(non->singles,sizeof(float));
  pos=(float *)calloc(non->singles,sizeof(float));

  printf("\n Note that the CD implementation assumes that the positions of\n");
  printf("the full system specified in the Position file is contained\n");
  printf("in a box as periodic boundary contitions are NOT applied.\n\n");

  // Loop over samples
  for (samples=non->begin;samples<non->end;samples++){

    // Calculate linear response    
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

      // Loop over polarizations of the initial excitation      
      for (x=0;x<3;x++){
        // Read mu(ti)
        if (read_mue(non,mu_eg,mu_traj,ti,x)!=1){
  	  printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	  printf("ITIME %d %d\n",ti,x);
	  exit(1);
        }
        // Initialize excitation on initial site
        clearvec(vecr,non->singles*non->singles);
        clearvec(veci,non->singles*non->singles);
        
        // Loop over initial sites
        for (j=0;j<N;j++){
          vecr[j+j*N]=mu_eg[j];
        }
        
        // Loop over delay
        for (t1=0;t1<non->tmax;t1++){
	  tj=ti+t1;
	  // Read Hamiltonian
	  if (read_He(non,Hamil_i_e,H_traj,tj)!=1){
	    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
            exit(1);
	  }

          // Loop over polarization values y for mu
          for (y=0;y<3;y++){
            // Exclude values taken up by the first interaction
            if (y!=x){	
              // Find corresponding value for the polarization used for the distance matrix
              z=3-x-y;
	      // Read mu(tj)
	      if (read_mue(non,mu_eg,mu_traj,tj,y)!=1){
	        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
	        printf("JTIME %d %d\n",tj,y);
	        exit(1);
              }
              // Read positions
              if (read_mue(non,pos,pos_traj,tj,z)!=1){
                printf("Position trajectory file to short, could not fill buffer!!!\n");
                printf("JTIME %d %d\n",tj,z);
                exit(1);
              }
 //             posj=pos[j];
              // Do projection on selected sites if asked
/*	      if (non->Npsites>0){
	        projection(mu_eg,non);
              }*/
              sign=0;
              // Determine the sign
              if (z==0 & y==1 & x==2){sign=1;}//{sign=1;}
              if (z==0 & y==2 & x==1){sign=-1;}//{}sign=-1;}
              if (z==1 & y==0 & x==2){sign=-1;}//{sign=-1;}
              if (z==2 & y==0 & x==1){sign=1;}//{sign=1;}
              if (z==2 & y==1 & x==0){sign=-1;}//{sign=-1;}
              if (z==1 & y==2 & x==0){sign=1;}//{sign=1;}
              if (sign==0){
                 printf(RED "Bug in CD routine.\n" RESET);
                 exit(1);
              }

              if (non->Npsites==0){
              /* Find response without projection */
                calc_CD1(re_S_1j,im_S_1j,t1,non,vecr,veci,mu_eg,pos,sign);
              } else if (non->Npsites<non->singles){
                projection(mu_eg,non);
              /* Find response with projection on single segment */
                calc_CD1(re_S_1j,im_S_1j,t1,non,vecr,veci,mu_eg,pos,sign);
              } else {
              /* Find response with projection on multiple segments */
                for (ip=0;ip<pro_dim;ip++){
                   multi_projection(mu_eg,mu_p,non,ip);
                   calc_CD1(re_S_1j+non->tmax*N*ip,im_S_1j+non->tmax*N*ip,t1,non,vecr,veci,mu_eg,pos,sign);
               }
              }
	      
              sign=0;
              // Determine the sign
              if (z==0 & y==1 & x==2){sign=1;}//{sign=1;}
              if (z==0 & y==2 & x==1){sign=-1;}//{}sign=-1;}
              if (z==1 & y==0 & x==2){sign=-1;}//{sign=-1;}
              if (z==2 & y==0 & x==1){sign=1;}//{sign=1;}
              if (z==2 & y==1 & x==0){sign=-1;}//{sign=-1;}
              if (z==1 & y==2 & x==0){sign=1;}//{sign=1;}
              if (sign==0){
                printf("Bug in CD routine.\n");
                exit(1);
              }

	      // Find response
//              calc_CD1(re_S_1j,im_S_1j,t1,non,vecr,veci,mu_eg,pos,sign);
	    }
          }
          
	  // Propagate vector
#pragma omp parallel for
          for (j=0;j<non->singles;j++){
	    if (non->propagation==1) propagate_vec_coupling_S(non,Hamil_i_e,vecr+j*N,veci+j*N,non->ts,1);
	    if (non->propagation==0){
	      if (non->thres==0 || non->thres>1){
	        propagate_vec_DIA(non,Hamil_i_e,vecr+j*N,veci+j*N,1);
	      } else {
	        elements=propagate_vec_DIA_S(non,Hamil_i_e,vecr+j*N,veci+j*N,1);
	        if (samples==non->begin){
	          if (t1==0){
		    if (x==0){
                      if (j==0){
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
        }
      }
    } // Cluster loop
  
    // Log time
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);
          
    time_now=log_time(time_now,log);
    fclose(log);
  }

  free(vecr);
  free(veci);
  free(mu_eg);
  free(Hamil_i_e);

  // Sum over results from different chromophores
  for (t1=0;t1<non->tmax;t1++){
    for (j=0;j<N;j++){
      for (ip=0;ip<pro_dim;ip++){
         re_S_1[t1+ip*non->tmax]+=re_S_1j[t1+j*non->tmax+ip*N*non->tmax];
         im_S_1[t1+ip*non->tmax]+=im_S_1j[t1+j*non->tmax+ip*N*non->tmax]; 
      }
    }
  } 

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

  fclose(mu_traj),fclose(H_traj),fclose(pos_traj);
  if (non->cluster!=-1){
    fclose(Cfile);
  }

  outone=fopen("TD_CD.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
/*    fprintf(outone,"%f %e %e\n",t1*non->deltat,re_S_1[t1]/samples,im_S_1[t1]/samples); */
      fprintf(outone,"%f ",t1*non->deltat);
      for (ip=0;ip<pro_dim;ip++){
         fprintf(outone,"%e %e ",re_S_1[t1+ip*non->tmax]/samples,im_S_1[t1+ip*non->tmax]/samples);
      }
      fprintf(outone,"\n");
  }
  fclose(outone);

  /* Do Forier transform and save */
  do_1DFFT(non,"CD.dat",re_S_1,im_S_1,samples);

  free(re_S_1),free(im_S_1);
  free(re_S_1j),free(im_S_1j);

  printf("----------------------------------------------\n");
  printf(" CD calculation succesfully completed\n");
  printf("----------------------------------------------\n\n");

  return;
}	

/* Do the final CD calculation */
/* Currently neglecting periodic boundary conditions */
void calc_CD1(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *mu,float *pos,int sign){
  int i,j;  
#pragma omp parallel for
  for (j=0;j<non->singles;j++){
    for (i=0;i<non->singles;i++){
      re_S_1[t1+j*non->tmax]+=sign*(pos[i]-pos[j])*mu[i]*cr[i+j*non->singles];
      im_S_1[t1+j*non->tmax]+=sign*(pos[i]-pos[j])*mu[i]*ci[i+j*non->singles];
    }
  }
  return;
}
