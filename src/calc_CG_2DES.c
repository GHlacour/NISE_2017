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
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "calc_CG_2DES.h"
#include <stdarg.h>
#include "project.h"
#include "propagate.h"
#include "eq_den.h"

/* This is the main routine which will control the flow of the
 * coarse grained 2DES calculation. It takes the general NISE
 * input as input.
 */

void calc_CG_2DES(t_non *non){
    float *re_doorway, *im_doorway; 
    float *re_window_SE, * im_window_SE;
    int pro_dim;
    re_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    re_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));    
    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
     (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      {
        CG_2DES_doorway(non, re_doorway, im_doorway);
        //CG_2DES_window_SE(non, re_window_SE, im_window_SE); 
      }  
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
     (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      {
        CG_2DES_window_SE(non, re_window_SE, im_window_SE); 
      }    
    free(re_doorway);
    free(im_doorway);
    free(re_window_SE);
    free(im_window_SE);
    return;
}

void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway){  /* what *non did?*/
  /* Initialize variables*/
  //float *re_doorway,*im_doorway; /* The dorrway part*/
  float *Hamil_i_e;
  float *mu_eg;
  float *vecr,*veci;
  float *mu_xyz;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;
  /* Floats */
  float shift1;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int alpha,beta,a_b,ti,tj,i,j; /*alpha and beta are corresponding the dipole for different segment,  a_b  = alpha*beta, which used in the writing file*/
  int t1,t2;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;
  int a,index,seg,site_num,seg_num;/*site num is the number of the site, which used to make sure the index number later*/

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;


  /* check projection */   
  if (non->Npsites!=non->singles){
      printf("Input segment number and the projection segment number are different.\n");
      printf("please check the input file or the projection file\n");
      exit(1);
  }
  /* projection */
  pro_dim=project_dim(non);
  /* Allocate memory*/  /*(tmax+1)*9*/
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));

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
   /*printf("A\n");*/
  /* Open file with cluster information if appicable */
  if (non->cluster!=-1){
    Cfile=fopen("Cluster.bin","rb");
    if (Cfile==NULL){
      printf("Cluster option was activated but no Cluster.bin file provided.\n");
      printf("Please, provide cluster file or remove Cluster keyword from\n");
      printf("input file.\n");
      exit(0);
    }
    Ncl=0; /* Counter for snapshots calculated*/
  } 

  vecr=(float *)calloc(non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_xyz=(float *)calloc(non->singles*3,sizeof(float));

  /* Here we walsnt to call the routine for checking the trajectory files*/
  control(non);
  itime =0;

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



  /* Read coupling, this is done if the coupling and transition-dipoles are *
   * time-independent and only one snapshot is stored */
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

    /* Reading in single fixed transition dipole vector matrix */
    for (alpha=0;alpha<3;alpha++){
      if (read_mue(non,mu_xyz+non->singles*alpha,mu_traj,0,alpha)!=1){
         printf("Dipole trajectory file to short, could not fill buffer!!!\n");
         printf("ITIME %d %d\n",0,alpha);
         exit(1);
      }
    }
  } 

    /* Loop over samples */
 for (samples=non->begin;samples<non->end;samples++){
      /* Calculate linear response */   
    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
        printf("Cluster trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d\n",ti);
        exit(1);
      }
     
        // Configuration belong to cluster
      if (non->cluster==cl){
        Ncl++;
      }
    }

      // Include snapshot if it is in the cluster or if no clusters are defined
    if (non->cluster==-1 || non->cluster==cl){   
      for (alpha=0;alpha<3;alpha++){
         /* Read mu(ti) */
        if (!strcmp(non->hamiltonian,"Coupling")){
          copyvec(mu_xyz+non->singles*alpha,vecr,non->singles);
        } else {
          if (read_mue(non,vecr,mu_traj,ti,alpha)!=1){
             printf("Dipole trajectory file to short, could not fill buffer!!!\n");
             printf("ITIME %d %d\n",ti,alpha);
             exit(1);
            }
          }       
          /* this is for time t1 to generate the vector for the dipole with 0 as the imagine part*/ 
          clearvec(veci,non->singles);
          /*this is just copy the input dipole to the real part of the vector*/
          /*This two step is necessary as the input dipole is real number, but after probagate it becomes complex number */
          //copyvec(vecr,mu_eg,non->singles);
          /* Loop over coherence time */

        
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
      
            for (beta=0;beta<3;beta++){
              /* Read mu(tj) */
              if (!strcmp(non->hamiltonian,"Coupling")){
                copyvec(mu_xyz+non->singles*beta,mu_eg,non->singles);
              } else {
                if (read_mue(non,mu_eg,mu_traj,tj,beta)!=1){
                  printf("Dipole trajectory file to short, could not fill buffer!!!\n");
                  printf("JTIME %d %d\n",tj,beta);
                  exit(1);
                }
              }
             
            /* Here calculate doorway function/
            /* Inner product for all sites*/
            /*here we need to make clear the seg_num and the site_num
            the number of the position should be decided by the segment number,
            we should sum all the value in one segment*/

             for (site_num=0;site_num<non->singles;site_num++){
                seg_num=non->psites[site_num];
                /*this equation is make sure the calculated data in the right position */
                index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                re_doorway[index]+=mu_eg[site_num]*vecr[site_num];
                im_doorway[index]+=mu_eg[site_num]*veci[site_num]; 

              }
          }  
          
          /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
          if (non->Npsites==non->singles){
            zero_coupling(Hamil_i_e,non);
          } else {
            printf("Segment number and the projection number are different");
            exit(1);
          }
          /* Propagate dipole moment */
          propagate_vector(non,Hamil_i_e,vecr,veci,1,samples,t1*alpha);
          
      }
    }
  }

    /* Update Log file with time and sample numner */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);        
    time_now=log_time(time_now,log);
    fclose(log);

  }
 
  /* The calculation is finished we can close all auxillary arrays before writing */
  /* output to file. */
  free(vecr);
  free(veci);
  free(mu_eg);
  free(mu_xyz);
  free(Hamil_i_e);

  /* The calculation is finished, lets write output */
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response 123!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  /* Print information on number of realizations included belonging to the selected */
  /* cluster and close the cluster file. (Only to be done if cluster option is active.) */
  if (non->cluster!=-1){
    printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
    fclose(Cfile);
  }

  /* Close Trajectory Files */
  fclose(mu_traj),fclose(H_traj);

  /* Save  the imaginary part for time domain response */
  outone=fopen("CG_2DES_doorway_im.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index =  seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",im_doorway[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }

  /* Save the real part for time domain response */
  outone=fopen("CG_2DES_doorway_re.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index =  seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",re_doorway[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }
 fclose(outone);
}

void CG_2DES_window_SE(t_non *non, float * re_window_SE, float * im_window_SE){
   /* Initialize variables*/
 /* The window part for SE*/
  float *Hamil_i_e;
  float *mu_eg;
  float *vecr,*veci;
  float *mu_xyz;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;
  /* Floats */
  float shift1;

  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int alpha,beta,a_b,ti,tj,i,j; /*alpha and beta are corresponding the dipole for different segment,  a_b  = alpha*beta, which used in the writing file*/
  int t1,t2;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;
  int a,b,index,seg,site_num,seg_num;/*site num is the number of the site, which used to make sure the index number later*/
// reserve memory for the density operator
  int N;
  N = non->singles;
  float *rho_l;
  rho_l=(float *)calloc(N*N,sizeof(float));

  //here the mid_vcr is just used  as the mid part for multiply the dipole with density operator
  float *mid_ver;
  mid_ver = (float *)calloc(N,sizeof(float));

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;
// reserve memory for transition dipole
  vecr=(float *)calloc(non->singles,sizeof(float));	
  veci=(float *)calloc(non->singles,sizeof(float));
  mu_eg=(float *)calloc(non->singles,sizeof(float));
  mu_xyz=(float *)calloc(non->singles*3,sizeof(float));

  /* check projection */   
  if (non->Npsites!=non->singles){
      printf("Input segment number and the projection segment number are different.\n");
      printf("please check the input file or the projection file\n");
      exit(1);
  }
  /* projection */
  pro_dim=project_dim(non);
  /* Allocate memory*/  /*(tmax+1)*9*/
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));

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
   /*printf("A\n");*/
  /* Open file with cluster information if appicable */
  if (non->cluster!=-1){
    Cfile=fopen("Cluster.bin","rb");
    if (Cfile==NULL){
      printf("Cluster option was activated but no Cluster.bin file provided.\n");
      printf("Please, provide cluster file or remove Cluster keyword from\n");
      printf("input file.\n");
      exit(0);
    }
    Ncl=0; /* Counter for snapshots calculated*/
  } 



  /* Here we walsnt to call the routine for checking the trajectory files*/
  control(non);
  itime =0;

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



  /* Read coupling, this is done if the coupling and transition-dipoles are *
   * time-independent and only one snapshot is stored */
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

    /* Reading in single fixed transition dipole vector matrix */
    for (alpha=0;alpha<3;alpha++){
      if (read_mue(non,mu_xyz+non->singles*alpha,mu_traj,0,alpha)!=1){
         printf("Dipole trajectory file to short, could not fill buffer!!!\n");
         printf("ITIME %d %d\n",0,alpha);
         exit(1);
      }
    }
  } 

    /* Loop over samples */
 for (samples=non->begin;samples<non->end;samples++){
      /* Calculate linear response */   
    ti=samples*non->sample;
    if (non->cluster!=-1){
      if (read_cluster(non,ti,&cl,Cfile)!=1){
        printf("Cluster trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d\n",ti);
        exit(1);
      }
     
        // Configuration belong to cluster
      if (non->cluster==cl){
        Ncl++;
      }
    }

      // Include snapshot if it is in the cluster or if no clusters are defined
    if (non->cluster==-1 || non->cluster==cl){   
      for (alpha=0;alpha<3;alpha++){
         /* Read mu(ti) */
        if (!strcmp(non->hamiltonian,"Coupling")){
          copyvec(mu_xyz+non->singles*alpha,vecr,non->singles);
        } else {
          if (read_mue(non,vecr,mu_traj,ti,alpha)!=1){
             printf("Dipole trajectory file to short, could not fill buffer!!!\n");
             printf("ITIME %d %d\n",ti,alpha);
             exit(1);
            }
        }       
          /* this is for time t1 to generate the vector for the dipole with 0 as the imagine part*/ 
          clearvec(veci,non->singles);
          /*this is just copy the input dipole to the real part of the vector*/
          /*This two step is necessary as the input dipole is real number, but after probagate it becomes complex number */
          //copyvec(vecr,mu_eg,non->singles);
          /* Loop over coherence time */

        
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
          //Here we generate the equilibrium density operator
          eq_den(Hamil_i_e,rho_l,N,non);
          // Multiply the density operator to dipole operator,vecr, as it is only the real number.
          for (a=0;a<N;a++){
            for (b=0;b<N;b++){
              mid_ver[a]+=rho_l[a+b*N]*vecr[b];
            }
          }
          // Update dipole operator
          for (a=0;a<N;a++){
            vecr[a]=mid_ver[a];
          }      

          for (beta=0;beta<3;beta++){
            /* Read mu(tj) */
            if (!strcmp(non->hamiltonian,"Coupling")){
              copyvec(mu_xyz+non->singles*beta,mu_eg,non->singles);
            } else {
              if (read_mue(non,mu_eg,mu_traj,tj,beta)!=1){
                printf("Dipole trajectory file to short, could not fill buffer!!!\n");
                printf("JTIME %d %d\n",tj,beta);
                exit(1);
              }
            }
            
            /* Here calculate doorway function/
            /* Inner product for all sites*/
            /*here we need to make clear the seg_num and the site_num
            the number of the position should be decided by the segment number,
            we should sum all the value in one segment*/

            for (site_num=0;site_num<non->singles;site_num++){
              seg_num=non->psites[site_num];
              /*this equation is make sure the calculated data in the right position */
              index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
              re_window_SE[index]+=mu_eg[site_num]*vecr[site_num];
              im_window_SE[index]+=mu_eg[site_num]*veci[site_num]; 
             }
          }  

          
          /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
          if (non->Npsites==non->singles){
            zero_coupling(Hamil_i_e,non);
          } else {
            printf("Segment number and the projection number are different");
            exit(1);
          }
          /* Propagate dipole moment */
          propagate_vector(non,Hamil_i_e,vecr,veci,1,samples,t1*alpha);
          
      }
    }
  }

    /* Update Log file with time and sample numner */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished sample %d\n",samples);        
    time_now=log_time(time_now,log);
    fclose(log);

  }
 
  /* The calculation is finished we can close all auxillary arrays before writing */
  /* output to file. */
  free(vecr);
  free(veci);
  free(mu_eg);
  free(mu_xyz);
  free(Hamil_i_e);
  free(rho_l);
  free(mid_ver);

  /* The calculation is finished, lets write output */
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Calculating Response!\n");
  fprintf(log,"Writing to file!\n");  
  fclose(log);

  /* Print information on number of realizations included belonging to the selected */
  /* cluster and close the cluster file. (Only to be done if cluster option is active.) */
  if (non->cluster!=-1){
    printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
    fclose(Cfile);
  }

  /* Close Trajectory Files */
  fclose(mu_traj),fclose(H_traj);

  /* Save  the imaginary part for time domain response */
  outone=fopen("CG_2DES_windows_SE_im.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index =  seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",im_window_SE[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }

  /* Save the real part for time domain response */
  outone=fopen("CG_2DES_windows_SE_re.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index = seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",re_window_SE[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }
 fclose(outone);
}


void CG_2DES_P_DA(t_non *non,float *P_DA);
void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_2DES,float *im_2DES);
void combine_CG_2DES(t_non *non,float *re_doorway,float *im_doorway,
    float *P_DA,float *re_window_GB,float *im_window_GB,
    float *re_window_SE,float *im_window_SE,float *re_window_EA,float *im_window_EA,
    float *re_2DES,float *im_2DES);