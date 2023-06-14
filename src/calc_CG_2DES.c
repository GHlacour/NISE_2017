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
#include "read_trajectory.h"
#include "eq_den.h"
#include "mpi.h"
#include "MPI_subs.h"


/* This is the main routine which will control the flow of the
 * coarse grained 2DES calculation. It takes the general NISE
 * input as input.
 */


void calc_CG_2DES(t_non *non){
    float *re_doorway, *im_doorway; 
    float *re_window_SE, * im_window_SE;
    float *re_window_GB, * im_window_GB;
    float *re_window_EA, * im_window_EA;
    float *re_2DES , *im_2DES;


    int pro_dim;
    pro_dim=project_dim(non);
    printf("Dimension %d\n",non->tmax*9*pro_dim);
    re_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    re_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));   
    re_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float)); 
    re_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));  
    printf("Print %d\n",non->tmax2);

  
    //re_2DES = (float *)calloc(non->tmax*non->tmax*non->tmax2,sizeof(float));
    //im_2DES = (float *)calloc(non->tmax*non->tmax*non->tmax2,sizeof(float));
    printf("Performing the CG_2DES calculation.\n"); 


    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway"))) {
        CG_2DES_doorway(non, re_doorway, im_doorway);
      }  
     if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ) {
        CG_2DES_window_SE(non, re_window_SE, im_window_SE); 
     } 
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_GB"))){
        CG_2DES_window_GB(non, re_window_GB, im_window_GB); 
      }
  
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_window_EA"))){
        CG_2DES_window_EA(non, re_window_EA, im_window_EA); 
      }

      /*if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
         (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")))  {*/
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_full_2DES_segments")))  {
        printf(" full  part");
        CG_full_2DES_segments(non,re_doorway,im_doorway,
                                  re_window_SE,im_window_SE,
                                  re_window_GB, im_window_GB,
                                  re_window_EA,im_window_EA);
     }

    free(re_doorway);
    free(im_doorway);
    free(re_window_SE);
    free(im_window_SE);
    free(re_window_GB);
    free(im_window_GB);
    free(re_window_EA);
    free(im_window_EA);  

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

  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

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

          //printf("door Hamil_i_e  "); 
          //printf("%f \n", Hamil_i_e[0]); 

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
                /*if (index>non->tmax*9*pro_dim){
                  printf("%d %d %d %d %d\n",seg_num,alpha,beta,t1,index);
                  exit(0);
                }
                if (index<0){
                  printf("Negative index\n");
                }*/
                //printf("%d %d\n",index,non->tmax*9*pro_dim);
                //printf("%f\n",re_doorway[index]);
                //printf("%f\n",im_doorway[index]);
                //printf("%d %f\n",site_num,vecr[site_num]);
                //index=0;
                re_doorway[index]+=mu_eg[site_num]*vecr[site_num];
                im_doorway[index]+=mu_eg[site_num]*veci[site_num]; 
                //printf("im_doorway[index] ");
                //printf("%f\n",im_doorway[index]);
              }
          }  

          /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
          if (non->Npsites==non->singles){
            zero_coupling(Hamil_i_e,non);
          } else {
            printf("Segment number and the projection number are different");
            exit(1);
          }       
          //printf("door way Hamil_i_e "); 
          //printf("%f \n", Hamil_i_e[0]);  
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
 printf("the doorway part run successfully \n ");  
}



//void CG_2DES_P_DA(t_non *non,float *P_DA,float* K, float* P0, int N){
void CG_2DES_P_DA(t_non *non,float *P_DA,float* K, int N){
    printf("Calculate population transfer!\n");
    float *eigK_re, *eigK_im; // eigenvalues of K
    float *evecL, *evecR; // eigenvectors of K
    float *ivecR, *ivecL; //inverse eigenvectors of K
    float *iP0;
    float *cnr;
    int a, b, c;
    FILE *outone;
    
    eigK_re = (float *)calloc(N*N,sizeof(float));
    eigK_im = (float *)calloc(N*N,sizeof(float));
    evecL = (float *)calloc(N*N,sizeof(float));
    evecR = (float *)calloc(N*N,sizeof(float));
    ivecL = (float *)calloc(N*N,sizeof(float));
    ivecR = (float *)calloc(N*N,sizeof(float));
    iP0 = (float *)calloc(N*N,sizeof(float));

    // Diagonalize K matrix
    diagonalize_real_nonsym(K, eigK_re, eigK_im, evecL, evecR, ivecL, ivecR, N);
    for (int a = 0; a<N; a++) {
        if (eigK_im[a]!=0) {
            printf("Transfer lifetime is not real!\n");
            exit(0);
        }
    }

    // Calculate P(t2) = expm(-K*t2)*P(0) = evecR*exp(Eig*t2)*ivecR*P(0)
    //Here P0 will be input outside
    /*for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                iP0[a + c * N] += ivecR[a + b * N] * P0[b + c * N];
            }
        }
    }*/   
    // iP0 = ivecR*P(0)

    // Loop over t2
    /*Here we assume the P0 is a N*N matrix*/
    for (int nt2 = 0; nt2<non->tmax2; nt2++) {
        cnr = (float *)calloc(N * N, sizeof(float));
        for (a = 0; a < N; a++) {
            for (b = 0; b < N; b++) {
                cnr[a + b * N] += exp(eigK_re[a] * nt2 * non->deltat) * ivecR[a + b*N];
            }
        }   // exp(Eig*t2)*iP0

        for (a = 0; a < N; a++) {
            for (b = 0; b < N; b++) {
                for (c = 0; c < N; c++) {
                    //P_DA[nt2 + (a + c * N)*non->tmax2] += evecR[a + b * N] * cnr[b + c * N];
                    P_DA[nt2*N+c*N*non->tmax2+a] += evecR[a + b * N] * cnr[b + c * N];
                }
            }
        }   // evecR*cnr
        free(cnr);
    }

    // Write to file
    outone=fopen("KPop.dat","w");
    for (int t2=0;t2<non->tmax2;t2+=non->dt2){
        fprintf(outone,"%f ",t2*non->deltat);
        for (int a=0;a<N;a++){
            for (int b=0;b<N;b++){
                fprintf(outone,"%f ",P_DA[t2*N+b*N*non->tmax2+a]);
                                          
            }
        }
        fprintf(outone,"\n"); 
    }
    fclose(outone);

    free(eigK_im);
    free(eigK_re);
    free(iP0);
    free(evecL);
    free(evecR);
    free(ivecL);
    free(ivecR);

    return;
    printf("the PD_t2 part run successfully \n ");  
}




void CG_2DES_window_SE(t_non *non, float *re_window_SE, float *im_window_SE){
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
  open_files(non,&H_traj,&mu_traj,&Cfile);

 Ncl=0; /* Counter for snapshots calculated*/

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


  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

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
	  eq_den(Hamil_i_e,rho_l,N,non);
          /*printf("Hamil_i_e  ");
          printf("%f \n", Hamil_i_e); */
          // Multiply the density operator to dipole operator,vecr, as it is only the real number.
          clearvec(mid_ver,non->singles);
          for (a=0;a<N;a++){
            for (b=0;b<N;b++){
              mid_ver[a]+=rho_l[a+b*N]*vecr[b];
            }
          }
          // Update dipole operator
          for (a=0;a<N;a++){
            vecr[a]=mid_ver[a];
          }


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
          propagate_vector(non,Hamil_i_e,vecr,veci,-1,samples,t1*alpha);
          
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
 printf("the SE part run successfully \n ");  
}


void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB){
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
  /* projection */
  pro_dim=project_dim(non);
  /* check projection */   
  if (non->Npsites!=non->singles){
      printf("Input segment number and the projection segment number are different.\n");
      printf("please check the input file or the projection file\n");
      exit(1);
  }

  /* Allocate memory*/  /*(tmax+1)*9*/
  nn2=non->singles*(non->singles+1)/2;
  Hamil_i_e=(float *)calloc(nn2,sizeof(float));

  /* Open Trajectory files */
  open_files(non,&H_traj,&mu_traj,&Cfile);
 Ncl=0; /* Counter for snapshots calculated*/
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


  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

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
          //eq_den(Hamil_i_e,rho_l,N,non);
          // Multiply the density operator to dipole operator,vecr, as it is only the real number.
          /*for (a=0;a<N;a++){
            for (b=0;b<N;b++){
              mid_ver[a]+=rho_l[a+b*N]*vecr[b];
            }
          }
          // Update dipole operator
          for (a=0;a<N;a++){
            vecr[a]=mid_ver[a];
          } */
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
              re_window_GB[index]+=mu_eg[site_num]*vecr[site_num];
              im_window_GB[index]+=mu_eg[site_num]*veci[site_num]; 
              //printf("re_window_GB "); 
              //printf("%f \n", re_window_GB); 
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
          propagate_vector(non,Hamil_i_e,vecr,veci,-1,samples,t1*alpha);
          
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
  outone=fopen("CG_2DES_windows_GB_im.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index =  seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",im_window_GB[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }

  /* Save the real part for time domain response */
  outone=fopen("CG_2DES_windows_GB_re.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){  
                index = seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",re_window_GB[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }
 fclose(outone);
 printf("the GB part run successfully \n ");  
}

void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA){
   /* Initialize variables*/
 /* The window part for SE*/
  float *Hamil_i_e,*Hamil_i_ee;
  float *mu_eg;
  float *vecr, *vecr1, *veci, *rho_i;   
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
  int alpha,beta,a_b,ti,tj,tk,i,j; /*alpha and beta are corresponding the dipole for different segment,  a_b  = alpha*beta, which used in the writing file*/
  int t1,t2;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;
  int a,b,index,seg,site_num,seg_num;/*site num is the number of the site, which used to make sure the index number later*/
// reserve memory for the density operator
  int N;
  N = non->singles;
  float *rho_l;
   /* Allocate memory*/  /*(tmax+1)*9*/
  nn2=non->singles*(non->singles+1)/2;
// reserve memory for transition dipole from first to second
  float *fr;
  float *fi;
  //here the mid_vcr is just used  as the mid part for multiply the dipole with density operator
//  float *mid_ver;
//  mid_ver = (float *)calloc(N,sizeof(float));

  /* projection */
  pro_dim=project_dim(non);
//Here we reserve memory for matrix multiplication
  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;

// reserve memory for transition dipole
  vecr=(float *)calloc(N*N,sizeof(float));	
  vecr1=(float *)calloc(N,sizeof(float));	
  veci=(float *)calloc(N*N,sizeof(float));
  mu_eg=(float *)calloc(N,sizeof(float));
  mu_xyz=(float *)calloc(N*3,sizeof(float));
  Hamil_i_e=  (float *)calloc(nn2,sizeof(float));
  Hamil_i_ee= (float *)calloc(nn2,sizeof(float));
  fr  = calloc(nn2*N, sizeof(float));
  fi  = calloc(nn2*N, sizeof(float));
  rho_i = (float *)calloc(N*N,sizeof(float));
  rho_l=(float *)calloc(N*N,sizeof(float));

  /* check projection */   
  if (non->Npsites!=non->singles){
      printf("Input segment number and the projection segment number are different.\n");
      printf("please check the input file or the projection file\n");
      exit(1);
  }

  /* Open Trajectory files */
  open_files(non,&H_traj,&mu_traj,&Cfile);
  Ncl=0; /* Counter for snapshots calculated*/
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


  /* Read coupling, this is done if the coupling and transition-dipoles are */
  /* time-independent and only one snapshot is stored */
  read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);
  
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
          copyvec(mu_xyz+non->singles*alpha,vecr1,non->singles);
        } else {
          if (read_mue(non,vecr1,mu_traj,ti,alpha)!=1){
             printf("Dipole trajectory file to short, could not fill buffer!!!\n");
             printf("ITIME %d %d\n",ti,alpha);
             exit(1);
            }
        }
        clearvec(rho_i,non->singles*non->singles); 
          //Here we generate the equilibrium density operator
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
        eq_den(Hamil_i_e,rho_l,N,non);       
        
        for (a=0;a<N;a++){
          dipole_double_CG2DES(non, vecr1, rho_l+a*N, rho_i+a*N, fr+a*nn2, fi+a*nn2); 
        }


          /* this is for time t1 to generate the vector for the dipole with 0 as the imagine part*/ 
          /*this is just copy the input dipole to the real part of the vector*/
          /*This two step is necessary as the input dipole is real number, but after probagate it becomes complex number */
          //copyvec(vecr,mu_eg,non->singles);
          /* Loop over coherence time */
        
         for (t1=0;t1<non->tmax;t1++){
            tj=ti+t1;
	    //printf("tj %d\n",tj);
          
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
            //dipole_double_CG2DES(non, mu_eg, vecr, veci, fr, fi);            
            /* Here calculate doorway function/
            /* Inner product for all sites*/
            /*here we need to make clear the seg_num and the site_num
            the number of the position should be decided by the segment number,
            we should sum all the value in one segment*/
            
            /* Multiply with mu_ef */
            /* Take fr and fi back to vecr and veci */
            clearvec(vecr,non->singles*non->singles);
            clearvec(veci,non->singles*non->singles);
            for (a=0;a<N;a++){
              //The dimention of vecr is only   vecr=(float *)calloc(N,sizeof(float));
                dipole_double_inverse_CG2DES(non, mu_eg, fr+a*nn2, fi+a*nn2, vecr+a*N, veci+a*N);
                //dipole_double_inverse_CG2DES(non, mu_eg, fr+a*nn2, fi+a*nn2, vecr, veci);
            }
	    
            /* Propagate back to time ti； tj=ti+t1; */
            for (tk=tj;tk>ti;tk--){
	     // printf("ti tj tk: %d %d %d\n",ti,tj,tk);
              /* Read Hamiltonian */
              /* Read hamiltonian at tk */
              if (!strcmp(non->hamiltonian,"Coupling")){
                if (read_Dia(non,Hamil_i_ee,H_traj,tk)!=1){
                  printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                  exit(1);
                }
              } else {
                if (read_He(non,Hamil_i_ee,H_traj,tk)!=1){
                printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                exit(1);
              }
              }

              /* Zero coupling in hamiltonian */
              /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
                zero_coupling(Hamil_i_ee,non);
	      
              /* Propagate single excited states vecr and veci backwards with propagate */
              /* Propagate dipole moment */
            for (a=0;a<N;a++){
               propagate_vector(non,Hamil_i_ee,vecr+a*N,veci+a*N,1,samples,tk*alpha);
	    //propagate_matrix(non,Hamil_i_ee,vecr,veci,1,samples,tk*alpha);
            }

            }

            /* Calculate window function by taking trace of vecr and veci */
            for (site_num=0;site_num<non->singles;site_num++){
              seg_num=non->psites[site_num];
              /*this equation is make sure the calculated data in the right position */
              index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
              re_window_EA[index]+=vecr[site_num*N+site_num];
              im_window_EA[index]+=veci[site_num*N+site_num]; 
              //printf("re_window_EA[index]");
              //printf("%d \n", index);
             }
          }  
          
          /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
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
            zero_coupling(Hamil_i_e,non);

          /* Propagate dipole moment */
          for (a=0;a<N;a++){
            propagate_vec_coupling_S_doubles_ES(non, Hamil_i_e, fr+a*nn2, fi+a*nn2, non->ts);   
          }
	 
	 }
	// exit(0);
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
  printf("this is a test for EA \n");  

  free(mu_xyz);
  free(Hamil_i_e);
  free(Hamil_i_ee);
  //free(mid_ver);
  free(vecr);
  free(veci);
  free(vecr1);
  free(fr);
  free(fi);
  free(rho_i);
  free(rho_l);
  free(mu_eg);

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
  printf("write EA IM ");
  /* Save  the imaginary part for time domain response */
  outone=fopen("CG_2DES_windows_EA_im.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){

      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){
                index =  seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",im_window_EA[index]/samples);
       }
     }
    }

  }
  printf("write EA RE ");
  fprintf(outone,"\n"); 
  /* Save the real part for time domain response */
  outone=fopen("CG_2DES_windows_EA_re.dat","w");
  for (t1=0;t1<non->tmax1;t1+=non->dt1){
    fprintf(outone,"%f ",t1*non->deltat);
    for (seg_num=0;seg_num<pro_dim;seg_num++){
      for (alpha=0;alpha<3;alpha++){
              for (beta=0;beta<3;beta++){  
                index = seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                //fprintf(outone,"%e %e ",re_doorway[index]/samples,im_doorway[index]/samples);
                fprintf(outone,"%e ",re_window_EA[index]/samples);
       }
     }
    }
  fprintf(outone,"\n"); 
  }
 fclose(outone);

 
    //printf("this is a test \n"); 
  //exit(1);
  printf("the EA part run successfully \n ");  
}


void CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
                                      float *re_window_SE,float *im_window_SE,
                                      float *re_window_GB, float *im_window_GB,
                                      float *re_window_EA,float *im_window_EA){
  float  *P_DA, *K, *PDA_t2;
  float *int_sna_t1_re,*int_sna_t1_im_NR,*int_sna_t1_im_R; 
  float *int_sna_t3_SE_re,*int_sna_t3_SE_im;
  float *int_sna_t3_GB_re,*int_sna_t3_GB_im;
  float *int_sna_t3_EA_re,*int_sna_t3_EA_im;
  float  *int_PDA;
  float up_ver1_re, up_ver1_im_NR,up_ver1_im_R;
  float *up_ver2_re,*up_ver2_im_NR,*up_ver2_im_R;
  float up_ver3_SE_re, up_ver3_SE_im_NR,up_ver3_SE_im_R;
  float up_ver3_GB_re, up_ver3_GB_im_NR,up_ver3_GB_im_R;
  float up_ver3_EA_re, up_ver3_EA_im_NR,up_ver3_EA_im_R;
  float *re_2DES_pa, *im_2DES_NR_pa, *im_2DES_R_pa;
  float *re_2DES_pe, *im_2DES_NR_pe, *im_2DES_R_pe;
  float *re_2DES_cr, *im_2DES_NR_cr, *im_2DES_R_cr;
  int pro_dim;
  int t1, t2, t3 ;
  int a,b,c,d,e,f;
  int seg_num_t1, seg_num_t3;
  FILE *outone,*log, *Rate;
  pro_dim=project_dim(non);
  //P0 =(float *)calloc(pro_dim*pro_dim,sizeof(float));
  K =(float *)calloc(pro_dim*pro_dim,sizeof(float));
  PDA_t2 = (float *)calloc(pro_dim*pro_dim,sizeof(float));
  P_DA = (float *)calloc(pro_dim*pro_dim*non->tmax2,sizeof(float)); 
  int_sna_t1_re = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t1_im_NR = (float *)calloc(pro_dim,sizeof(float));
  int_sna_t1_im_R =  (float *)calloc(pro_dim,sizeof(float));
  int_sna_t3_SE_re = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t3_SE_im = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t3_GB_re = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t3_GB_im = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t3_EA_re = (float *)calloc(pro_dim,sizeof(float)); 
  int_sna_t3_EA_im = (float *)calloc(pro_dim,sizeof(float)); 
  int_PDA = (float *)calloc(pro_dim,sizeof(float));
  up_ver2_re = (float *)calloc(pro_dim,sizeof(float));
  up_ver2_im_NR = (float *)calloc(pro_dim,sizeof(float));
  up_ver2_im_R = (float *)calloc(pro_dim,sizeof(float));



  re_2DES_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_NR_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_R_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  re_2DES_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_NR_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_R_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  re_2DES_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_NR_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  im_2DES_R_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
  int index;
  int sampleCount;
  sampleCount = 3; // which numbe should I set this number 
  float **re_2DES_pa_sum, **im_2DES_NR_pa_sum, **im_2DES_R_pa_sum;
  float **re_2DES_pe_sum, **im_2DES_NR_pe_sum, **im_2DES_R_pe_sum;
   float **re_2DES_cr_sum, **im_2DES_NR_cr_sum, **im_2DES_R_cr_sum;

  //re_2DES_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_NR_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_R_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //re_2DES_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_NR_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_R_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //re_2DES_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_NR_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));
  //im_2DES_R_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1*non->tmax2,sizeof(float),sizeof(float*));


  re_2DES_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_NR_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_R_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  re_2DES_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_NR_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_R_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  re_2DES_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_NR_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
  im_2DES_R_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));





  int molpol, px[4]; 
  polar(px,molpol);


  /*Read the rate matrix*/
  Rate=fopen("RateMatrix.dat","r");
  if (Rate==NULL){
    printf("RateMatrix file not found!\n");
    exit(1);
  }
  //printf("%f ", Rate[0]);
  //printf("%f ", Rate[2]);

  for (a=0;a<pro_dim*pro_dim;a++){
    if (fscanf(Rate,"%f",&K[a])!=1){
      printf("Error in reading in ratematrix!\n");
      exit(0);
    }
    printf("%f ",K[a]);
  }

  CG_2DES_P_DA(non,P_DA, K, pro_dim);
  /* tmax and tmax1,tmax2,tmax3 are essential the same */
  /*This part is calculate xxyy,yyxx,zzxx , and xxxx yyyy zzzz which in total 9*/
  /*for (t1=0;t1<non->tmax;t1++)*/



  for (t1=0; t1<non->tmax1; t1+=1){
    //printf("%d" ,  t1);
    for (a=0;a<9;a+=4){
      for (seg_num_t1=0;seg_num_t1<pro_dim;seg_num_t1++){
        int_sna_t1_re[seg_num_t1] = re_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_NR[seg_num_t1] = im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_R[seg_num_t1] = -1*im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1]; 

      }  
      for (t2=0; t2<non->tmax2; t2+=1){
            for (int a=0;a<pro_dim;a++){
                for (int b=0;b<pro_dim;b++){
                  PDA_t2[pro_dim*b+a]=P_DA[t2*pro_dim+pro_dim*b*non->tmax2+a];

                  //P_DA[nt2*N+c*N*non->tmax2+a] += evecR[a + b * N] * cnr[b + c * N];
              
              }
            }
          
            /*dimention 1*N* N*N *N*1*/
            /*first calculate the part1: dimention: 1*N* N*N */
            for (c=0;c<pro_dim;c++){
              up_ver1_re=0; 
              up_ver1_im_NR=0;
              up_ver1_im_R=0;
                 
              for (d=0;d<pro_dim;d++){
                int_PDA[d] = PDA_t2[pro_dim*d+c];
              }
              for  (e=0;e<pro_dim;e++){
                up_ver1_re += int_sna_t1_re[e]*int_PDA[e];
                up_ver1_im_NR+=int_sna_t1_im_NR[e]*int_PDA[e];
                up_ver1_im_R+=int_sna_t1_im_R[e]*int_PDA[e];
              }
            up_ver2_re[c]=up_ver1_re;
            up_ver2_im_NR[c]=up_ver1_im_NR;
            up_ver2_im_R[c]=up_ver1_im_R;
          }  


          for (t3=0; t3<non->tmax1; t3+=1){
              up_ver3_SE_re=0;
              up_ver3_SE_im_NR=0;
              up_ver3_SE_im_R=0;
              up_ver3_GB_re=0;
              up_ver3_GB_im_NR=0;
              up_ver3_GB_im_R=0;
              up_ver3_EA_re=0;
              up_ver3_EA_im_NR=0;
              up_ver3_EA_im_R=0;  
              for (b=0;b<9;b+=4){
                  for (seg_num_t3=0;seg_num_t3<pro_dim;seg_num_t3++){
                    int_sna_t3_SE_re[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_SE_im[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_re[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_im[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_re[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_im[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];     
                  //printf("im");
                  //printf("%f" , int_sna_t3_SE_im[seg_num_t3]);
                  //printf("re");
                  //printf("%f" ,  int_sna_t3_SE_re[seg_num_t3]);



                  }
   

              /*second calculate the part2: dimention: 1*N* N*1 */    
                for  (f=0;f<pro_dim;f++){
                  up_ver3_SE_re += up_ver2_re[f]*int_sna_t3_SE_re[f];
                  up_ver3_SE_im_NR += up_ver2_im_NR[f]*int_sna_t3_SE_im[f];
                  up_ver3_SE_im_R += up_ver2_im_R[f]*int_sna_t3_SE_im[f];
                  up_ver3_GB_re += up_ver2_re[f]*int_sna_t3_GB_re[f];
                  up_ver3_GB_im_NR += up_ver2_im_NR[f]*int_sna_t3_GB_im[f];
                  up_ver3_GB_im_R += up_ver2_im_R[f]*int_sna_t3_GB_im[f];
                  up_ver3_EA_re += up_ver2_re[f]*int_sna_t3_EA_re[f];
                  up_ver3_EA_im_NR += up_ver2_im_NR[f]*int_sna_t3_EA_im[f];
                  up_ver3_EA_im_R += up_ver2_im_R[f]*int_sna_t3_EA_im[f];                  
                }
          if (a==b){
            /* Parallel polarization:*/
            float polWeight = polarweight(0, 0);
            //printf("test part  ");
            //printf("%f ",re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3] );
            //printf("\n");
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 0);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 0);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;


          } else{
           /* Parallel polarization:*/
            float polWeight = polarweight(0, 3);

            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 1);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 3);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            }
           }
        }
      }
    }
  }


  /*This part calculate the xyxy,yxyx and yxxy,xyyx*/
  for (t1=0; t1<non->tmax1; t1+=1){
    for (a=1;a<4;a+=2){
      for (seg_num_t1=0;seg_num_t1<pro_dim;seg_num_t1++){
        int_sna_t1_re[seg_num_t1] = re_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_NR[seg_num_t1] = im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_R[seg_num_t1] = -1*im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];  
      }  
      for (t2=0; t2<non->tmax2; t2+=1){
            for (int a=0;a<pro_dim;a++){
                for (int b=0;b<pro_dim;b++){
                  PDA_t2[pro_dim*b+a]=P_DA[t2*pro_dim+pro_dim*b*non->tmax2+a];
                  //P_DA[nt2*N+c*N*non->tmax2+a] += evecR[a + b * N] * cnr[b + c * N];
              
              }
            }
            /*dimention 1*N* N*N *N*1*/
            /*first calculate the part1: dimention: 1*N* N*N */
            for (c=0;c<pro_dim;c++){
              up_ver1_re=0; 
              up_ver1_im_NR=0;
              up_ver1_im_R=0;
                 
              for (d=0;d<pro_dim;d++){
                int_PDA[d] = PDA_t2[pro_dim*d+c];
              }
              for  (e=0;e<pro_dim;e++){
                up_ver1_re += int_sna_t1_re[e]*int_PDA[e];
                up_ver1_im_NR+=int_sna_t1_im_NR[e]*int_PDA[e];
                up_ver1_im_R+=int_sna_t1_im_R[e]*int_PDA[e];
              }
            up_ver2_re[c]=up_ver1_re;
            up_ver2_im_NR[c]=up_ver1_im_NR;
            up_ver2_im_R[c]=up_ver1_im_R;
          }  
          for (t3=0; t3<non->tmax1; t3+=1){
              up_ver3_SE_re=0;
              up_ver3_SE_im_NR=0;
              up_ver3_SE_im_R=0;
              up_ver3_GB_re=0;
              up_ver3_GB_im_NR=0;
              up_ver3_GB_im_R=0;
              up_ver3_EA_re=0;
              up_ver3_EA_im_NR=0;
              up_ver3_EA_im_R=0;  
              for (b=1;b<4;b+=2){
                  for (seg_num_t3=0;seg_num_t3<pro_dim;seg_num_t3++){
                    int_sna_t3_SE_re[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_SE_im[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_re[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_im[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_re[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_im[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];                    
                  }
              /*second calculate the part2: dimention: 1*N* N*1 */    
                for  (f=0;f<pro_dim;f++){
                  up_ver3_SE_re += up_ver2_re[f]*int_sna_t3_SE_re[f];
                  up_ver3_SE_im_NR += up_ver2_im_NR[f]*int_sna_t3_SE_im[f];
                  up_ver3_SE_im_R += up_ver2_im_R[f]*int_sna_t3_SE_im[f];
                  up_ver3_GB_re += up_ver2_re[f]*int_sna_t3_GB_re[f];
                  up_ver3_GB_im_NR += up_ver2_im_NR[f]*int_sna_t3_GB_im[f];
                  up_ver3_GB_im_R += up_ver2_im_R[f]*int_sna_t3_GB_im[f];
                  up_ver3_EA_re += up_ver2_re[f]*int_sna_t3_EA_re[f];
                  up_ver3_EA_im_NR += up_ver2_im_NR[f]*int_sna_t3_EA_im[f];
                  up_ver3_EA_im_R += up_ver2_im_R[f]*int_sna_t3_EA_im[f];                  
                }
          if (a==b ){
            /*Here only for xyxy*/
            /* Parallel polarization:*/
           float polWeight = polarweight(0, 9);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 9);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 9);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;


          } else{
            /*Here only for xyyx*/
           /* Parallel polarization:*/
           float polWeight = polarweight(0, 15);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 15);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 15);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            }
           }
        }
      }
    }
  }


   /*This part calculate the xzxz zxzx and zxxz xzzx*/

  for (t1=0; t1<non->tmax1; t1+=1){
    for (a=2;a<7;a+=4){
      for (seg_num_t1=0;seg_num_t1<pro_dim;seg_num_t1++){
        int_sna_t1_re[seg_num_t1] = re_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_NR[seg_num_t1] = im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_R[seg_num_t1] = -1*im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];  
      }  
      for (t2=0; t2<non->tmax2; t2+=1){
            for (int a=0;a<pro_dim;a++){
                for (int b=0;b<pro_dim;b++){
                  PDA_t2[pro_dim*b+a]=P_DA[t2*pro_dim+pro_dim*b*non->tmax2+a];
                  //P_DA[nt2*N+c*N*non->tmax2+a] += evecR[a + b * N] * cnr[b + c * N];
              
              }
            }
            /*dimention 1*N* N*N *N*1*/
            /*first calculate the part1: dimention: 1*N* N*N */
            for (c=0;c<pro_dim;c++){
              up_ver1_re=0; 
              up_ver1_im_NR=0;
              up_ver1_im_R=0;
                 
              for (d=0;d<pro_dim;d++){
                int_PDA[d] = PDA_t2[pro_dim*d+c];
              }
              for  (e=0;e<pro_dim;e++){
                up_ver1_re += int_sna_t1_re[e]*int_PDA[e];
                up_ver1_im_NR+=int_sna_t1_im_NR[e]*int_PDA[e];
                up_ver1_im_R+=int_sna_t1_im_R[e]*int_PDA[e];
              }
            up_ver2_re[c]=up_ver1_re;
            up_ver2_im_NR[c]=up_ver1_im_NR;
            up_ver2_im_R[c]=up_ver1_im_R;
          }  
          for (t3=0; t3<non->tmax1; t3+=1){
              up_ver3_SE_re=0;
              up_ver3_SE_im_NR=0;
              up_ver3_SE_im_R=0;
              up_ver3_GB_re=0;
              up_ver3_GB_im_NR=0;
              up_ver3_GB_im_R=0;
              up_ver3_EA_re=0;
              up_ver3_EA_im_NR=0;
              up_ver3_EA_im_R=0;  
              for (b=2;b<7;b+=4){
                  for (seg_num_t3=0;seg_num_t3<pro_dim;seg_num_t3++){
                    int_sna_t3_SE_re[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_SE_im[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_re[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_im[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_re[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_im[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];                    
                  }
  
              /*second calculate the part2: dimention: 1*N* N*1 */    
                for  (f=0;f<pro_dim;f++){
                  up_ver3_SE_re += up_ver2_re[f]*int_sna_t3_SE_re[f];
                  up_ver3_SE_im_NR += up_ver2_im_NR[f]*int_sna_t3_SE_im[f];
                  up_ver3_SE_im_R += up_ver2_im_R[f]*int_sna_t3_SE_im[f];
                  up_ver3_GB_re += up_ver2_re[f]*int_sna_t3_GB_re[f];
                  up_ver3_GB_im_NR += up_ver2_im_NR[f]*int_sna_t3_GB_im[f];
                  up_ver3_GB_im_R += up_ver2_im_R[f]*int_sna_t3_GB_im[f];
                  up_ver3_EA_re += up_ver2_re[f]*int_sna_t3_EA_re[f];
                  up_ver3_EA_im_NR += up_ver2_im_NR[f]*int_sna_t3_EA_im[f];
                  up_ver3_EA_im_R += up_ver2_im_R[f]*int_sna_t3_EA_im[f];                  
                }
          if (a==b ){
            /*Here only for xzxz and zxzx*/
            /* Parallel polarization:*/
           float polWeight = polarweight(0, 9);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 9);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 9);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;


          } else{
            /*Here only for zxxz xzzx*/
           /* Parallel polarization:*/


         float polWeight = polarweight(0, 15);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 15);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 15);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            }
           }
        }
      }
    }
  }

 /*This part calculate the yzyz zyzy  and zyyz yzzy*/
  for (t1=0; t1<non->tmax1; t1+=1){
    for (a=5;a<8;a+=2){
      for (seg_num_t1=0;seg_num_t1<pro_dim;seg_num_t1++){
        int_sna_t1_re[seg_num_t1] = re_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_NR[seg_num_t1] = im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];
        int_sna_t1_im_R[seg_num_t1] = -1*im_doorway[ seg_num_t1*9*non->tmax+a*non->tmax+t1];  
      }  
      for (t2=0; t2<non->tmax2; t2+=1){
            for (int a=0;a<pro_dim;a++){
                for (int b=0;b<pro_dim;b++){
                  PDA_t2[pro_dim*b+a]=P_DA[t2*pro_dim+pro_dim*b*non->tmax2+a];
                  //P_DA[nt2*N+c*N*non->tmax2+a] += evecR[a + b * N] * cnr[b + c * N];
              
              }
            }
            /*dimention 1*N* N*N *N*1*/
            /*first calculate the part1: dimention: 1*N* N*N */
            for (c=0;c<pro_dim;c++){
              up_ver1_re=0; 
              up_ver1_im_NR=0;
              up_ver1_im_R=0;
                 
              for (d=0;d<pro_dim;d++){
                int_PDA[d] = PDA_t2[pro_dim*d+c];
              }
              for  (e=0;e<pro_dim;e++){
                up_ver1_re += int_sna_t1_re[e]*int_PDA[e];
                up_ver1_im_NR+=int_sna_t1_im_NR[e]*int_PDA[e];
                up_ver1_im_R+=int_sna_t1_im_R[e]*int_PDA[e];
              }
            up_ver2_re[c]=up_ver1_re;
            up_ver2_im_NR[c]=up_ver1_im_NR;
            up_ver2_im_R[c]=up_ver1_im_R;
          }  


          for (t3=0; t3<non->tmax1; t3+=1){
              up_ver3_SE_re=0;
              up_ver3_SE_im_NR=0;
              up_ver3_SE_im_R=0;
              up_ver3_GB_re=0;
              up_ver3_GB_im_NR=0;
              up_ver3_GB_im_R=0;
              up_ver3_EA_re=0;
              up_ver3_EA_im_NR=0;
              up_ver3_EA_im_R=0;  
              for (b=5;b<8;b+=2){
                  for (seg_num_t3=0;seg_num_t3<pro_dim;seg_num_t3++){
                    int_sna_t3_SE_re[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_SE_im[seg_num_t3] = re_window_SE[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_re[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_GB_im[seg_num_t3] = re_window_GB[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_re[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];
                    int_sna_t3_EA_im[seg_num_t3] = re_window_EA[seg_num_t3*9*non->tmax+b*non->tmax+t3];                    
                  }
  
              /*second calculate the part2: dimention: 1*N* N*1 */    
                for  (f=0;f<pro_dim;f++){
                  up_ver3_SE_re += up_ver2_re[f]*int_sna_t3_SE_re[f];
                  up_ver3_SE_im_NR += up_ver2_im_NR[f]*int_sna_t3_SE_im[f];
                  up_ver3_SE_im_R += up_ver2_im_R[f]*int_sna_t3_SE_im[f];
                  up_ver3_GB_re += up_ver2_re[f]*int_sna_t3_GB_re[f];
                  up_ver3_GB_im_NR += up_ver2_im_NR[f]*int_sna_t3_GB_im[f];
                  up_ver3_GB_im_R += up_ver2_im_R[f]*int_sna_t3_GB_im[f];
                  up_ver3_EA_re += up_ver2_re[f]*int_sna_t3_EA_re[f];
                  up_ver3_EA_im_NR += up_ver2_im_NR[f]*int_sna_t3_EA_im[f];
                  up_ver3_EA_im_R += up_ver2_im_R[f]*int_sna_t3_EA_im[f];                  
                }
          if (a==b ){
            /*Here only for yzyz zyzy*/
            /* Parallel polarization:*/
        float polWeight = polarweight(0, 9);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 9);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 9);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
          } else{
            /*Here only for zyyz yzzy*/
           /* Parallel polarization:*/

         float polWeight = polarweight(0, 15);
            re_2DES_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pa[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pa[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pa[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
            /*Perpendicular polarization*/
            polWeight = polarweight(1, 15);

            re_2DES_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_pe[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_pe[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_pe[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;
           
            /*Cross polarization*/
            polWeight = polarweight(2, 15);
            re_2DES_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_re*polWeight;
            re_2DES_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_re*polWeight;
            re_2DES_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_re*polWeight;
            im_2DES_NR_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_NR*polWeight;
            im_2DES_NR_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_NR*polWeight;
            im_2DES_NR_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_NR*polWeight;
            im_2DES_R_cr[t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_SE_im_R*polWeight;
            im_2DES_R_cr[non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_GB_im_R*polWeight;
            im_2DES_R_cr[2*non->tmax1*non->tmax2*non->tmax3+t1+t3*non->tmax1+t2*non->tmax1*non->tmax3]+=up_ver3_EA_im_R*polWeight;


            }
           }
        }


      }
    }
  }




  for (t1=0; t1<non->tmax1; t1+=1){
    //for (t2=0; t2<non->tmax2; t2+=1){
      t2 = 0;
      for (t3=0; t3<non->tmax3; t3+=1){
        //a represent GB SE EA 
        for (a=0; a<3; a+=1){
          index = t1+t3*non->tmax1+t2*non->tmax3*non->tmax1;
          //exit(0);
          //printf("%d" ,index);
          //re_2DES_pa_sum[t3][t1+t2*non->tmax1] +=re_2DES_pa[a*index];
          //im_2DES_NR_pa_sum[t3][t1+t2*non->tmax1]  +=im_2DES_NR_pa[a*index];
          //im_2DES_R_pa_sum[t3][t1+t2*non->tmax1]  +=im_2DES_R_pa[a*index];
          //re_2DES_pe_sum[t3][t1+t2*non->tmax1]  +=re_2DES_pe[a*index];
          //im_2DES_NR_pe_sum[t3][t1+t2*non->tmax1]  +=im_2DES_NR_pe[a*index];
          //im_2DES_R_pe_sum[t3][t1+t2*non->tmax1]  +=im_2DES_R_pe[a*index];
          //re_2DES_cr_sum[t3][t1+t2*non->tmax1]  +=re_2DES_cr[a*index];
          //im_2DES_NR_cr_sum[t3][t1+t2*non->tmax1]  +=im_2DES_NR_cr[a*index];
          //im_2DES_R_cr_sum[t3][t1+t2*non->tmax1]  +=im_2DES_R_cr[a*index];
          //printf("%f",re_2DES_pa[a*index]);
          //printf("\n");
          re_2DES_pa_sum[t3][t1]     +=re_2DES_pa[a*index];
          im_2DES_NR_pa_sum[t3][t1]  +=im_2DES_NR_pa[a*index];
          im_2DES_R_pa_sum[t3][t1]   +=im_2DES_R_pa[a*index];
          re_2DES_pe_sum[t3][t1]     +=re_2DES_pe[a*index];
          im_2DES_NR_pe_sum[t3][t1]  +=im_2DES_NR_pe[a*index];
          im_2DES_R_pe_sum[t3][t1]   +=im_2DES_R_pe[a*index];
          re_2DES_cr_sum[t3][t1]     +=re_2DES_cr[a*index];
          im_2DES_NR_cr_sum[t3][t1]  +=im_2DES_NR_cr[a*index];
          im_2DES_R_cr_sum[t3][t1]   +=im_2DES_R_cr[a*index];

        }
      //}
    }
  }
  /* Print 2D */

  print2D("RparII.dat", re_2DES_pa_sum, im_2DES_NR_pa_sum, non, sampleCount);
  print2D("RparI.dat",  re_2DES_pa_sum, im_2DES_R_pa_sum,  non, sampleCount);
  print2D("RperII.dat", re_2DES_pe_sum, im_2DES_NR_pe_sum, non, sampleCount);
  print2D("RperI.dat",  re_2DES_pe_sum, im_2DES_R_pe_sum,  non, sampleCount);
  print2D("RcroII.dat", re_2DES_cr_sum, im_2DES_NR_cr_sum, non, sampleCount);
  print2D("RcroI.dat",  re_2DES_cr_sum, im_2DES_R_cr_sum,  non, sampleCount);
  printf("----------------------------------------\n");
  printf(" 2DES calculation succesfully completed\n");


  free(re_2DES_pa_sum) ;
  free(im_2DES_NR_pa_sum) ;
  free(im_2DES_R_pa_sum) ;
  free(re_2DES_pe_sum);
  free(im_2DES_NR_pe_sum) ;
  free(im_2DES_R_pe_sum);
  free(re_2DES_cr_sum);
  free(im_2DES_NR_cr_sum);
  free(im_2DES_R_cr_sum);




  free(int_sna_t1_re),  free(int_sna_t1_im_NR),  free(int_sna_t1_im_R);
  free(int_sna_t3_SE_re),  free(int_sna_t3_SE_im), free(int_sna_t3_GB_re);   
  free(int_sna_t3_GB_im) , free(int_sna_t3_EA_re),  free(int_sna_t3_EA_im); 
  free(PDA_t2),  free(int_PDA), free(up_ver2_re);
  free(up_ver2_im_NR),  free(up_ver2_im_R);
  free(re_2DES_pa) ,free(im_2DES_NR_pa),  free(im_2DES_R_pa);
  free(re_2DES_pe),  free(im_2DES_NR_pe) ,free(im_2DES_R_pe);
  free(re_2DES_cr),  free(im_2DES_NR_cr),  free(im_2DES_R_cr);
  fclose(Rate);

  printf("the signal part can run successfully \n");

}









