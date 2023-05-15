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
    float *re_2DES_pa, *im_2DES_NR_pa, *im_2DES_R_pa;
    float *re_2DES_pe, *im_2DES_NR_pe, *im_2DES_R_pe;
    float *re_2DES_cr, *im_2DES_NR_cr, *im_2DES_R_cr; 
    float **re_2DES_pa_sum, **im_2DES_NR_pa_sum, **im_2DES_R_pa_sum;
    float **re_2DES_pe_sum, **im_2DES_NR_pe_sum, **im_2DES_R_pe_sum;
    float **re_2DES_cr_sum, **im_2DES_NR_cr_sum, **im_2DES_R_cr_sum; 

    int pro_dim;
    pro_dim=project_dim(non);
    re_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_doorway   = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    re_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_SE = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));   
    re_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_GB = (float *)calloc(non->tmax*9*pro_dim,sizeof(float)); 
    re_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));
    im_window_EA = (float *)calloc(non->tmax*9*pro_dim,sizeof(float));  
    re_2DES_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_NR_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_R_pa = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    re_2DES_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_NR_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_R_pe = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    re_2DES_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_NR_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    im_2DES_R_cr = (float *)calloc(3*non->tmax*non->tmax*non->tmax2,sizeof(float));
    re_2DES_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_NR_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_R_pa_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    re_2DES_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_NR_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_R_pe_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    re_2DES_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_NR_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));
    im_2DES_R_cr_sum = (float **)calloc2D(non->tmax3,non->tmax1,sizeof(float),sizeof(float*));

  
    //re_2DES = (float *)calloc(non->tmax*non->tmax*non->tmax2,sizeof(float));
    //im_2DES = (float *)calloc(non->tmax*non->tmax*non->tmax2,sizeof(float));

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
      if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
     (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      {
        CG_2DES_window_GB(non, re_window_GB, im_window_GB); 
      }
            if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
     (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      {
        CG_2DES_window_EA(non, re_window_EA, im_window_EA); 
      }
                    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
         (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      { 
      CG_full_2DES_segments(non, re_2DES_pa,im_2DES_NR_pa,im_2DES_R_pa,
                          re_2DES_pe, im_2DES_NR_pe, im_2DES_R_pe,
                          re_2DES_cr, im_2DES_NR_cr, im_2DES_R_cr);
      }
                    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
         (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))
      {
      combine_CG_2DES(non,re_2DES_pa_sum,im_2DES_NR_pa_sum,im_2DES_R_pa_sum,
                          re_2DES_pe_sum,im_2DES_NR_pe_sum,im_2DES_R_pe_sum,
                          re_2DES_cr_sum,im_2DES_NR_cr_sum,im_2DES_R_cr_sum);
      }
                    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
         (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES")))


    free(re_doorway);
    free(im_doorway);
    free(re_window_SE);
    free(im_window_SE);
    free(re_window_GB);
    free(im_window_GB);
    free(re_window_EA);
    free(im_window_EA);
    free(re_2DES_pa);
    free(im_2DES_NR_pa);
    free(im_2DES_R_pa);
    free(re_2DES_pe);
    free(im_2DES_NR_pe);
    free(im_2DES_R_pe);
    free(re_2DES_cr);
    free(im_2DES_NR_cr);
    free(im_2DES_R_cr);
    free(re_2DES_pa_sum) ;
    free(im_2DES_NR_pa_sum) ;
    free(im_2DES_R_pa_sum) ;
    free(re_2DES_pe_sum);
    free(im_2DES_NR_pe_sum) ;
    free(im_2DES_R_pe_sum);
    free(re_2DES_cr_sum);
    free(im_2DES_NR_cr_sum);
    free(im_2DES_R_cr_sum);    

    return;
}

//void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway);

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
}



void CG_2DES_P_DA(t_non *non,float *P_DA,float* K, float* P0, int N){
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
    
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                iP0[a + c * N] += ivecR[a + b * N] * P0[b + c * N];
            }
        }
    }    // iP0 = ivecR*P(0)

    // Loop over t2
    /*Here we assume the P0 is a N*N matrix*/
    for (int nt2 = 0; nt2<non->tmax2; nt2++) {
        cnr = (float *)calloc(N * N, sizeof(float));
        for (a = 0; a < N; a++) {
            for (b = 0; b < N; b++) {
                cnr[a + b * N] += exp(eigK_re[a] * nt2 * non->deltat) * iP0[a + b*N];
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
                fprintf(outone,"%f ",P_DA[t2+(N*a+b)*non->tmax2]);
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
          eq_den(Hamil_i_e,rho_l,N,non);
            //printf("rho_l "); 
            //printf("%f \n", rho_l[0]);

          /*printf("Hamil_i_e  "); 
          printf("%f \n", Hamil_i_e); */
          // Multiply the density operator to dipole operator,vecr, as it is only the real number.
          for (a=0;a<N;a++){
          /*printf("a "); 
          printf("%f \n", a); 
          printf("N "); 
          printf("%f \n", N); */
            for (b=0;b<N;b++){
              mid_ver[a]+=rho_l[a+b*N]*vecr[b];
                //printf(" SE vecr[b] "); 
                //printf("%f \n", vecr[a]); 
                //printf("rho_l "); 
                //printf("%f \n", rho_l[a+b*N]); 
                /*printf("a "); 
                printf("%f \n", a);
                printf("b "); 
                printf("%f \n", b);*/  
          /*printf("vecr[b] "); 
          printf("%f \n", vecr); 
          printf("mid_ver "); 
          printf("%f \n", mid_ver);*/                     
            }
          }
          // Update dipole operator
          for (a=0;a<N;a++){
            vecr[a]=mid_ver[a];
          }  

          //printf("%f \n", vecr); 
          for (beta=0;beta<3;beta++){
            /* Read mu(tj) */
            if (!strcmp(non->hamiltonian,"Coupling")){
              copyvec(mu_xyz+non->singles*beta,mu_eg,non->singles);
            } else {
              if (read_mue(non,mu_eg,mu_traj,tj,beta)!=1){
          //Here we generate the equilibrium density operator
            
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
}





void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA){
   /* Initialize variables*/
 /* The window part for SE*/
  float *Hamil_i_e,*Hamil_i_ee;
  float *mu_eg;
  float *vecr,*veci,*rho_i;   
  float *mu_xyz;

  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;
  /* Floats */
  float shift1;
  float *mid_rho,*mid_u,*up_rho_u;


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
  rho_l=(float *)calloc(N*N,sizeof(float));
// reserve memory for transition dipole from first to second
  float* fr = calloc(nn2*N, sizeof(float));
  float* fi = calloc(nn2*N, sizeof(float));
   rho_i = (float *)calloc(N,sizeof(float));
  //here the mid_vcr is just used  as the mid part for multiply the dipole with density operator
  float *mid_ver;
  mid_ver = (float *)calloc(N,sizeof(float));


//Here we reserve memory for matrix multiplication
  mid_rho=(float *)calloc(N,sizeof(float));
  mid_u=(float *)calloc(N,sizeof(float));
  up_rho_u = (float *)calloc(N*N,sizeof(float));
  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);
  shift1=(non->max1+non->min1)/2;
  printf("Frequency shift %f.\n",shift1);
  non->shifte=shift1;


// reserve memory for transition dipole
  vecr=(float *)calloc(N,sizeof(float));	
  veci=(float *)calloc(N,sizeof(float));
  mu_eg=(float *)calloc(N,sizeof(float));
  mu_xyz=(float *)calloc(N*3,sizeof(float));


  /* check projection */   
  if (non->Npsites!=non->singles){
      printf("Input segment number and the projection segment number are different.\n");
      printf("please check the input file or the projection file\n");
      exit(1);
  }
  /* projection */
  pro_dim=project_dim(non);

  Hamil_i_e=  (float *)calloc(nn2,sizeof(float));
  Hamil_i_ee= (float *)calloc(nn2,sizeof(float));

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
          eq_den(Hamil_i_e,rho_l,N,non);  
          clearvec(rho_i,non->singles);
          //printf ("rho_i");
          //printf("%f \n", rho_i[2]);        

          for (a=0;a<N;a++){
            //for (b=0;b<N;b++){
              //index = a+b*N
              //mid_rho[b]+=rho_l[a+b*N];
            //}
            //dipole_double_CG2DES(non, vecr, mid_rho, rho_i, fr[a], fi[a]);
            dipole_double_CG2DES(non, vecr, rho_l+a*N, rho_i, fr+a*nn2, fi+a*nn2); 
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
            //dipole_double_CG2DES(non, mu_eg, vecr, veci, fr, fi);            
            /* Here calculate doorway function/
            /* Inner product for all sites*/
            /*here we need to make clear the seg_num and the site_num
            the number of the position should be decided by the segment number,
            we should sum all the value in one segment*/
            
            /* Multiply with mu_ef */
            /* Take fr and fi back to vecr and veci */
            for (a=0;a<N;a++){
              //The dimention of vecr is only   vecr=(float *)calloc(N,sizeof(float));	
                dipole_double_inverse_CG2DES(non, mu_eg, fr+a*nn2, fi+a*nn2, vecr+a*N, veci+a*N);
                //dipole_double_inverse_CG2DES(non, mu_eg, fr+a*nn2, fi+a*nn2, vecr, veci);
            }
            /* Propagate back to time ti； tj=ti+t1; */
            for (tk=tj;tk>ti;tk--){
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
              if (non->Npsites==non->singles){
                zero_coupling(Hamil_i_ee,non);
              } else {
                printf("Segment number and the projection number are different");
                exit(1);
              }              
              /* Propagate vecr and veci backwards with propagate */
              /* Propagate dipole moment */
            for (a=0;a<N;a++){
               propagate_vector(non,Hamil_i_ee,vecr+a*N,veci+a*N,1,samples,tk*alpha);
            }
            }

            /* Calculate window function by taking trace of vecr and veci */
            for (site_num=0;site_num<non->singles;site_num++){
              seg_num=non->psites[site_num];
              /*this equation is make sure the calculated data in the right position */
              index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
              re_window_EA[index]+=vecr[site_num*N+site_num];
              im_window_EA[index]+=veci[site_num*N+site_num]; 
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
          for (a=0;a<N;a++){
            propagate_vec_coupling_S_doubles_ES(non, Hamil_i_e, fr+a*nn2, fi+a*nn2, non->ts);    
                  
          }
     
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

 

  free(mu_xyz);
  free(Hamil_i_e);
  free(Hamil_i_ee);
  free(rho_l);
  free(mid_ver);
  free(fr);
  free(fi);
  free(vecr);
  free(veci);
  free(mu_eg);
    //printf("this is a test \n"); 
  //exit(1);
}
//void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB)
//void CG_2DES_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
//void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_2DES_pa,float *im_2DES_NR_pa,float *im_2DES_R_pa,
                          float *re_2DES_pe,float *im_2DES_NR_pe,float *im_2DES_R_pe,
                          float *re_2DES_cr,float *im_2DES_NR_cr,float *im_2DES_R_cr){

                          }
void combine_CG_2DES(t_non *non,float **re_2DES_pa_sum,float **im_2DES_NR_pa_sum,float **im_2DES_R_pa_sum,
                          float **re_2DES_pe_sum,float **im_2DES_NR_pe_sum,float **im_2DES_R_pe_sum,
                          float **re_2DES_cr_sum,float **im_2DES_NR_cr_sum,float **im_2DES_R_cr_sum){

                          }








