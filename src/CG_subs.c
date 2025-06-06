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

/* Frequently used indexing function */
/* It is needed for the response functions, which depends both on segment number,
   two cartesian directions (alpha, beta) and a time (t1) limited by the coherence time.*/
   int CG_index(t_non *non,int seg_num,int alpha,int beta,int t1){
    int index;
    index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
    return index;
}

/* Normalize response function with respect to the number of samples */
void normalize_DW(t_non *non,float *re,float *im,int samples){
  int MAX_NUMBERS,number;
  int pro_dim;

  pro_dim=project_dim(non);
  MAX_NUMBERS = pro_dim*9*non->tmax;
  for (number=0;number<MAX_NUMBERS;number++){
    im[number]=im[number]/samples;
    re[number]=re[number]/samples;
  }
}

/* This function writed the final response (doorway/window) functions to file. The function below
   reads these files. */
void write_response_to_file(t_non *non,char fname[],float *im,float *re,int tmax){
  FILE *file_handle;
  int MAX_NUMBERS,number;
  int pro_dim;
  file_handle=fopen(fname,"w");
  pro_dim=project_dim(non);
  MAX_NUMBERS = pro_dim*9*non->tmax;
  for (number=0;number<MAX_NUMBERS;number++){
    fprintf(file_handle,"%e %e ",im[number],re[number]);
  }
  fclose(file_handle);
}

/* Read the absorption/emission function from file. It follows the format
   from writing in the function above. */
void read_doorway_window_from_file(t_non *non,char fname[],float *im,float *re,int tmax){
  FILE *file_handle;
  int pro_dim,i;
  int num_count = 0;
  float number;
  int MAX_NUMBERS;
  pro_dim=project_dim(non);
  MAX_NUMBERS = pro_dim*9*non->tmax;
  file_handle=fopen(fname,"r");
  if (file_handle == NULL) {
    printf("Error opening the file %s.\n",fname);
    exit(0);
  }
  
  // Read numbers line by line
  for (i=0;i<MAX_NUMBERS;i++){
    fscanf(file_handle,"%f %f",&im[i],&re[i]);
  }
    fclose(file_handle);
}
  
/* Calculate the doorway functions */
void CG_doorway(t_non *non,float *re_doorway,float *im_doorway){
  /* Initialize variables*/
  float *Hamil_i_e;
  float *mu_eg;
  float *vecr,*veci;
  float *mu_xyz;
  
  /* File handles */
  FILE *H_traj,*mu_traj;
  FILE *C_traj;
  FILE *outone,*log;
  FILE *Cfile;
  FILE *WTime;
  /* Floats */
  float shift1;
  
  /* Integers */
  int nn2;
  int itime,N_samples;
  int samples;
  int alpha,beta,a_b,ti,tj,i,j; 
  /* alpha and beta are corresponding the dipole for different times */
  int t1,t2;
  int elements;
  int cl,Ncl;
  int pro_dim,ip;
  int a,index,seg,site_num,seg_num;
  /* site_num is the number of the site */ 
  /* seg_num is used to number segments */
  
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
  fprintf(log,"Calculating doorway functions.\n");
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
  
      /* Configuration belong to cluster */
      if (non->cluster==cl){
        Ncl++;
      }
    }
  
    /* Include snapshot if it is in the cluster or if no clusters are defined */
    if (non->cluster==-1 || non->cluster==cl){   
      for (alpha=0;alpha<3;alpha++){
        /* Read mu(ti) */
          read_dipole(non,mu_traj,vecr,mu_xyz,alpha,ti);         
          clearvec(veci,non->singles);
           
          /* Loop over coherence time */
        for (t1=0;t1<non->tmax;t1++){
          tj=ti+t1;
          /* Read Hamiltonian */
            read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
  
          for (beta=0;beta<3;beta++){
              /* Read mu(tj) */
              read_dipole(non,mu_traj,mu_eg,mu_xyz,beta,tj);
  
              /* Here calculate doorway function */
              /* Inner product for all sites */
              /* Here we define the seg_num and the site_num the number
                 * depending on the segment number,
                 * we sum all contributions from one segment*/
  
              for (site_num=0;site_num<non->singles;site_num++){
                  seg_num=non->psites[site_num];
                  /* This equation is make sure the calculated data in the right
                     * position */
                   //index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                      index=CG_index(non,seg_num,alpha,beta,t1);
                  re_doorway[index]+=mu_eg[site_num]*vecr[site_num];
                  im_doorway[index]+=mu_eg[site_num]*veci[site_num]; 
              }
          }  
  
          /* Do projection and make sure the segment number equal to the
           * segment number in the projection file. */   
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
    fprintf(log,"Finished sample (doorway function) %d\n",samples);        
    time_now=log_time(time_now,log);
    fclose(log);
  }
  
  /* Normalize doorway/window function with respect to number of samples */
  normalize_DW(non,im_doorway,re_doorway,samples);
  /* Save  the imaginary part for time domain response */
  write_response_to_file(non,"doorway.dat",im_doorway,re_doorway,non->tmax1);
  /* The calculation is finished we can close all auxillary arrays before
   * writing output to file. */
  free(vecr);
  free(veci);
  free(mu_eg);
  free(mu_xyz);
  free(Hamil_i_e);
  
  /* Print information on number of realizations included belonging to the
   * selected cluster and close the cluster file. (Only to be done if cluster
   * option is active.) */
  if (non->cluster!=-1){
    printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
    fclose(Cfile);
  }
  /* Close Trajectory Files */
  fclose(mu_traj),fclose(H_traj);
  
  printf("The doorway part successfully completed!\n");  
}
  
/* Determine the transfer propability matrix */
void CG_P_DA(t_non *non,float *P_DA,int N){
  float *eigK_re, *eigK_im; // eigenvalues of K
  float *evecL, *evecR; // eigenvectors of K
  float *ivecR, *ivecL; //inverse eigenvectors of K
  float *cnr;
  float sum_eig_im;
  int a, b, c,i,j;
  float *K;
  int nt2;
  float _Complex *ivecL_com, *ivecR_com;
  float _Complex *ivecL_com_inv, *ivecR_com_inv;
  float _Complex *cnr_com;
  float _Complex *P_DA_com;
  FILE *outone;
  FILE *Rate;
  FILE *WTime;
  float factor;
  nt2 =non->tmax2;
  factor =  ( nt2 * non->deltat)/1000; /* The rate matrix is in ps-1 */
  sum_eig_im=0;
  printf("Calculating population transfer matrix for t2 %f fs!\n",factor*1000);
  eigK_re = (float *)calloc(N,sizeof(float));
  eigK_im = (float *)calloc(N,sizeof(float));
  evecL = (float *)calloc(N*N,sizeof(float));
  evecR = (float *)calloc(N*N,sizeof(float));
  ivecL = (float *)calloc(N*N,sizeof(float));
  ivecR = (float *)calloc(N*N,sizeof(float));
  cnr = (float *)calloc(N * N, sizeof(float));
  K=(float *)calloc(N*N,sizeof(float));
  cnr_com =(float _Complex *)calloc(N * N, sizeof(float _Complex));
  ivecL_com = (float _Complex *)calloc(N * N, sizeof(float _Complex));
  ivecR_com = (float _Complex *)calloc(N * N, sizeof(float _Complex));
  ivecL_com_inv = (float _Complex *)calloc(N * N, sizeof(float _Complex));
  ivecR_com_inv = (float _Complex *)calloc(N * N, sizeof(float _Complex));
  P_DA_com = (float _Complex *)calloc(N * N, sizeof(float _Complex));

  clearvec(P_DA,N*N);
  /* Open the rate matrix file */
  Rate=fopen("RateMatrix.dat","r");
  if (Rate==NULL){
    printf("RateMatrix file not found!\n");
    exit(1);
  }
  /* Read rate matrix */
  printf("\nUsing the Rate matrix:\n");
  for (a=0;a<N*N;a++){
    if (fscanf(Rate,"%f",&K[a])!=1){
      printf("Error in reading in rate matrix!\n");
      exit(0);
    }
    printf("%f ",K[a]);
  }
  printf("\n\nCompleted reading the rate matrix.\n");
  
  /* Diagonalize K matrix */
  diagonalize_real_nonsym(K, eigK_re, eigK_im, evecL, evecR, ivecL, ivecR, N);
  /* Check if the eigenvalues contain imaginary parts */
  for (int a = 0; a<N; a++) {
      sum_eig_im += fabs(eigK_im[a]);
  }
  
  /* Calculate the transfer matrix */
  if (sum_eig_im == 0.0){ 
    /* Calculate the inverse part */
    inversie_real_matrix(eigK_re, eigK_im, evecL, evecR, ivecL, ivecR, N);
    /* Calculate P(t2) = expm(-K*t2)*P(0) = evecR*exp(Eig*t2)*ivecR*P(0) */
    clearvec(cnr,N*N); /* Empty auxillary vector */ 
    for (a = 0; a < N; a++) {
      for (b = 0; b < N; b++) {
        cnr[a + b * N] += exp(eigK_re[a]/1000 * nt2 * non->deltat) * ivecL[a + b*N];
      }
    }   // exp(Eig*t2)*iP0
  
    for (a = 0; a < N; a++) {
      for (b = 0; b < N; b++) {
        for (c = 0; c < N; c++) {
          P_DA[c+N*a] += evecL[a + b * N] * cnr[b + c * N];
        }
      }
    }   // evecR*cnr
  } else {
    /* Warn user if there are imaginary components */
    printf(YELLOW "The rate matrix has imaginary parts!\n" RESET);
    printf("Sum of imaginary contributions: %f\n", sum_eig_im);
  
    /* This is for complex eigenvalues and vectors */
    inversie_complex_matrix(eigK_re, eigK_im, evecL, evecR, ivecL_com_inv, ivecR_com_inv, N);
    /* Here we reconstruct the matrix when complex */
    /* Actually we only need the left one */
    for (int i = 0; i < N; i++){
      // If the current and next eigenvalues form a complex conjugate pair
      if (i < N - 1 && eigK_im[i] != 0.0 && eigK_im[i + 1] == -eigK_im[i]) {
        for (int j = 0; j < N; j++) {
          // u(j) = VL(:,j) + i*VL(:,j+1)
          ivecL_com[i * N + j] = evecL[i * N + j]-evecL[(1 + i) * N + j]*_Complex_I;
          ivecL_com[(1+i) * N + j] = evecL[i * N + j]+evecL[(1 + i) * N + j]*_Complex_I;         
        }
        // Skip the next eigenvector (since it's part of the complex conjugate pair)
        i++;
      } else {
        for (int j = 0; j < N; j++) {
          ivecL_com[i * N + j] = evecL[i * N + j]+0*_Complex_I;
        }
      }
    }
  
    /* Calculate P(t2) = expm(-K*t2)*P(0) = evecR*exp(Eig*t2)*ivecR*P(0) */
    clearvec(cnr,N*N); /* Empty auxillary vector */ 
    for (a = 0; a < N; a++) {
      for (b = 0; b < N; b++) {
        cnr_com[a + b * N] += exp((eigK_re[a])*factor) *((creal(ivecL_com_inv[a + b*N])*cos((eigK_im[a])*factor)-sin((eigK_im[a])*factor)*cimag(ivecL_com_inv[a + b*N]))+(cimag(ivecL_com_inv[a + b*N])*cos((eigK_im[a])*factor)+sin((eigK_im[a])*factor)*creal(ivecL_com_inv[a + b*N]))*_Complex_I);     
      }
    }  
    for (a = 0; a < N; a++) {
      for (b = 0; b < N; b++) {
        for (c = 0; c < N; c++) { 
          P_DA_com[c+a*N] += (creal(ivecL_com[a + b * N])*creal(cnr_com[b + c * N])-cimag(ivecL_com[a + b * N])*cimag(cnr_com[b + c * N]))+(creal(ivecL_com[a + b * N])*cimag(cnr_com[b + c * N])+creal(cnr_com[b + c * N])*cimag(ivecL_com[a + b * N])) *_Complex_I;
        }
      }
    }   // evecL*cnr
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        if (fabs(cimag(P_DA_com[i*N+j])) < 0.00001) {
          P_DA[i*N+j]=creal(P_DA_com[i*N+j]);
        } else {
          printf("\nThe imaginary part is none zero\n");
        }      
      }
    }
  }
  
  /* Write to file */
  outone=fopen("KPop.dat","w");
  fprintf(outone,"%f ",nt2*non->deltat);
  for (int a=0;a<N;a++){
      for (int b=0;b<N;b++){
          fprintf(outone,"%f ",P_DA[a*N+b]);                               
      }
  }
  fprintf(outone,"\n"); 
  fclose(outone);
  
  free(eigK_im);
  free(eigK_re);
  free(evecL);
  free(evecR);
  free(ivecL);
  free(ivecR);
  free(P_DA_com);
  free(ivecL_com);
  free(ivecR_com);
  free(ivecL_com_inv);
  free(ivecR_com_inv);
  free(cnr_com);
  
  printf("The waiting time propagation successfully completed!\n");  
  return;
}
  
  /* Calcualte doorway function for stimulated emission */
  void CG_window_SE(t_non *non, float *re_window_SE, float *im_window_SE){
    /* Initialize variables*/
    /* The window part for SE*/
    float *Hamil_i_e;
    float *mu_eg;
    float *vecr,*veci;
    float *mu_xyz;
    float shift1;
    float *rho_l;
    float *mid_ver;
  
    /* File handles */
    FILE *H_traj,*mu_traj;
    FILE *C_traj;
    FILE *outone,*log;
    FILE *Cfile;
  
    /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int alpha,beta,a_b,ti,tj,i,j; /* alpha and beta are corresponding the dipole for different segment,  a_b  = alpha*beta, which used in the writing file */
    int t1,t2;
    int elements;
    int cl,Ncl;
    int pro_dim,ip;
    int a,b,index,seg,site_num,seg_num; /* site num is the number of the site, which used to make sure the index number later */
    int N;
  
    N = non->singles;
    // reserve memory for the density operator
    rho_l=(float *)calloc(N*N,sizeof(float));
    // here the mid_vcr is just used  as the mid part for multiply the dipole with density operator
    mid_ver = (float *)calloc(N,sizeof(float));
    // reserve memory for transition dipole
    vecr=(float *)calloc(non->singles,sizeof(float));	
    veci=(float *)calloc(non->singles,sizeof(float));
    mu_eg=(float *)calloc(non->singles,sizeof(float));
    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
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
    fprintf(log,"Calculating SE window functions.\n");
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
            read_dipole(non,mu_traj,vecr,mu_xyz,alpha,ti);
          clearvec(veci,non->singles);
  
            /* Read Hamiltonian */
          read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
              
            eq_den(Hamil_i_e,rho_l,N,non);
          //  write_ham_to_file("Density.dat",rho_l,N);
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
              read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
  
            for (beta=0;beta<3;beta++){
              /* Read mu(tj) */
              read_dipole(non,mu_traj,mu_eg,mu_xyz,beta,tj);
                                     
              /* Here calculate window function
                 Inner product for all sites
                 here we need to make clear the seg_num and the site_num the
                     number of the position should be decided by the segment
                     number, we should sum all the value in one segment */
  
              for (site_num=0;site_num<non->singles;site_num++){
                seg_num=non->psites[site_num];
   
                /* this equation is make sure the calculated data in
                       the right position */
                     //index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                    index=CG_index(non,seg_num,alpha,beta,t1);
                re_window_SE[index]+=mu_eg[site_num]*vecr[site_num];
                im_window_SE[index]+=mu_eg[site_num]*veci[site_num];          
              }
            }  
  
            /* Zero the coupling between different segments.*/   
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
      fprintf(log,"Finished sample (window SE) %d\n",samples);        
      time_now=log_time(time_now,log);
      fclose(log);
    }

    /* Normalize doorway/window function with respect to number of samples */
    normalize_DW(non,im_window_SE,re_window_SE,samples);
    /* Save the time domain response */
    write_response_to_file(non,"windows_SE.dat",im_window_SE,re_window_SE,non->tmax1); 
    /* Close all auxillary arrays before writing */
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
    fprintf(log,"Finished Calculating SE window functions!\n");
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
    printf("The window function SE part successfully completed!\n");  
  }
  
  /* Calculate the window function for ground state bleach */
  void CG_window_GB(t_non *non,float *re_window_GB,float *im_window_GB){
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
    float *rho_l;
  
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
    
    int N;
    N = non->singles;
    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;
    rho_l=(float *)calloc(N*N,sizeof(float));
    // Reserve memory for transition dipole
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
    fprintf(log,"Calculating GB window functions.\n");
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
            read_dipole(non,mu_traj,vecr,mu_xyz,alpha,ti);
          clearvec(veci,non->singles);
          /* Loop over coherence time */
          for (t1=0;t1<non->tmax;t1++){
            tj=ti+t1;
            /* Read Hamiltonian */
              read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
            for (beta=0;beta<3;beta++){
              /* Read mu(tj) */
                read_dipole(non,mu_traj,mu_eg,mu_xyz,beta,tj);
              
              /* Here calculate window function */
              /* Inner product for all sites */
              /* here we need to make clear the seg_num and the site_num
                     the number of the position should be decided by the segment
                     number, we should sum all the value in one segment */
  
              for (site_num=0;site_num<non->singles;site_num++){
                seg_num=non->psites[site_num];
                /* this equation is make sure the calculated data in the right position */
                //index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                    index=CG_index(non,seg_num,alpha,beta,t1);
                re_window_GB[index]+=mu_eg[site_num]*vecr[site_num];
                im_window_GB[index]+=mu_eg[site_num]*veci[site_num]; 
              }
            }   
            /* Zero coupling between different coulpmings.*/   
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
      fprintf(log,"Finished sample (window GB) %d\n",samples);        
      time_now=log_time(time_now,log);
      fclose(log);
    }
   
    /* Normalize doorway/window function with respect to number of samples */
    normalize_DW(non,im_window_GB,re_window_GB,samples);
    /* Save the time domain response */
    write_response_to_file(non,"windows_GB.dat",im_window_GB,re_window_GB,non->tmax1); 
    /* Close all auxillary arrays before writing */
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
    printf("The GB window function part successfully completed!\n");  
  }
  
  /* Calculate the window function for the excited state absorption */
  void CG_window_EA(t_non *non,float *re_window_EA,float *im_window_EA){
    /* Initialize variables*/
    /* The window part for SE*/
    float *Hamil_i_e,*Hamil_i_ee;
    float *mu_eg;
    float *vecr, *vecr1, *veci, *rho_i;   
    float *mu_xyz;
    float *rho_l;
  
    /* File handles */
    FILE *H_traj,*mu_traj;
    FILE *C_traj;
    FILE *outone,*log;
    FILE *Cfile;
    /* Floats */
    float shift1;
    float *fr;
    float *fi;
    /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int alpha,beta,a_b,ti,tj,tk,i,j; /*alpha and beta are corresponding the dipole for different segment,  a_b  = alpha*beta, which used in the writing file*/
    int t1,t2;
    int elements;
    int cl,Ncl;
    int pro_dim,ip;
    int N;
    int a,b,index,seg,site_num,seg_num;/*site num is the number of the site, which used to make sure the index number later*/
    
    N = non->singles;
     /* Allocate memory*/  /*(tmax+1)*9*/
    nn2=non->singles*(non->singles+1)/2;
    /* projection */
    pro_dim=project_dim(non);
    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;
  
    // Reserve memory for transition dipole
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
    fprintf(log,"Calculating EA window functions.\n");
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
            read_dipole(non,mu_traj,vecr1,mu_xyz,alpha,ti);
          clearvec(rho_i,non->singles*non->singles); 
  
          // Here we generate the equilibrium density operator
            /* Read Hamiltonian */
            read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
          eq_den(Hamil_i_e,rho_l,N,non);       
          
          for (a=0;a<N;a++){
            dipole_double_CG2DES(non, vecr1, rho_l+a*N, rho_i+a*N, fr+a*nn2, fi+a*nn2); 
          }
          /* Loop over coherence time */
          for (t1=0;t1<non->tmax;t1++){
            tj=ti+t1;
            
            for (beta=0;beta<3;beta++){
              /* Read mu(tj) */
                  read_dipole(non,mu_traj,mu_eg,mu_xyz,beta,tj);
              
                  /* Here calculate window function */
              /* Inner product for all sites */
              /* here we need to make clear the seg_num and the site_num
                     the number of the position should be decided by the segment
                     number, we should sum all the value in one segment */
              /* Multiply with mu_ef Take fr and fi back to vecr and veci */
              clearvec(vecr,non->singles*non->singles);
              clearvec(veci,non->singles*non->singles);
              for (a=0;a<N;a++){
                // The dimention of vecr is only   vecr=(float *)calloc(N,sizeof(float));
                dipole_double_inverse_CG2DES(non, mu_eg, fr+a*nn2, fi+a*nn2, vecr+a*N, veci+a*N);
              }
              /* Propagate back to time ti； tj=ti+t1; */
              for (tk=tj;tk>ti;tk--){
                /* Read hamiltonian at tk */
                  read_Hamiltonian(non,Hamil_i_ee,H_traj,tk);
  
                /* Zero coupling in Hamiltonian (TLC check this!) */
                zero_coupling(Hamil_i_ee,non);
                /* Propagate single excited states vecr and veci backwards with propagate */
                /* Propagate dipole moment */
  #pragma omp parallel for
                for (a=0;a<N;a++){
                  propagate_vector(non,Hamil_i_ee,vecr+a*N,veci+a*N,-1,samples,tk*alpha);
                }
  //		    propagate_matrix(non,Hamil_i_ee,vecr,veci,-1,samples,tk*alpha);
              }
              /* Calculate window function by taking trace of vecr and veci */
              for (site_num=0;site_num<non->singles;site_num++){
                seg_num=non->psites[site_num];
                /* this equation is make sure the calculated data in the right position */
                //index=seg_num*9*non->tmax+alpha*3*non->tmax+beta*non->tmax+t1;
                    index=CG_index(non,seg_num,alpha,beta,t1);
                re_window_EA[index]+=vecr[site_num*N+site_num];
                im_window_EA[index]+=veci[site_num*N+site_num]; 
              }
            }  
            /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
              /* Read Hamiltonian */
              read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
              
              /* Zero coupling between different segments */
            zero_coupling(Hamil_i_e,non);
            /* Propagate dipole moment */
              /* Currently only coupling propagation done! TLC */
  #pragma omp parallel for
            for (a=0;a<N;a++){
              propagate_vec_coupling_S_doubles_ES(non, Hamil_i_e, fr+a*nn2, fi+a*nn2, non->ts);   
            }
          }
        }
      }
      /* Update Log file with time and sample number */
      log=fopen("NISE.log","a");
      fprintf(log,"Finished sample (window EA) %d\n",samples);        
      time_now=log_time(time_now,log);
      fclose(log);
    }

    /* Normalize doorway/window function with respect to number of samples */
    normalize_DW(non,im_window_EA,re_window_EA,samples);
    /* Save  the imaginary part for time domain response */
    write_response_to_file(non,"windows_EA.dat",im_window_EA,re_window_EA,non->tmax1);
    /* Free all auxillary arrays */
    free(mu_xyz);
    free(Hamil_i_e);
    free(Hamil_i_ee);
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
    printf("The EA window function part successfully completed!\n");
    return;  
  }
  