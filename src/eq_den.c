#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "propagate.h"
#include "absorption.h"
//#include "luminescence.h"
#include "1DFFT.h"
#include "project.h"
#include "eq_den.h"


/*Here we defined the equilibrium density operator*/
void eq_den(float *Hamiltonian_i, float* rho_l, int N, t_non *non){
  int index,N;
  float *H,*e,*e_1,;
  float *rho_r;
  //float *rho_l;
  float *diag_sum;  /*Here is store the diagonal element for one segment*/
  N=non->singles;
  H=(float *)calloc(N*N,sizeof(float));
  e=(float *)calloc(N,sizeof(float));
  e_1=(float *)calloc(N,sizeof(float));
  rho_r=(float *)calloc(N*N,sizeof(float));
  //rho_l=(float *)calloc(N*N,sizeof(float));
  int a,b,c,pro_dim;
  int site_num,site_num_1,site_num_2,seg_num;/*site num is the number of the site, which used to make sure the index number later*/

  float kBT=non->temperature*0.6950859; // Here 0.6950859 is the K_b in the unit of cm-1
  /*Here I have one question, the input1D file does not contain temperature, 
  where we should include it */
  
  float u,i_u; 
  /* projection */
  pro_dim=project_dim(non);
  diag_sum = (float *)calloc(pro_dim,sizeof(float));

  /* Do projection and make sure the segment number equal to the segment number in the projection file.*/   
  if (non->Npsites==non->singles){
    zero_coupling(Hamil_i_e,non);
  } else {
    printf("Segment number and the projection number are different");
    exit(1);
  }

  // Build Hamiltonian, convert the triangular Hamiltonian to the  square Hamiltonian
  for (a=0;a<N;a++){
    /* Fill diagonal elements */
    H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; 
    /* Fill couplings */
    for (b=a+1;b<N;b++){
      H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
      H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
    }
  }
/*Here we diagonalize the square Hamiltonian, the H is replaced by the eigenvector in eigen basis, e is the eigenvalue */

  diagonalizeLPD(H,e,N);
 
  // Exponentiate [u=exp(-H*kBT)]
  for (a=0;a<N;a++){
    e_1[a]=exp(-e[a]*kBT);
    u=u+e_1[a];
  }
  /*Here calculate the inverse of u*/
  i_u=1.0/u;

  // Transform to site basis, H*u*H
  /*Here we first calculate the right side u*H, which output rho_r */
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      rho_r[b+a*N]+=H[b+a*N]*e_1[b]*i_u;
    }
  }  
  /*Secondly, we calculate the left side H*u_r, which output rho_l */
  for (a=0;a<N;a++){
    for (b=0;b<N;b++){
      for (c=0;c<N;c++){
        rho_l[a+c*N]+=H[b+a*N]*rho_r[b+c*N];
      }
    }
  }
/* The weights in the site basis were calculated, 
Here we should not combine the two loop in one loop as it would results in the heary calcuation.*/

/*Here we need to normalize rho, exp(-KB*H)/[Trexp(-KB*H)] within every segment*/
/*We first sum the diagonal element in one segment*/
  for (site_num=0;site_num<non->singles;site_num++){
    seg_num=non->psites[site_num];
    diag_sum[seg_num]+= rho_l[site_num][site_num];
  }

    /*Secondly, we normalize the density matrix within one segment*/
    //for (seg_num=0;seg_num<pro_dim;seg_num++){
  for (site_num_1=0;site_num_1<non->singles;site_num_1++){
    seg_num=non->psites[site_num_1];
    for (site_num_2=0;site_num_2<non->singles;site_num_2++){
      if (seg_num==non->psites[site_num_2]){
        rho_l[seg_num_1][seg_num_2]=rho_l[seg_num_2][seg_num_1]=rho_l[seg_num_1][seg_num_2]/diag_sum[seg_num];
      }
    }
  }
  

  free(H);
  free(e_1);
  free(e);
  free(rho_r);
  //free(rho_l);
  return ;
}


