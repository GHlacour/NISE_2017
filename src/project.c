#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "project.h"
#include "randomlib.h"
#include "util/asprintf.h"

// This subroutine nullify all elements of vector for non selected sites
void projection(float* phi, t_non* non) {
    int i;
    for (i = 0; i < non->singles; i++) {
        phi[i] = phi[i] * non->psites[i];
    }
    return;
}

// This subroutine nullify all elements of vector for non selected sites in
// multi projection
void multi_projection(float *phi_in,float *phi_out,t_non *non,int ip){
    int i;
    for (i=0;i<non->singles;i++){
    if (non->psites[i]==ip){
       phi_out[i]=phi_in[i];
    } else {
           phi_out[i]=0.0f;
    }
    }
}

/* This subroutine will nullify all the excitonic couplings between different segments */
void zero_coupling(float *Hamil_i_e, t_non *non){
   /* Defining a double loop because we are checking the couplings between the pairs of molecules */
    int i,j;
    for (i=0;i<non->singles;i++){
        for (j=i+1;j<non->singles;j++){
            if (non->psites[i]!=non->psites[j]){
               Hamil_i_e[Sindex(i,j,non->singles)]==0;
            }
        }
    }
}

// Analyse projection input
int project_dim(t_non* non){
   int N,i;
   int max;
   max=0;
   if (non->Npsites==0){
      return 1;
   }
   if (non->Npsites<non->singles){
      return 1;
   }
   if (non->Npsites==non->singles){
      N=non->singles;
      for (i=0;i<non->singles;i++){
      if (non->psites[i]>=N){
             printf(RED "More projection segments than number of sites not allowed\n" RESET);
         exit(0);
      }
      if (non->psites[i]>max){
             max=non->psites[i];
      }
      }
      //printf("Identified %d projection segments\n",max+1);
      return max+1;
   }
printf(RED "Something went wrong in projection input analysis.\n");
   printf("You should not be able to end up here! Contact developers!\n" RESET);
   exit(0);
}

