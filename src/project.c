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

/*This subroutine will nullify all the excitonic couplings between different segments*/
void multi_projection_Hamiltonian(float *Hamil_i_e, t_non *non){
   /*Double loop: checking the couplings between the pairs of molecules */
    int i;
    int j;
    for (i=0;i<non->singles;i++){
        for (j=i+1;j<non->singles;j++){
            if (non->psites[i] != non->psites[j]){
               Hamil_i_e[Sindex(i,j,non->singles)]=0.0;
            } 
    	  }
    }
    return;
}

/*This subroutine will nullify all the excitonic couplings in the same segment */
void multi_projection_Coupling(float *Hamil_i_e, t_non *non){
   /*Double loop: checking the couplings between the pairs of molecules */
    int i;
    int j;
    for (i=0;i<non->singles;i++){
        for (j=i+1;j<non->singles;j++){
            if (non->psites[i] == non->psites[j]){
               Hamil_i_e[Sindex(i,j,non->singles)]=0.0;
            } 
    	  }
    }
    return;
}

/* This subroutine will nullify all the excitonic couplings between different segments */
void zero_coupling(float *Hamil_i_e, t_non *non){
   /* Defining a double loop because we are checking the couplings between the pairs of molecules */
    int i,j;
    for (i=0;i<non->singles;i++){
        for (j=i+1;j<non->singles;j++){
            if (non->psites[i]!=non->psites[j]){
               Hamil_i_e[Sindex(i,j,non->singles)]=0;
            }
        }
    }
   return;
}

/* Analyse projection input */
int project_dim(t_non* non){
   int N,i;
   int max;
   max=0;
   /* No projection specified all is one segment */
   if (non->Npsites==0){
      return 1;
   }
   /* Exactly one segment with Npsites is specified */
   if (non->Npsites<non->singles){
      return 1;
   }
   /* All sites are distributed in segments */
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
   if (non->Npsites>non->singles){
    printf(RED "More segments than sites was specified.\n");
    printf("Please, specify as many segments as sites or less!\n" RESET);
    exit(0);
   }
   printf(RED "Something went wrong in projection input analysis.\n");
   printf("You should not be able to end up here! Contact developers!\n" RESET);
   exit(0);
}

/* Find degeneracies of segments */
void project_degeneracies(t_non* non,int *degen,int segments){
    int i;
    /* Clear array */
    for (i=0;i<segments;i++){
        degen[i]=0;
    }
    /* Count members of each segment */
    for (i=0;i<non->singles;i++){
	    degen[non->psites[i]]++;
    }
    return;    
}

/* Find the indices in site basis of all sites belonging to the indicated segment */
int find_H_indices_segment(int *psites, int *H_indices_si,int si, t_non *non){
    // find the indices in the full system hamiltonian
    // for the indicated segment si such that this only needs to be done once
    int i, N, n_i;
    n_i = 0;
    N=non->singles;
    clearvec_int(H_indices_si,N);

    // overwrite the input H_index array, for the segment concerned
    for (i=0;i<N;i++){
       if (psites[i]==si){
           H_indices_si[n_i]=i;
	    n_i++;
        }
    }

    // return the number of molecules in this segment for convenient looping
    return n_i;
}

/* Isolate the segment hamiltonian in upper triangle format */
void isolate_segment_Hamiltonian_triu(float *Hamiltonian_full_triu, float *Hamiltonian_segment_triu, int *H_indices_si, int N_i, t_non *non){
    int N, H_a, H_b, i, j;
    int  H_triu_full_idx, H_triu_full_idx_inter, H_triu_si_idx;
    H_triu_si_idx = 0;
    H_triu_full_idx = 0;
    N = non->singles;
    clearvec(Hamiltonian_segment_triu,(N_i+1)*N_i/2);
 
    /* make use of the formula k = N(N-1)/2 - (N-i)(N-1-i)/2+j to convert to upper triangle index*/
    for (i=0;i<N_i;i++){
        H_a = H_indices_si[i];
	    H_triu_full_idx_inter = N*(N-1)/2 - (N-H_a)*(N-1-H_a)/2; 
        for (j=i;j<N_i;j++){
            H_b = H_indices_si[j];
	        H_triu_full_idx = H_triu_full_idx_inter + H_b;
	        Hamiltonian_segment_triu[H_triu_si_idx] = Hamiltonian_full_triu[H_triu_full_idx];
	        // printf("si idx %d\n",H_triu_si_idx);
	        H_triu_si_idx++;
	        //Hamiltonian_segment[i*N_i+j] = Hamiltonian_full[H_a*N+H_b];
        }
    }
    return;
}

/* Isolate the segment hamiltonian in full matrix format */
void isolate_segment_Hamiltonian(float *Hamiltonian_full, float *Hamiltonian_segment, int *H_indices_si, int N_i, t_non *non){
    int N, H_a, H_b, i, j;
    N = non->singles;
    clearvec(Hamiltonian_segment,N_i*N_i);
 
    for (i=0;i<N_i;i++){
        H_a = H_indices_si[i];
        for (j=0;j<N_i;j++){
            H_b = H_indices_si[j];
            Hamiltonian_segment[i*N_i+j] = Hamiltonian_full[H_a*N+H_b];
        }
    }
    return;
}

/* Isolate the inter-segment couplling Hamiltonian in full matrix format */
void isolate_coupling_block(float *J_full, float *J_ij, int N_i, int N_j, int *H_indices_si, int *H_indices_sj, t_non *non){
    int N, H_a, H_b, i, j;
    N = non->singles;
    clearvec(J_ij, N_i*N_j);

    for (i=0;i<N_i;i++){
        H_a = H_indices_si[i];
        for (j=0;j<N_j;j++){
            H_b = H_indices_sj[j];
            J_ij[i*N_j+j] = J_full[H_a*N+H_b];
	    // printf("Jij %f\n",J_ij[i*N_j+j]);
        }
    }
    return;
}

/* Find the size of the largest segment in site counts. */
int find_max_segment_size(int *psites, t_non *non){
    // analyse the projection to find the size of the largest segment
    int i, N, n_max;
    int *segment_sizes;
    N=non->singles;
    n_max = 0;
    // this array is longer than necessary in certain cases
    segment_sizes = (int *)calloc(N,sizeof(int));
    for (i=0;i<N;i++){
        segment_sizes[psites[i]] += 1;
    }
    for (i=0;i<N;i++){
        if (segment_sizes[i] > n_max){
	    n_max = segment_sizes[i];
	    }
    }

    free(segment_sizes);
    // return the max segment size
    return n_max;
}