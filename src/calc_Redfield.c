#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"
#include "propagate.h"
#include "read_trajectory.h"

/* This implementation follows the paper of Mino Yang and Graham Fleming, Chemical Physics 282 (2002) 163–180 Eq. 33 */
void calc_Redfield(t_non *non){
    FILE *input;
    FILE *output;

    float *SD_matrix;
    float *avHam;
    float *omega;
    float dummy;
    float *c,*e;
    float *Redfield;
    float deltaE;
    float cont,product;
    float sigma;
    float domega;

    int N; /* Number of chromophores */
    int TT; /* Total timesteps */
    int i,j,k,l;
    int segments;
    int warnings;

    N=non->singles;
    TT=non->length;

    warnings=0;

    /* Reserve Memory */
    SD_matrix=(float *)calloc(N*TT,sizeof(float));
    omega=(float *)calloc(TT,sizeof(float));
    avHam=(float *)calloc(N*(N+1)/2,sizeof(float));
    c=(float *)calloc(N*N,sizeof(float));
    e=(float *)calloc(N,sizeof(float));
    Redfield=(float *)calloc(N*N,sizeof(float));

    /* Read in Spectral Density */
    input = fopen("SpectralDensity.dat", "r");
    if (input == NULL) {
       printf("Problem opening SpectralDensity.dat file.\n");
       printf("Did you run Technique Correlate?\n");
       exit(1);
    }

    i = 0;
    while (!feof(input)) {
       int result = fscanf(input, "%f", &omega[i]);
       if (result != 1) break;  // Break if we couldn't read a float

       for (j = 0; j < N; j++) {
          result = fscanf(input, "%f", &SD_matrix[j*TT + i]);
          if (result != 1) {
             printf("Unexpected end of file or format error\n");
             printf("While reading SpectralDensity.dat.\n");
             fclose(input);
             exit(1);
          }
       }
       i++;

       if (i >= TT) {
          printf("Warning: Maximum array size reached\n");
          printf("Was a different trajectory length specified for the\n");
          printf("Correlation and Redfield techniques?\n");
          printf("Largest frequency read from SpectralDensity file\n");
          printf("is %f cm-1.\n",omega[i-1]);
          break;
      }
    }

    fclose(input);
    
    domega=omega[1]-omega[0];
    printf("\nRead spectral density from the SpectralDensity.dat file.\n");
    printf("The frequency gap used is %f cm-1.\n",domega);

    sigma=non->inhomogen;
    if (sigma<=0) sigma=domega*10;

    printf("\n");
    printf("The Spectral Density will be smoothened with a %f cm-1\n",sigma);
    printf("wide Gaussian function. This can be changed with the Inhomogeneous keyword.\n\n");
    if (sigma<5*domega){
        printf("Using smoothning smaller or close to the resolution of the spectral density.\n");
        printf("may result in inaccurate results. Please, consider a larger value!\n");
    }

    /* Read in Average Hamiltonian */
    input=fopen("Av_Hamiltonian.txt","r");
    if (input==NULL){
        printf("Problem opening Av_Hamiltonian.txt file.\n");
        printf("Did you run Technique Analyse?\n");
        exit(1);
    }
    fscanf(input,"%f",&dummy);
    for (i=0;i<N;i++){
        for (j=i;j<N;j++){
            fscanf(input,"%f",&avHam[Sindex(i,j,N)]);
        }
    }
    fclose(input);
    
    /* Apply projection if defined */
    segments=project_dim(non);
    printf("Found %d segments.\n",segments);
    if (segments>1) multi_projection_Hamiltonian(avHam,non); /* Not tested */
    /* Diagonalize Hamiltonian */
    build_diag_H(avHam,c,e,N);

    printf("Using %f K as temperature.\n",non->temperature);

    /* Calculate Redfield rates */
    for (i=0;i<N;i++){
#pragma omp parallel for shared(warnings) private(product,cont,deltaE)
        for (j=0;j<N;j++){
            if (i!=j){
                deltaE=e[j]-e[i];
                if (fabs(deltaE)>omega[(int)(TT/2)]) {
                    if (warnings==0){
                        printf("Energy gap larger than spectral density cutoff deteted!\n");
                        printf("Please, veryfy that the spectral density is low enough for cutoff.\n");
                        printf("Create a trajectory with smaller time intervals if needed.\n\n");
                    }
                    warnings++;
                }
                /* Loop over bath spectral density */
                for (k=0;k<N;k++){
                    product=c[k*N+i]*c[k*N+j]*c[k*N+i]*c[k*N+j]; // i and j are eigenstates, k is site
                    /* Loop over positive frequencies */
                    for (l=0;l<TT;l++){
                        cont=exp(-0.5*(deltaE-omega[l])*(deltaE-omega[l])/sigma/sigma)/sigma/sq2pi;
                        Redfield[i*N+j]+=cont*product*SD_matrix[k*TT+l];
                    }
                    /* Loop over negative frequencies */
                    for (l=1;l<TT;l++){
                        cont=exp(-0.5*(deltaE+omega[l])*(deltaE+omega[l])/sigma/sigma)/sigma/sq2pi;
                        /* Suppression of upwards energy transfer */
                        cont=cont*exp(-omega[l]/non->temperature/k_B);
                        Redfield[i*N+j]+=cont*product*SD_matrix[k*TT+l];
                    }
                } 
            }
        }
        printf("Finished transfer from average eigenstate %d\n",i);
    }

   if (warnings>0){
       printf("There were %d instances of too large gap detected in total.\n\n",warnings);
   } 

    /* Convert from cm-1 to ps and Set diagonal rates */
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            if (i!=j){
                Redfield[i*N+j]=Redfield[i*N+j]*icm2ifs*1000*twoPi*twoPi; // Where does the last 2pi come from?
                Redfield[j*N+j]-=Redfield[i*N+j];
            }
        }
    }

    /* Write Redfield rates to file */
    output=fopen("RedfieldRates.dat","w");
    if (output==NULL){
        printf("Problem opening RedfieldRates.dat file.\n");
        printf("Disk full?\n");
        exit(1);
    }

    for (i=0;i<N;i++){
        //fprintf(output,"%f ",e[i]);
        for (j=0;j<N;j++){
            fprintf(output,"%f ",Redfield[i*N+j]); 
        }
        fprintf(output,"\n");
    }
    fclose(output);


    /* Write Redfield states to file */
    output=fopen("RedfieldStates.dat","w");
    if (output==NULL){
        printf("Problem opening RedfieldStates.dat file.\n");
        printf("Disk full?\n");
        exit(1);
    }
    
    for (i=0;i<N;i++){
        fprintf(output,"%f ",e[i]);
        for (j=0;j<N;j++){
            fprintf(output,"%f ",c[j*N+i]); // j runs over sites, i over eigenstates
        }
        fprintf(output,"\n");
    }
    fclose(output);

    /*Free used memory */
    free(SD_matrix),free(avHam),free(omega),free(c),free(e),free(Redfield);
    printf("----------------------------------------------\n");
    printf(" Redfield calculation succesfully completed\n");
    printf("----------------------------------------------\n\n");

    return;
}
