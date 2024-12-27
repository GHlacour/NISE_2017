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

    int N; /* Number of chromophores */
    int TT; /* Total timesteps */
    int i,j,k,l;
    int segments;

    N=non->singles;
    TT=non->length;
    sigma=non->inhomogen;
    if (sigma<=0) sigma=1;

    printf("\n");
    printf("The Spectral Density will be smoothened with a %f cm-1\n",sigma);
    printf("wide gaussian function. Change with Inhomogeneous keyword.\n\n");
    
    /* Reserve Memory */
    SD_matrix=(float *)calloc(N*TT,sizeof(float));
    omega=(float *)calloc(TT,sizeof(float));
    avHam=(float *)calloc(N*(N+1)/2,sizeof(float));
    c=(float *)calloc(N*N,sizeof(float));
    e=(float *)calloc(N,sizeof(float));
    Redfield=(float *)calloc(N*N,sizeof(float));

    /* Read in Spectral Density */
    input=fopen("SpectralDensity.dat","r");
    if (input==NULL){
        printf("Problem opening SpectralDensity.dat file.\n");
        printf("Did you run Technique Correlate?\n");
        exit(1);
    }    
    for (i=0;i<TT;i++){
        fscanf(input,"%f",&omega[i]);
        for (j=0;j<N;j++){
            fscanf(input,"%f",&SD_matrix[j*TT+i]);
        }
    }
    fclose(input);
    

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
        for (j=0;j<N;j++){
            if (i!=j){
                deltaE=e[i]-e[j];
                /* Loop over bath spectral density */
                for (k=0;k<N;k++){
                    product=c[i*N+k]*c[j*N+k]*c[i*N+k]*c[j*N+k];
                    /* Loop over positive frequencies */
                    for (l=0;l<TT;l++){
                        cont=exp(-0.5*(deltaE-omega[l])*(deltaE-omega[l])/sigma/sigma)/sigma/sq2pi;
                        Redfield[i*N+j]+=cont*product*SD_matrix[k*TT+l];
                    }
                    /* Loop over negative frequencies */
                    for (l=1;l<TT;l++){
                        cont=exp(-0.5*(deltaE+omega[l])*(deltaE+omega[l])/sigma/sigma)/sigma/sq2pi;
                        cont=cont*exp(omega[l]/non->temperature/k_B);
                        Redfield[i*N+j]+=cont*product*SD_matrix[k*TT+l];
                    }
                } 
            }
        }
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
            fprintf(output,"%f ",c[i*N+j]);
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
