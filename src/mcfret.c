#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"

void mcfret(t_non *non){
    /* Define variables and arrays */
    int N;
    float *re_S_Abs,*im_S_Abs;
    float *re_S_Emi,*im_S_Emi;

    re_S_Abs=(float *)calloc(N*N*non->tmax,sizeof(float));
    im_S_Abs=(float *)calloc(N*N*non->tmax,sizeof(float));
    re_S_Emi=(float *)calloc(N*N*non->tmax,sizeof(float));
    im_S_Emi=(float *)calloc(N*N*non->tmax,sizeof(float));

    /* Call the MCFRET Routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Autodetect")) || (!strcmp(non->technique, "MCFRET-Absorption"))
        || (!strcmp(non->technique, "MCFRET-Emission")) || (!strcmp(non->technique, "MCFRET-Coupling")) || (!strcmp(non->technique,
        "MCFRET-Rate")) || (!strcmp(non->technique, "MCFRET-Analyse")) ) {
        printf("Perfroming the MCFRET calculation.\n");
    }

    /* Call the MCFRET Response for Absorption matrix */
    if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Absorption")  ){
        printf("Calculation Absortion Matrix for MCFRET.\n");
        mcfret_response_function(non,re_S_Abs,im_S_Abs,0);
    }

    /* Call the MCFRET Response for Absorption matrix */
    if (!strcmp(non->technique, "MCFRET") || !strcmp(non->technique, "MCFRET-Emission")  ){
        printf("Calculation Emssion Matrix for MCFRET.\n");
        mcfret_response_function(non,re_S_Emi,im_S_Emi,1);
    }

    free(re_S_Abs),free(im_S_Abs);
    free(re_S_Emi),free(im_S_Emi);
    return;
}

void mcfret_autodetect(t_non *non, float treshold);

void mcfret_response_function(t_non *non, float *re_S_1, float *im_S_1,int emission){

    return;
}

void mcfret_coupling(t_non *non);

void mcfret_rate(t_non *non);

void mcfret_validate(t_non *non);
void mcfret_analyse(t_non *non);
