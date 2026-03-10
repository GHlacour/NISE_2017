#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "types.h"
#include "NISE_subs.h"
#include "project.h"
#include "randomlib.h"
#include "util/asprintf.h"

// Read projection input
void read_shift(t_non *non){
    int i;
    FILE *shiftFile;
    
    // Check if file name is given
    if (non->singleShiftFName[0] == '\0') {
        non->SingleShiftSites = 0;
        return;
    }

    // Open single shift file if present
    shiftFile=fopen(non->singleShiftFName,"r");
    if (shiftFile == NULL) {
        printf(RED "Single shift file %s not found!\n" RESET,non->singleShiftFName);
        return;
    }
    
    // Read shifts
    if (fscanf(shiftFile,"%d",&non->SingleShiftSites)!=1){
        printf(RED "Could not read number of single shifts site from %s\n" RESET,non->singleShiftFName);
        exit(-1);
    }

    // Allocate memory
    non->SingleShiftSite=(int *)calloc(non->SingleShiftSites,sizeof(int));
    non->SingleShift=(float *)calloc(2*non->SingleShiftSites,sizeof(float));

    // Read shifts
    for (i=0;i<non->SingleShiftSites;i++){
        fscanf(shiftFile,"%d %f %f",&non->SingleShiftSite[i],&non->SingleShift[i*2],&non->SingleShift[i*2+1]);
    }

    fclose(shiftFile);
    return;
}

// Free memory allocated for single shifts
void free_shift(t_non *non){
    if (non->SingleShiftSites>0){
        free(non->SingleShiftSite);
        free(non->SingleShift);
    }
    return;
}

// Apply single site shifts to Hamiltonian
void apply_singleshift(t_non *non,float *Hamil_i_e){
    int i;
    int index;
    for (i=0;i<non->SingleShiftSites;i++){
        index=Sindex(non->SingleShiftSite[i],non->SingleShiftSite[i], non->singles);
        Hamil_i_e[index]=Hamil_i_e[index]*non->SingleShift[2*i]+non->SingleShift[2*i+1];
    }
    return;
}