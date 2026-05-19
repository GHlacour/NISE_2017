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

    // If there are shifts difined read them and check if they are valid
    if (non->SingleShiftSites>0){
        printf("\nApplying single shifts to %d sites.\n",non->SingleShiftSites);
    } else {
        printf("No single shifts defined.\n");
    }
    // Read shifts
    for (i=0;i<non->SingleShiftSites;i++){
        fscanf(shiftFile,"%d %f %f",&non->SingleShiftSite[i],&non->SingleShift[i*2],&non->SingleShift[i*2+1]);
        // Check if site number exists
        if (non->SingleShiftSite[i]<0 || non->SingleShiftSite[i]>=non->singles){
            printf(RED "Invalid site index for single shift: %d. Must be between 0 and %d.\n" RESET,non->SingleShiftSite[i],non->singles-1);
            exit(-1);
        }
        // Check if scaling factor is valid
        if (non->SingleShift[i*2]<=0){
            printf(RED "Invalid scaling factor for single shift: %f. Must be positive.\n" RESET,non->SingleShift[i*2]);
            exit(-1);
        }
        // Check if offset is valid
        if (fabs(non->SingleShift[i*2+1])>non->max1-non->min1){
            printf(YELLOW "Invalid offset for single shift: %f. Must be within the range +-%f.\n" RESET,non->SingleShift[i*2+1],non->max1-non->min1);
        }
        if (non->printLevel<2){
            printf("Applying single shift to site %d: multiply by %f and add %f.\n",non->SingleShiftSite[i],non->SingleShift[i*2],non->SingleShift[i*2+1]);
        }
    }
    // Add an extra line for better readability if single shifts are applied
    if (non->SingleShiftSites>0){
        printf("\n");
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
