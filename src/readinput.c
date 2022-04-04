#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types.h"
#include "readinput.h"

/* Read the input file */
void readInput(int argc, char* argv[], t_non* non) {
    char inputFName[256];
    FILE* inputFile;
    char* pStatus;
    char Buffer[256];
    size_t LabelLength;
    char* pValue;
    int control;
    char prop[256];

    // Defaults
    non->interpol = 1;
    non->begin = 0;
    non->end = 0;
    non->ts = 5;
    non->anharmonicity = 0;
    non->couplingcut = 0;
    non->temperature = 300;
    non->cluster = -1; // Average over all snapshots no clusters
    non->fft = 0;
    non->printLevel = 1; // Set to standard print level
    non->tmax1=256;
    non->tmax2=0;
    non->tmax3=256;
    non->homogen=0.0;
    non->inhomogen=0.0;
    sprintf(non->basis, "Local");
    sprintf(non->hamiltonian, "Full");

    if (argc < 2) {
        printf(RED "Specify input file name on command line!\n");
        printf("Program terminated!\n" RESET);
        exit(-1);
    }
    else {
        strcpy(&inputFName[0], argv[1]);
        printf("Using input file '%s'.\n", inputFName);
    }

    // Open input file
    inputFile = fopen(inputFName, "r");
    if (inputFile == NULL) {
        printf("File not found!\n");
        exit(-1);
    }

    control = 0;

    // Read input data
    do {
        pStatus = fgets(&Buffer[0], sizeof(Buffer), inputFile);
        if (pStatus == NULL) {
            break;
        }

        // Compute LabelLength
        LabelLength = strcspn(&Buffer[0], " ");

        // Propagation keyword
        if (keyWordS("Propagation", Buffer, prop, LabelLength) == 1) continue;

        // Hamiltonian file keyword
        if (keyWordS("Hamiltonianfile", Buffer, non->energyFName, LabelLength) == 1) continue;

        // Dipole file keyword
        if (keyWordS("Dipolefile", Buffer, non->dipoleFName, LabelLength) == 1) continue;

        // Transisition polarizability file keyword
        if (keyWordS("Alphafile", Buffer, non->alphaFName, LabelLength) == 1) continue;

        // Anharmonicity file keyword
        if (keyWordS("Anharmonicfile", Buffer, non->anharFName, LabelLength) == 1) continue;

        // Overtone dipole file keyword
        if (keyWordS("Overtonedipolefile", Buffer, non->overdipFName, LabelLength) == 1) continue;

        // Position file keyword
        if (keyWordS("Positionfile", Buffer, non->positionFName, LabelLength) == 1) continue;
        // PDB file keyword
        if (keyWordS("PDBfile", Buffer, non->pdbFName, LabelLength) == 1) continue;

        // Coupling file keyword
        if (keyWordS("Couplingfile", Buffer, non->couplingFName, LabelLength) == 1) continue;

        // Read Trajectory length
        if (keyWordI("Length", Buffer, &non->length, LabelLength) == 1) continue;

        // Read Begin point of Trajectory (for perfect parallel splitting)
        //    if (keyWordI("BeginSample",Buffer,&non->begin,LabelLength)==1) continue;

        // Read Samplerate
        if (keyWordI("Samplerate", Buffer, &non->sample, LabelLength) == 1) continue;

        // Read Beginpoint
        if (keyWordI("BeginPoint", Buffer, &non->begin, LabelLength) == 1) continue;

        // Read Endpoint
        if (keyWordI("EndPoint", Buffer, &non->end, LabelLength) == 1) continue;

        // Read Lifetime
        if (keyWordF("Lifetime", Buffer, &non->lifetime, LabelLength) == 1) continue;

        // Read Homogeneous Lifetime
        if (keyWordF("Homogeneous", Buffer, &non->homogen, LabelLength) == 1) continue;

        // Read Inhomogeneous Lifetime
        if (keyWordF("Inhomogeneous", Buffer, &non->inhomogen, LabelLength) == 1) continue;

        // Read timestep
        if (keyWordF("Timestep", Buffer, &non->deltat, LabelLength) == 1) continue;

        // Read timestep
        if (keyWordF("Anharmonicity", Buffer, &non->anharmonicity, LabelLength) == 1) continue;

        // Read length for possible zeropadding
        if (keyWordI("FFT", Buffer, &non->fft, LabelLength) == 1) continue;

        // Read Dephasingtime
        //    if (keyWordF("DephasingTime",Buffer,&non->dephasing,LabelLength)==1) continue;

        // Read treshhold for sparce matrices
        if (keyWordF("Threshold", Buffer, &non->thres, LabelLength) == 1) continue;

        // Read coupling cutoff for sparce matrices with coupling prop. scheme
        if (keyWordF("Couplingcut", Buffer, &non->couplingcut, LabelLength) == 1) continue;

        // Read Temperature (For Luminescence
        if (keyWordF("Temperature", Buffer, &non->temperature, LabelLength) == 1) continue;

        // Read maxtimes
        if (keyWord3I("RunTimes", Buffer, &non->tmax1, &non->tmax2, &non->tmax3, LabelLength) == 1) continue;

        // Read integration steps
        if (keyWordI("Integrationsteps", Buffer, &non->is, LabelLength) == 1) continue;

        // Read trotter steps
        if (keyWordI("Trotter", Buffer, &non->ts, LabelLength) == 1) continue;

        // Read interpolation steps
        if (keyWordI("Interpolation", Buffer, &non->interpol, LabelLength) == 1) continue;

        // Read cluster info if applicable
        if (keyWordI("Cluster", Buffer, &non->cluster, LabelLength) == 1) continue;

        // Read single excited states
        if (keyWordI("Singles", Buffer, &non->singles, LabelLength) == 1) continue;

        // Read double excited states
        if (keyWordI("Doubles", Buffer, &non->doubles, LabelLength) == 1) continue;

        // Read minfrequencies
        if (keyWord3F("MinFrequencies", Buffer, &non->min1, &non->min2, &non->min3, LabelLength) == 1) continue;

        // Read maxfrequencies
        if (keyWord3F("MaxFrequencies", Buffer, &non->max1, &non->max2, &non->max3, LabelLength) == 1) continue;

        // Read static
        if (keyWord3F("Static", Buffer, &non->statstart, &non->statend, &non->statstep, LabelLength) == 1) continue;

        // Read technique
        if (keyWordS("Technique", Buffer, non->technique, LabelLength) == 1) continue;

        // Read Hamiltonian Type
        if (keyWordS("HamiltonianType", Buffer, non->hamiltonian, LabelLength) == 1) continue;

        // Read basis (for population calculation)
        if (keyWordS("Basis", Buffer, non->basis, LabelLength) == 1) continue;

        // Read projection (for calculating spectra from selected sites)
        if (keyWordProject("Projection", Buffer, LabelLength, &non->Npsites, inputFile, non->singles, non) == 1)
            continue;

        // Read double excited states
        if (keyWordI("PrintLevel", Buffer, &non->printLevel, LabelLength) == 1) continue;
       

    }
    while (1 == 1);
    fclose(inputFile);

    non->dt1 = 1, non->dt2 = 1, non->dt3 = 1;
    // Set length of linear response function
    non->tmax = non->tmax1;
    /* Is the length large enough? */
    if (non->length < non->tmax1 + non->tmax2 + non->tmax3) {
        printf(RED "The trajectory length is too small!\n");
        printf("It must be longer than %d snapshots.\n" RESET, non->tmax1 + non->tmax2 + non->tmax3);
        exit(0);
    }

    // Check RunTimes keyword setting
    if (non->tmax1 == 0) {
        printf(RED "First runtime variable is zero.\n");
        printf("You need to specify the RunTimes keyword!\n" RESET);
        exit(0);
    }

    // Decide propagation scheme
    non->propagation = 0;
    if (!strcmp(prop, "Coupling")) {
        non->propagation = 1;
        printf("\nUsing propagation scheme 'Coupling'!\n");
        printf("Coupling cutoff %f effective during t1 and t3.\n\n",
               non->couplingcut);
    }
    if (!strcmp(prop, "Diagonal")) {
        non->propagation = 2;
        printf("\nUsing propagation with full diagonalization!\n\n");
        printf(RED "Presently NOT implemented. Use sparse with no cutoff!\n" RESET);
        exit(0);
    }

    if (non->propagation == 0) {
        printf("Rescaling threshold with factor %g. (dt/hbar)**2\n",
               (non->deltat * icm2ifs * twoPi / non->ts) * (non->deltat * icm2ifs * twoPi / non->ts));
        non->thres *= (non->deltat * icm2ifs * twoPi / non->ts) * (non->deltat * icm2ifs * twoPi / non->ts);
        if (non->thres > 0.1) {
            printf("Unrealistic value for threshold!\n");
            exit(0);
        }
        printf("Scaled value for threshold: %g.\n", non->thres);
        printf("Will neglect elements of time-evolution operator\n");
        printf("smaller than this value.\n");
    }

    if ((!strcmp(non->technique, "2DIR")) || (!strcmp(non->technique, "GBIR")) || (!strcmp(non->technique, "SEIR")) ||
     (!strcmp(non->technique, "EAIR")) || (!strcmp(non->technique, "noEAIR")) || (!strcmp(non->technique, "2DUVvis")) ||
     (!strcmp(non->technique, "EAUVvis")) || (!strcmp(non->technique, "noEAUVvis")) || (!strcmp(non->technique, "2DSFG")) ||
     (!strcmp(non->technique, "SEUVvis")) || (!strcmp(non->technique, "GBUVvis"))) {
        printf("\nThe waiting time will be %f fs.\n\n", non->tmax2 * non->deltat);
    }

    // Prepare static spec
    non->statsteps = rint(fabs((non->statend - non->statstart) / non->statstep));

    return;
}

// Read string input
int keyWordS(char* keyWord, char* Buffer, char* value, size_t LabelLength) {
    char* pValue;
    char dummy[256];
    if (!strncmp(&Buffer[0], &keyWord[0], LabelLength)) {
        printf("%s:", keyWord);
        pValue = &Buffer[LabelLength];
        while (*pValue == ' ') {
            pValue++;
        }

        sscanf(Buffer, "%s %s", dummy, value);
        printf(" %s\n", value);
        return 1;
    }
    return 0;
}

// Read integer input
int keyWordI(char* keyWord, char* Buffer, int* ivalue, size_t LabelLength) {
    char* pValue;
    char dummy[256];
    if (!strncmp(&Buffer[0], &keyWord[0], LabelLength)) {
        printf("%s:", keyWord);
        pValue = &Buffer[LabelLength];
        while (*pValue == ' ') {
            pValue++;
        }

        *ivalue = atoi(pValue);
        printf(" %d\n", *ivalue);
        return 1;
    }
    return 0;
}

// Read triple integer input
int keyWord3I(char* keyWord, char* Buffer, int* i1, int* i2, int* i3, size_t LabelLength) {
    char* pValue;
    char dummy[256];
    if (!strncmp(&Buffer[0], &keyWord[0], LabelLength)) {
        printf("%s:", keyWord);
        pValue = &Buffer[LabelLength];
        while (*pValue == ' ') {
            pValue++;
        }

        sscanf(Buffer, "%s %d %d %d", dummy, i1, i2, i3);
        printf(" %d %d %d\n", *i1, *i2, *i3);
        return 1;
    }
    return 0;
}

// Read float input
int keyWordF(char* keyWord, char* Buffer, float* ivalue, size_t LabelLength) {
    char* pValue;
    char dummy[256];
    if (!strncmp(&Buffer[0], &keyWord[0], LabelLength)) {
        printf("%s:", keyWord);
        pValue = &Buffer[LabelLength];
        while (*pValue == ' ') {
            pValue++;
        }

        *ivalue = atof(pValue);
        printf(" %f\n", *ivalue);
        return 1;
    }
    return 0;
}

// Read triple double input
int keyWord3F(char* keyWord, char* Buffer, float* f1, float* f2, float* f3, size_t LabelLength) {
    char* pValue;
    char dummy[256];
    if (!strncmp(&Buffer[0], &keyWord[0], LabelLength)) {
        printf("%s:", keyWord);
        pValue = &Buffer[LabelLength];
        while (*pValue == ' ') {
            pValue++;
        }

        sscanf(Buffer, "%s %f %f %f", dummy, f1, f2, f3);
        printf(" %f %f %f\n", *f1, *f2, *f3);
        return 1;
    }
    return 0;
}

// Read projection input
int keyWordProject(char* keyWord, char* Buffer, size_t LabelLength, int* singles, FILE* inputFile, int N, t_non* non) {
    char* pValue;
    char dummy[256];
    char Buf[256];
    int i, NN, j;
    char* pStatus;
    FILE* projectFile;
    char projectFName[256];

    if (!strncmp(&Buffer[0], &keyWord[0], 6)) {
        printf("%s:\n", keyWord);

        // Read projection information
        // Select
        pStatus = fgets(&Buf[0], sizeof(Buf), inputFile);
        printf("%s", Buf);
        /* Proceed here if site numbers are provided directly in input file */
        if (!strncmp(&Buf[0], "Sites", 5)) {
            sscanf(Buf, "%s %d", dummy, &NN);
            //      printf("%d\n",NN);
            non->psites = (int *)calloc(N, sizeof(int));
            for (i = 0; i < NN; i++) {
                if (NN != N) {
                    fscanf(inputFile, "%d ", &j);
                    if (non->psites[j] ==1){
                       printf(RED "\n\nSite %d was defined twice for projection!\n",j);
                       printf("Check your list of sites in the input file.\n");
                       printf("Aborting calculation.\n" RESET);
                       exit(0);
                    }
                    non->psites[j] = 1;
                    printf("%d ", j);
                }
                else {
                    non->psites[i] = 1;
                    //	  printf("MS %d %d\n",i,modify->select[i]);
                }
            }
        /* Proceed here if site numbers are provided in a separate file */
        } else if (!strncmp(&Buf[0], "Projectfile", 11)) {
          pValue = &Buf[11];
          while (*pValue == ' ') {
            pValue++;
          }
          /* Read file name */
          sscanf(Buf, "%s %s", dummy, projectFName);
          printf("Reading sites to project on from %s\n", projectFName);
          /* Read projection file */
          projectFile=fopen(projectFName,"r");
          if (projectFile == NULL) {
             printf(RED "Projectfile %s not found!\n" RESET,projectFName);
             exit(-1);
          }

          fscanf(projectFile, "%d ", &NN);
          printf("Projecting on the following %d sites:\n",NN);
          non->psites = (int *)calloc(N, sizeof(int));
          for (i = 0; i < NN; i++) {
              if (NN != N) {
                  fscanf(projectFile, "%d ", &j);
                  if (non->psites[j] ==1){
                     printf(RED "\n\nSite %d was defined twice for projection!\n",j);
                     printf("Check your list of sites in the file: %s \n",projectFName);
                     printf("Aborting calculation.\n" RESET);
                     exit(0);
                  }
                  non->psites[j] = 1;
                  printf("%d ", j);
              }
              else {
                  non->psites[i] = 1;
     //             printf("MS %d %d\n",i,modify->select[i]);
              }
          }       
          fclose(projectFile);
        } else {
            printf(RED "Neither Sites nor Projectfile keyword found after Project!\n" RESET);
            exit(-1);
        }
        printf("\n");

        *singles = NN;
        return 1;
    }
    return 0;
}
