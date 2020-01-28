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
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "MPI_subs.h"
#include <stdarg.h>
#include "mpi.h"

// Print results to the corresponding files
void print2D(char* filename, float** arrR, float** arrI, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        const int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            arrR[t3][t1] /= sampleCount;
            arrI[t3][t1] /= sampleCount;
            // Divide boarder points with 2 for FFT
            if (t3==0) arrR[t3][t1] /= 2, arrI[t3][t1] /= 2;
            fprintf(out, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                arrR[t3][t1], arrI[t3][t1]);
        }
    }
    fclose(out);
}

void calculateWorkset(t_non* non, int** workset, int* sampleCount, int* clusterCount) {
    // Open clustering file if necessary
    FILE* Cfile;
    if (non->cluster != -1) {
        Cfile = fopen("Cluster.bin", "rb");
        if (Cfile == NULL) {
            printf("Cluster option was activated but no Cluster.bin file provided.\n");
            printf("Please, provide cluster file or remove Cluster keyword from\n");
            printf("input file.\n");
            exit(0);
        }
    }
    *clusterCount = 0;

    // Determine number of samples
    const int totalSampleCount = (non->length - non->tmax1 - non->tmax2 - non->tmax3 - 1) / non->sample + 1;
    if (totalSampleCount > 0) {
        printf("Making %d samples!\n", totalSampleCount);
    }
    else {
        printf("Insufficient data to calculate spectrum.\n");
        printf("Please, lower max times or provide longer\n");
        printf("trajectory.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (non->end == 0) non->end = totalSampleCount;
    *sampleCount = non->end - non->begin;
    if (*sampleCount == 0) {
        // Avoid dividing by zero
        *sampleCount = 1;
    }

    // Allocate workset array
    *workset = calloc(2 * 21 * *sampleCount, sizeof(int)); // two integers per work item (sample + polDir), 21 polDirs per sample

    // Loop over samples, fill work set array of things to do
    int currentWorkItem = 0;
    for(int currentSample = non->begin; currentSample < non->end; currentSample++) {
        // Check clustering
        if(non->cluster != -1) {
            int tj = currentSample * non->sample + non->tmax1;

            int currentCluster;
            if (read_cluster(non, tj, &currentCluster, Cfile) != -1) {
                printf("Cluster trajectory file to short, could not fill buffer!!!\n");
                printf("ITIME %d\n", tj);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // If we do not need to consider current cluster, skip calculation
            if(non->cluster != currentCluster) {
                log_item("Skipping sample %d, incorrect cluster!\n", currentSample);
                (*sampleCount)--;
                continue;
            }

            (*clusterCount)++;
        }

        // Set work items
        for(int molPol = 0; molPol < 21; molPol++) {
            (*workset)[currentWorkItem * 2] = currentSample;
            (*workset)[currentWorkItem * 2 + 1] = molPol;
            currentWorkItem++;
        }
    }

     // Now we have a workset array with all sample/molPol combinations that need to be calculated.
    // Furthermore, sampleCount is adjusted to reflect the actual number of samples to consider when clustering

    if (non->cluster != -1) {
        fclose(Cfile);
    }
}

// Function to wait until all MPI requests in the given array have finished
void asyncWaitForMPI(MPI_Request requests[], int requestCount, int initialWaitingTime, int maxWaitingTime) {
    int waitingTime = initialWaitingTime;
    int finished = 0;
    MPI_Testall(requestCount, requests, &finished, MPI_STATUSES_IGNORE);
    while(!finished) {
        #ifdef _WIN32
            Sleep(waitingTime);
        #else
            nanosleep((const struct timespec[]){{0, waitingTime * 1000000L}}, NULL);
        #endif

        if(waitingTime < maxWaitingTime) waitingTime *= 2;

        MPI_Testall(requestCount, requests, &finished, MPI_STATUSES_IGNORE);
    }
}

