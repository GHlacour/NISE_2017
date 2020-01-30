#ifndef _MPI_SUBS_
#define _MPI_SUBS_

#include <mpi.h>
void print2D(char* filename, float** arrR, float** arrI, t_non* non, int sampleCount);
void calculateWorkset(t_non* non, int** workset, int* sampleCount, int* clusterCount);
void asyncWaitForMPI(MPI_Request requests[], int requestCount, int initialWaitingTime, int maxWaitingTime);

#endif // _MPI_SUBS_
