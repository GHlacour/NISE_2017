#ifndef _CALC_2DES_
#define _CALC_2DES_
#include <mpi.h>

void calc_2DES(t_non* non, int parentRank, int parentSize, int subRank, int subSize, MPI_Comm subComm, MPI_Comm rootComm);

#endif // _CALC_2DES_
