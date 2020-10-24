#ifndef _CALC_2DIR_
#define _CALC_2DIR_
#include <mpi.h>

void calc_2DIR(t_non* non, int parentRank, int parentSize, int subRank, int subSize, MPI_Comm subComm, MPI_Comm rootComm);

#endif // _CALC_2DIR_
