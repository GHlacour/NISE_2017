#ifndef _CALC_2DIRraman_
#define _CALC_2DIRraman_
#include <mpi.h>

void calc_2DIRraman(t_non* non, int parentRank, int parentSize, int subRank, int subSize, MPI_Comm subComm, MPI_Comm rootComm);
void propagate_0(t_non* non, float* Hamil_i_e, int currentSample, int molPol, int t3, float* fr, float* fi,float** ft1r, float** ft1i,float** t1nr, float** t1ni, float* rightnr, float* rightni, float** rightrr, float** rightri,float* Anh);
void propagate_1(t_non* non, float* Hamil_i_e, float* fr, float* fi,float** ft1r, float** ft1i,float** t1nr, float** t1ni, float* rightnr, float* rightni, float** rightrr, float** rightri,float* Anh);

#endif // _CALC_2DIRraman_
