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
#include "calc_CG_2DES.h"
#include <stdarg.h>
#include "project.h"

/* This is the main routine which will control the flow of the
 * coarse grained 2DES calculation. It takes the general NISE
 * input as input.
 */

void calc_CG_2DES(t_non *non){
    printf("Hello world!\n");
    return;
}
void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway);
void CG_2DES_P_DA(t_non *non,float *P_DA);
void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_2DES_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_2DES,float *im_2DES);
void combine_CG_2DES(t_non *non,float *re_doorway,float *im_doorway,
    float *P_DA,float *re_window_GB,float *im_window_GB,
    float *re_window_SE,float *im_window_SE,float *re_window_EA,float *im_window_EA,
    float *re_2DES,float *im_2DES);
