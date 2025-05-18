#ifndef _calc_FD_CG_2DES_
#define _calc_FD_CG_2DES_

void calc_FD_CG_2DES(t_non *non);
void call_final_FD_CG_2DES(
  t_non *non,int pro_dim,float *re_doorway,float *im_doorway,
  float *re_window_SE, float *im_window_SE,float *re_window_GB, float *im_window_GB,
  float *re_window_EA, float *im_window_EA, double *Q1, double *Q2);
void FD_CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
    float *re_window_SE,float *im_window_SE,float *re_window_GB, float *im_window_GB,
    float *re_window_EA,float *im_window_EA,float *P_DA,int N, char *waittime,int wfile, double *Q1, double *Q2);
void read_in_QF_for_FD_CG_2DES(t_non *non, double *Q1, double *Q2, int segments);
#endif /* _calc_CG_2DES_ */     
