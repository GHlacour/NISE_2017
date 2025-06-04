#ifndef _calc_CG_2DES_
#define _calc_CG_2DES_

void calc_CG_2DES(t_non *non);
void call_final_CG_2DES(
  t_non *non,int pro_dim,float *re_doorway,float *im_doorway,
  float *re_window_SE, float *im_window_SE,float *re_window_GB, float *im_window_GB,
  float *re_window_EA, float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
    float *re_window_SE,float *im_window_SE,float *re_window_GB, float *im_window_GB,
    float *re_window_EA,float *im_window_EA,float *P_DA,int N, char *waittime,int wfile);
#endif /* _calc_CG_2DES_ */     
