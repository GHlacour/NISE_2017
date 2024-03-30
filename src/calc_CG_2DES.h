#ifndef _calc_CG_2DES_
#define _calc_CG_2DES_

void calc_CG_2DES(t_non *non);
void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway);
void CG_2DES_P_DA(t_non *non,float *P_DA, int N);
void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_2DES_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_doorway,float *im_doorway,
                                      float *re_window_SE,float *im_window_SE,
                                      float *re_window_GB, float *im_window_GB,
                                      float *re_window_EA,float *im_window_EA,
				      float *P_DA,int N);
#endif /* _calc_CG_2DES_ */     
