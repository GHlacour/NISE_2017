#ifndef _calc_CG_2DES_
#define _calc_CG_2DES_
void calc_CG_2DES(t_non *non);
void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway);
//void CG_2DES_P_DA(t_non *non,float *P_DA);
void CG_2DES_P_DA(t_non *non,float *P_DA,float* K, float* P0, int N);
void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_2DES_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_2DES_pa,float *im_2DES_NR_pa,float *im_2DES_R_pa,
                          float *re_2DES_pe,float *im_2DES_NR_pe,float *im_2DES_R_pe,
                          float *re_2DES_cr,float *im_2DES_NR_cr,float *im_2DES_R_cr);
void combine_CG_2DES(t_non *non,float *re_2DES_pa_sum,float *im_2DES_NR_pa_sum,float *im_2DES_R_pa_sum,
                          float *re_2DES_pe_sum,float *im_2DES_NR_pe_sum,float *im_2DES_R_pe_sum,
                          float *re_2DES_cr_sum,float *im_2DES_NR_cr_sum,float *im_2DES_R_cr_sum);
#endif /* _calc_CG_2DES_ */