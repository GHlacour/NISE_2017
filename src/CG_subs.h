#ifndef _CG_subs_
#define _CG_subs_

int CG_index(t_non *non,int seg_num,int alpha,int beta,int t1);
void normalize_DW(t_non *non,float *re,float *im,int samples);
void write_response_to_file(t_non *non,char fname[],float *im,float *re,int tmax);
void read_doorway_window_from_file(t_non *non,char fname[],float *im,float *re,int tmax);
void CG_doorway(t_non *non,float *re_doorway,float *im_doorway);
void CG_P_DA(t_non *non,float *P_DA, int N);
void CG_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
void CG_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
#endif /* _CG_subs_ */ 