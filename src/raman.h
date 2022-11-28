#ifndef _RAMAN_
#define _RAMAN_
void raman(t_non *non);
void calc_Ram_VV(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha);
void calc_Ram_VH(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha);

#endif // _RAMAN_
