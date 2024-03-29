#ifndef _SFG_
#define _SFG_
void sfg(t_non *non);
// We need SSP, and PPP
void calc_SFG_SSP(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha);
void calc_SFG_PPP(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci,float *alpha);

#endif // _SFG_
