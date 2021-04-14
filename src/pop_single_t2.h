#ifndef _POP_T2_
#define _POP_T2_

void pop_print(char* filename, float* pop, t_non* non, int sampleCount);
void propagate_NISE(t_non* non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci);
void pop_single_t2(t_non* non);

#endif // _POP_T2_