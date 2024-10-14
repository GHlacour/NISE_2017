#ifndef _CORRELATE_
#define _CORRELATE_
void subtractMean(float* signal, int N);
void calculateCorrelation(float *input1, float *input2, float *output,float *SD, int N);

void calc_Correlation(t_non *non);
#endif