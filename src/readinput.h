#ifndef _READINPUT_
#define _READINPUT_

void readInput(int argc,char *argv[],t_non *non);
int keyWordS(char *keyWord,char *Buffer,char *value,size_t LabelLength);
int keyWordI(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength);
int keyWord3I(char *keyWord,char *Buffer,int *i1,int *i2,int *i3,size_t LabelLength);
int keyWordF(char *keyWord,char *Buffer,float *ivalue,size_t LabelLength);
int keyWord3F(char *keyWord,char *Buffer,float *f1,float *f2,float *f3,size_t LabelLength);
int keyWordProject(char *keyWord,char *Buffer,size_t LabelLength,int *singles,FILE *inputFile,int N,t_non *non);
#endif // _READINPUT_
