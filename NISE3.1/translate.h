#ifndef _TRANSLATE_
#define _TRANSLATE_

//float sqrt2=1.41421356237;

typedef struct {
  char inputEnergy[256];
  char outputEnergy[256];
  char inputDipole[256];
  char outputDipole[256];
  char inputAlpha[256];
  char outputAlpha[256];
  char inputAnh[256];
  char outputAnh[256];
  char inputOver[256];
  char outputOver[256];
  char inputDipolex[256];
  char outputDipolex[256];
  char inputDipoley[256];
  char outputDipoley[256];
  char inputDipolez[256];
  char outputDipolez[256];
  int singles,doubles;
  char inputFormat[256];
  char outputFormat[256];
  float anharmonicity;
  float dipole12;
  int length;
  char skipDoubles[256];
  int modify;
} t_trans;

typedef struct {
  int singles,doubles;
  int *label;
  float *shift;
  int *select;
} t_modify;

typedef struct {
  FILE *IE,*ID,*OE,*OD; // GRO formats
  FILE *IA,*IO,*OA,*OO; // SKI format
  FILE *IALP,*OALP;
  FILE *IDx,*IDy,*IDz; // MIT format
  FILE *ODx,*ODy,*ODz;
} t_files;

typedef struct {
  float *He,*Hf;
  float *Anh,*Over;
  float *mu_ge,*mu_ef;
  float *alpha;
  int t;
} t_ham;

void readInp(t_trans *tdat,t_ham *ham,t_files *FH);
void writeOut(t_trans *tdat,t_ham *ham,t_files *FH,int snapshot);
void constructHf(t_trans *tdat,t_ham *ham);
void transinput(int argc,char *argv[],t_trans *tdat,t_modify *modify);
void openInOutput(t_trans *tdat,t_files *HANDLES);
void closeInOutput(t_files *files,t_trans *tdat);
void initializeHam(t_trans *tdat,t_ham *ham);
void freeHam(t_ham *ham);
//int Sindex(int a,int b,int N);
int keyWordModify(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength,t_modify *modify,FILE *inputFile,int N,t_trans *tdat);
void modifyHam(t_trans *tdat,t_modify *modify,t_ham *ham);
void revMod(t_trans *tdat,t_modify *modify);

#endif // _TRANSLATE_
