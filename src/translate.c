#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "translate.h"
#include "types.h"
#include "NISE_subs.h"
#include "readinput.h"

/* This code translate the Hamiltonian between different formats */
/* Existing formats:
   GROBIN
   GROASC
   MITASC
   SPECTRON
   MITTXT
   SKIBIN    // New GROBIN including diagonal anharmonicity files
*/
/* The code allows keeping only a subset of states and it allows
   isotope labeling. Input and output format can be the same. */

int main(int argc,char *argv[]){
  t_trans *tdat;
  t_files *files;
  t_ham *ham;
  t_modify *modify;
  int i;

  tdat=(t_trans *)calloc(1,sizeof(t_trans));
  files=(t_files *)calloc(1,sizeof(t_files));
  modify=(t_modify *)calloc(1,sizeof(t_modify));
  ham=(t_ham *)calloc(1,sizeof(t_ham));
  transinput(argc,argv,tdat,modify);
  //  printf("1\n");
  openInOutput(tdat,files);
  //  printf("2\n");
  initializeHam(tdat,ham);
  //  printf("3\n");
  for (i=0;i<tdat->length;i++){
    //    printf("R %d\n",i);
    readInp(tdat,ham,files);
    // Modify one exciton Hamiltonian if needed
    //    printf("A\n");
    if (tdat->modify==1) modifyHam(tdat,modify,ham);
    //    Construct doubly excited
    //    printf("C %d\n",i);
    //    printf("B\n");
    if (strcmp(tdat->skipDoubles,"Doubles")){
      constructHf(tdat,ham);
    }
    //    printf("W %d\n",i);
    writeOut(tdat,ham,files,i);
    //    printf("C\n");
    if (tdat->modify==1){
      //      freeHam(ham);
      //      printf("D\n");
      revMod(tdat,modify);
      //      printf("E\n");
      //      initializeHam(tdat,ham);
    }
  }
  //  printf("X\n");
  //  closeInOutput(files,tdat);
  freeHam(ham);
  return 0;
}

void readInp(t_trans *tdat,t_ham *ham,t_files *FH){
  int i,j,k,l,N;
  int control;
  int t;
  char buffer[4096];
  char fbuf[256];
  char dum[16];
  char *pStatus;
  size_t fLength;
  float f;
  int index,x;
  // Read GROBIN format
  if (strcmp(tdat->inputFormat,"GROBIN")==0){
    fread(&ham->t,sizeof(int),1,FH->IE);
    N=tdat->singles*(tdat->singles+1)/2;
    fread(&ham->He[0],sizeof(float),N,FH->IE);
    //    for (k=0;k<tdat->singles;k++){
      //      for (l=k;l<tdat->singles;l++){
      //      l=k;
      //	index=Sindex(k,l,tdat->singles);
	//	printf("%d %d %d %f\n",k,l,index,ham->He[index]);
	//      }
    //    }
    N=tdat->doubles*(tdat->doubles+1)/2;
    //    printf("N %d",tdat->doubles);
    fread(&ham->Hf[0],sizeof(float),N,FH->IE);

    fread(&ham->t,sizeof(int),1,FH->ID);
    N=tdat->singles*3;
    fread(&ham->mu_ge[0],sizeof(float),N,FH->ID);
    N=tdat->singles*tdat->doubles*3;
    fread(&ham->mu_ef[0],sizeof(float),N,FH->ID);
    if (FH->IALP!=NULL){
      fread(&ham->t,sizeof(int),1,FH->IALP);
      N=tdat->singles*3;
      fread(&ham->alpha[0],sizeof(float),N,FH->IALP);
    }
  }

  // Read SKIBIN format
  if (strcmp(tdat->inputFormat,"SKIBIN")==0){
    fread(&ham->t,sizeof(int),1,FH->IE);
    N=tdat->singles*(tdat->singles+1)/2;
    fread(&ham->He[0],sizeof(float),N,FH->IE);

    fread(&ham->t,sizeof(int),1,FH->ID);
    N=tdat->singles*3;
    fread(&ham->mu_ge[0],sizeof(float),N,FH->ID);
    
    fread(&ham->t,sizeof(int),1,FH->IA);
    N=tdat->singles;
    fread(&ham->Anh[0],sizeof(float),N,FH->IA);

    fread(&ham->t,sizeof(int),1,FH->IO);
    N=tdat->singles*3;
    fread(&ham->Over[0],sizeof(float),N,FH->IO);
  }

  // Read MITASC
  if (strcmp(tdat->inputFormat,"MITASC")==0){
    pStatus=fgets(&buffer[0],sizeof(buffer),FH->IE);
    j=0;
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<tdat->singles;l++){
	fLength = strcspn(&buffer[j], "\t");
	strcpy(fbuf,"                ");
	strncpy(fbuf,&buffer[j],fLength);
	f=atof(fbuf);
	if (k<=l){
	  index=l+tdat->singles*k-(k*(k+1)/2);
	  ham->He[index]=f;
	}
	j=j+fLength+1;    
      }
      j++;
    }
    pStatus=fgets(&buffer[0],sizeof(buffer),FH->IDx);
    j=0;
    for (k=0;k<tdat->singles;k++){
      fLength = strcspn(&buffer[j], "\t");
      strcpy(fbuf,"                ");
      strncpy(fbuf,&buffer[j],fLength);
      f=atof(fbuf);
      ham->mu_ge[k]=f;
      j=j+fLength+1;
    }
    pStatus=fgets(&buffer[0],sizeof(buffer),FH->IDy);
    j=0;
    for (k=0;k<tdat->singles;k++){
      fLength = strcspn(&buffer[j], "\t");
      strcpy(fbuf,"                ");
      strncpy(fbuf,&buffer[j],fLength);
      f=atof(fbuf);
      ham->mu_ge[k+tdat->singles]=f;
      j=j+fLength+1;
    }
    pStatus=fgets(&buffer[0],sizeof(buffer),FH->IDz);
    j=0;
    for (k=0;k<tdat->singles;k++){
      fLength = strcspn(&buffer[j], "\t");
      strcpy(fbuf,"                ");
      strncpy(fbuf,&buffer[j],fLength);
      f=atof(fbuf);
      ham->mu_ge[k+tdat->singles*2]=f;
      j=j+fLength+1;
    }
  }

  // Read GROASC
  if (strcmp(tdat->inputFormat,"GROASC")==0){
    fscanf(FH->IE,"%d",&t);
    for (k=0;k<tdat->singles;k++){
      for (l=k;l<tdat->singles;l++){
	fscanf(FH->IE,"%f",&f);
	if (k<=l){
	  index=l+tdat->singles*k-(k*(k+1)/2);
	  ham->He[index]=f;
	}	
      }
    }
    fscanf(FH->ID,"%d",&t);
    for (x=0;x<3;x++){
      for (k=0;k<tdat->singles;k++){
	fscanf(FH->ID,"%f",&f);
	index=x*tdat->singles+k;
	ham->mu_ge[index]=f;
      }
    }
    if (FH->IALP!=NULL){
      fscanf(FH->IALP,"%d",&t);
      for (x=0;x<3;x++){
	for (k=0;k<tdat->singles;k++){
	  fscanf(FH->IALP,"%f",&f);
	  index=x*tdat->singles+k;
	  ham->alpha[index]=f;
	}
      }
    }
  }

  // Read MITTXT
  if (strcmp(tdat->inputFormat,"MITTXT")==0){
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<tdat->singles;l++){
	fscanf(FH->IE,"%f",&f);
	if (k<=l){
	  index=l+tdat->singles*k-(k*(k+1)/2);
	  ham->He[index]=f;
	}	
      }
    }
    for (k=0;k<tdat->singles;k++){
      for (x=0;x<3;x++){
	fscanf(FH->ID,"%f",&f);
	index=x*tdat->singles+k;
	ham->mu_ge[index]=f;
      }
    }
  }

  // Read SPECTRON
  if (strcmp(tdat->inputFormat,"SPECTRON")==0){
    fscanf(FH->IE,"%s %d",dum,&t);
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<=k;l++){
	fscanf(FH->IE,"%f",&f);
	if (k>=l){
	  index=k+tdat->singles*l-(l*(l+1)/2);
	  ham->He[index]=f;
	}	
      }
    }
    fscanf(FH->ID,"%s %d",dum,&t);
    for (k=0;k<tdat->singles;k++){
      for (x=0;x<3;x++){
	fscanf(FH->ID,"%f",&f);
	index=x*tdat->singles+k;
	ham->mu_ge[index]=f;
      }
    }
  }
  return;
}

void writeOut(t_trans *tdat,t_ham *ham,t_files *FH,int snapshot){
  int i,j,k,l,N;
  int t;
  int index,x;
  // Write GROBIN format
  if (strcmp(tdat->outputFormat,"GROBIN")==0){
    
    //    printf("X %f\n",ham->He[1]);
    
    fwrite(&ham->t,sizeof(int),1,FH->OE);
    N=tdat->singles*(tdat->singles+1)/2;
    fwrite(&ham->He[0],sizeof(float),N,FH->OE);
    //    printf("X %f\n",ham->He[1]);
    if (strcmp(tdat->skipDoubles,"Doubles")){
      N=tdat->doubles*(tdat->doubles+1)/2;
      //    printf(" %d\n",tdat->doubles);
      fwrite(&ham->Hf[0],sizeof(float),N,FH->OE);
    }
    //    for (i=0;i<tdat->doubles;i++){
    //      printf("%d %f %f\n",i,ham->Hf[Sindex(i,0,tdat->doubles)],ham->mu_ef[Sindex(i,0,tdat->doubles)]);
    //    }

    fwrite(&ham->t,sizeof(int),1,FH->OD);
    N=tdat->singles*3;
    fwrite(&ham->mu_ge[0],sizeof(float),N,FH->OD);
    if (strcmp(tdat->skipDoubles,"Doubles")){
      N=tdat->singles*tdat->doubles*3;
      fwrite(&ham->mu_ef[0],sizeof(float),N,FH->OD);
    }
    if (FH->IALP!=NULL){
      fwrite(&ham->t,sizeof(int),1,FH->OALP);
      N=tdat->singles*3;
      fwrite(&ham->alpha[0],sizeof(float),N,FH->OALP);
    }
  }

  // Write SKIBIN format
  if (strcmp(tdat->outputFormat,"SKIBIN")==0){    
    fwrite(&ham->t,sizeof(int),1,FH->OE);
    N=tdat->singles*(tdat->singles+1)/2;
    fwrite(&ham->He[0],sizeof(float),N,FH->OE);
    
    fwrite(&ham->t,sizeof(int),1,FH->OD);
    N=tdat->singles*3;
    fwrite(&ham->mu_ge[0],sizeof(float),N,FH->OD);

    fwrite(&ham->t,sizeof(int),1,FH->OA);
    N=tdat->singles;
    fwrite(&ham->Anh[0],sizeof(float),N,FH->OA);
    
    fwrite(&ham->t,sizeof(int),1,FH->OO);
    N=tdat->singles*3;
    fwrite(&ham->Over[0],sizeof(float),N,FH->OO);
  }


  // Write MITASC
  if (strcmp(tdat->outputFormat,"MITASC")==0){
    //    printf("x\n");
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<tdat->singles;l++){
	index=Sindex(k,l,tdat->singles);
	fprintf(FH->OE,"%f\t",ham->He[index]);
      }
    }
    fprintf(FH->OE,"\n");
    for (k=0;k<tdat->singles;k++){
      fprintf(FH->ODx,"%f\t",ham->mu_ge[0*tdat->singles+k]);
      fprintf(FH->ODy,"%f\t",ham->mu_ge[1*tdat->singles+k]);
      fprintf(FH->ODz,"%f\t",ham->mu_ge[2*tdat->singles+k]);
    }
    fprintf(FH->ODx,"\n");
    fprintf(FH->ODy,"\n");
    fprintf(FH->ODz,"\n");
  }
  // Write GROASC
  if (strcmp(tdat->outputFormat,"GROASC")==0){
    fprintf(FH->OE,"%d ",0);
    for (k=0;k<tdat->singles;k++){
      for (l=k;l<tdat->singles;l++){
	index=Sindex(k,l,tdat->singles);
	fprintf(FH->OE,"%f ",ham->He[index]);
	//	if (k==l) printf("%d %d %d %f\n",k,l,index,ham->He[index]);
      }
    }
    fprintf(FH->OE,"\n");
    fprintf(FH->OD,"%d ",0);
    for (x=0;x<3;x++){
      for (k=0;k<tdat->singles;k++){
	fprintf(FH->OD,"%f ",ham->mu_ge[x*tdat->singles+k]);
      }
    }
    fprintf(FH->OD,"\n");

    if (FH->IALP!=NULL){
      fprintf(FH->OALP,"%d ",0);
      for (x=0;x<3;x++){
	for (k=0;k<tdat->singles;k++){
	  fprintf(FH->OALP,"%f ",ham->alpha[x*tdat->singles+k]);
	}
      }
      fprintf(FH->OALP,"\n");
    }
  }

  // Write MITTXT
  if (strcmp(tdat->outputFormat,"MITTXT")==0){
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<tdat->singles;l++){
	index=Sindex(k,l,tdat->singles);
	fprintf(FH->OE,"%f ",ham->He[index]);
	//	if (k==l) printf("%d %d %d %f\n",k,l,index,ham->He[index]);
      }
      fprintf(FH->OE,"\n");
    }
    fprintf(FH->OE,"\n");
    for (k=0;k<tdat->singles;k++){
      for (x=0;x<3;x++){
	fprintf(FH->OD,"%f ",ham->mu_ge[x*tdat->singles+k]);
      }
      fprintf(FH->OD,"\n");
    }
    fprintf(FH->OD,"\n");
  }
  
  // Write SPECTRON
  if (strcmp(tdat->outputFormat,"SPECTRON")==0){
    fprintf(FH->OE,"SNAPSHOT %d\n",snapshot+1);
    for (k=0;k<tdat->singles;k++){
      for (l=0;l<=k;l++){
	index=Sindex(k,l,tdat->singles);
	fprintf(FH->OE,"%f ",ham->He[index]);
	//	if (k==l) printf("%d %d %d %f\n",k,l,index,ham->He[index]);
      }
      fprintf(FH->OE,"\n");
    }
    //    fprintf(FH->OE,"\n");

    fprintf(FH->OD,"SNAPSHOT %d\n",snapshot+1);
    for (k=0;k<tdat->singles;k++){
      for (x=0;x<3;x++){
	fprintf(FH->OD,"%f ",ham->mu_ge[x*tdat->singles+k]);
      }
      fprintf(FH->OD,"\n");
    }
    //    fprintf(FH->OD,"\n");
  }

  return;
}

void constructHf(t_trans *tdat,t_ham *ham){
  int *iA,*iB;
  int i,j,k,x;
  int Nf,Ne;
  float A,w;
  Ne=tdat->singles,Nf=tdat->doubles;
  A=tdat->anharmonicity;

  if (Nf!=Ne*(Ne+1)/2 && (tdat->modify!=1)) return; // Doubles already exist
  Nf=Ne*(Ne+1)/2;

  // Doubly excited indexing
  iA=(int *)calloc(Nf,sizeof(int));
  iB=(int *)calloc(Nf,sizeof(int));
  k=0;
  for (i=0;i<Ne;i++){
    for (j=i;j<Ne;j++){
      iA[k]=i,iB[k]=j;
      k++;
    }
  }// 9/1-2007 debug
  /*  for (i=0;i<Ne;i++){
    for (j=i;j<Ne;j++){
      k=i+j*Ne-j*(j+1)/2;
      iA[k]=i,iB[k]=j;
    }
    }*/


  // Construct doubly excited Hamiltonian
  for (i=0;i<Nf;i++){
    for (j=i;j<Nf;j++){
      w=0;
      // Set diagonal values
      if (i==j){
	if (iA[i]==iB[i]){
	  w=2*ham->He[Sindex(iA[i],iA[i],Ne)]-A;
	  //	  printf("%d %d %f\n",iA[i],iB[i],w);
	} else {
	  w=ham->He[Sindex(iA[i],iA[i],Ne)]+ham->He[Sindex(iB[i],iB[i],Ne)];
	  //	  printf("%d %d %f\n",iA[i],iB[i],w);
	}
      } else {
	// Set offdiagonal values
	// AA-AB or AA-BA
	if (iA[i]==iB[i]){
	  if (iA[j]==iA[i] || iB[j]==iA[i]){
	    w=sqrt2*ham->He[Sindex(iA[j],iB[j],Ne)];
	  }
	}
	// AB-AA or BA-AA
	if (iA[j]==iB[j]){
	  if (iA[i]==iA[j] || iB[i]==iA[j]){
	    w=sqrt2*ham->He[Sindex(iA[i],iB[i],Ne)];
	  }
	}
	if (iA[i]!=iB[i] && iA[j]!=iB[j]){
	  // AB-AC
	  if (iA[i]==iA[j]){
	    w=ham->He[Sindex(iB[i],iB[j],Ne)];
	    // BA-CA
	  } else if (iB[i]==iB[j]) {
	    w=ham->He[Sindex(iA[i],iA[j],Ne)];
	    // AB-CA
	  } else if (iA[i]==iB[j]){
	    w=ham->He[Sindex(iA[j],iB[i],Ne)];
	    // BA-AC
	  } else if (iA[j]==iB[i]){
	    w=ham->He[Sindex(iA[i],iB[j],Ne)];
	  }
	}
	//	printf("%d %d %d %d %f\n",iA[i],iB[i],iA[j],iB[j],w);
      }
      ham->Hf[Sindex(i,j,Nf)]=w;
    }
  }

  // Doubly excited dipole
  for (x=0;x<3;x++){ 
    for (i=0;i<Ne;i++){
      for (j=0;j<Nf;j++){	 
	w=0;
	// i->2i
	if (iA[j]==iB[j] && i==iA[j]){
	  w=sqrt2*ham->mu_ge[x*Ne+i];
	  // i->i j or i->j i
	} else if (iA[j]!=iB[j]) {
	  if (iA[j]==i) w=ham->mu_ge[Ne*x+iB[j]];
	  if (iB[j]==i) w=ham->mu_ge[Ne*x+iA[j]];
	}
	ham->mu_ef[Nf*Ne*x+Nf*i+j]=w;
      }
    }
  }
  // End doubly excited dipoles
  return;
}

void transinput(int argc,char *argv[],t_trans *tdat,t_modify *modify){
  char inputFName[256];
  FILE *inputFile;
  char *pStatus;
  char Buffer[256];
  size_t LabelLength;
  char *pValue;
  int control;
  if (argc<2){
    printf("Specify input file name on command line!\n");
    printf("Program terminated!\n");
    exit(-1);
  } else {
    strcpy(&inputFName[0], argv[1]);
    printf("Using input file '%s'.\n",inputFName);
  }

  // Open input file
  inputFile=fopen(inputFName, "r");
  if (inputFile == NULL) {
    printf("File not found!\n");
    exit(-1);
  }

  // Defaults
  tdat->anharmonicity=16;
  tdat->dipole12=1.41421356237;
  tdat->doubles=-1;
  tdat->length=1;
  tdat->modify=0; // No modifications

  control=0;
  // Read input data
  do {
    pStatus = fgets(&Buffer[0],sizeof(Buffer),inputFile);
    if (pStatus == NULL) {
      break;
    }
    
    // Compute LabelLength
    LabelLength = strcspn(&Buffer[0], " ");

    if (keyWordS("InputEnergy",Buffer,tdat->inputEnergy,LabelLength)==1) continue;

    if (keyWordS("OutputEnergy",Buffer,tdat->outputEnergy,LabelLength)==1) continue;

    if (keyWordS("InputDipole",Buffer,tdat->inputDipole,LabelLength)==1) continue;

    if (keyWordS("OutputDipole",Buffer,tdat->outputDipole,LabelLength)==1) continue;
    if (keyWordS("InputAlpha",Buffer,tdat->inputAlpha,LabelLength)==1) continue;

    if (keyWordS("OutputAlpha",Buffer,tdat->outputAlpha,LabelLength)==1) continue;
    if (keyWordS("InputAnharm",Buffer,tdat->inputAnh,LabelLength)==1) continue;

    if (keyWordS("OutputAnharm",Buffer,tdat->outputAnh,LabelLength)==1) continue;

    if (keyWordS("InputOverto",Buffer,tdat->inputOver,LabelLength)==1) continue;

    if (keyWordS("OutputOverto",Buffer,tdat->outputOver,LabelLength)==1) continue;
    if (keyWordS("InputDipoleX",Buffer,tdat->inputDipolex,LabelLength)==1) continue;

    if (keyWordS("OutputDipoleX",Buffer,tdat->outputDipolex,LabelLength)==1) continue;
    if (keyWordS("InputDipoleY",Buffer,tdat->inputDipoley,LabelLength)==1) continue;

    if (keyWordS("OutputDipoleY",Buffer,tdat->outputDipoley,LabelLength)==1) continue;
    if (keyWordS("InputDipoleZ",Buffer,tdat->inputDipolez,LabelLength)==1) continue;

    if (keyWordS("OutputDipoleZ",Buffer,tdat->outputDipolez,LabelLength)==1) continue;

    if (keyWordI("Singles",Buffer,&tdat->singles,LabelLength)==1) continue;
    if (keyWordI("Length",Buffer,&tdat->length,LabelLength)==1) continue;
    if (keyWordI("Doubles",Buffer,&tdat->doubles,LabelLength)==1) continue;
    
    if (keyWordF("Anharmonicity",Buffer,&tdat->anharmonicity,LabelLength)==1) continue;
    if (keyWordS("InputFormat",Buffer,tdat->inputFormat,LabelLength)==1) continue;

    if (keyWordS("OutputFormat",Buffer,tdat->outputFormat,LabelLength)==1) continue;
   // Skip Doubles
    if (keyWordS("Skip",Buffer,tdat->skipDoubles,LabelLength)==1) continue;
    // Modifications
    if (keyWordModify("Modify",Buffer,&tdat->modify,LabelLength,modify,inputFile,tdat->singles,tdat)==1) continue;
    
  } while (1==1);
  fclose(inputFile);

  // Standard peptide
  if (tdat->doubles==-1) tdat->doubles=(tdat->singles*(tdat->singles+1))/2;
  
  // Test format validity
  if (strcmp(tdat->inputFormat,"GROBIN")==0) control++;
  if (strcmp(tdat->inputFormat,"GROASC")==0) control++;
  if (strcmp(tdat->inputFormat,"MITASC")==0) control++;
  if (strcmp(tdat->inputFormat,"SPECTRON")==0) control++;
  if (strcmp(tdat->inputFormat,"MITTXT")==0) control++;
  if (strcmp(tdat->inputFormat,"SKIBIN")==0) control++;
  if (control!=1){
    printf("Input format %s unknown.\n",tdat->inputFormat);
    exit(0);
  }
  if (strcmp(tdat->outputFormat,"GROBIN")==0) control+=2;
  if (strcmp(tdat->outputFormat,"GROASC")==0) control+=2;
  if (strcmp(tdat->outputFormat,"MITASC")==0) control+=2;
  if (strcmp(tdat->outputFormat,"SPECTRON")==0) control+=2;
  if (strcmp(tdat->outputFormat,"MITTXT")==0) control+=2;
  if (strcmp(tdat->outputFormat,"SKIBIN")==0) control+=2;
  if (control!=3){
    printf("Output format %s unknown.\n",tdat->outputFormat);
    exit(0);
  }

  return;
}

void openInOutput(t_trans *tdat,t_files *HANDLES){ 

  // Input

  // GROBIN format
  if (strcmp(tdat->inputFormat,"GROBIN")==0){
    HANDLES->IE=fopen(tdat->inputEnergy,"rb");
    HANDLES->ID=fopen(tdat->inputDipole,"rb");
    HANDLES->IALP=fopen(tdat->inputAlpha,"rb");
  }

  // SKIBIN format
  if (strcmp(tdat->inputFormat,"SKIBIN")==0){
    HANDLES->IE=fopen(tdat->inputEnergy,"rb");
    HANDLES->ID=fopen(tdat->inputDipole,"rb");
    HANDLES->IA=fopen(tdat->inputAnh,"rb");
    HANDLES->IO=fopen(tdat->inputOver,"rb");
  }

  //  printf("T6\n");
  // GROASC/SPECTRON format
  if (strcmp(tdat->inputFormat,"GROASC")==0 || (strcmp(tdat->inputFormat,"SPECTRON")==0) || (strcmp(tdat->inputFormat,"MITTXT")==0)){
    HANDLES->IE=fopen(tdat->inputEnergy,"r");
    HANDLES->ID=fopen(tdat->inputDipole,"r");
    if (strcmp(tdat->inputFormat,"GROASC")==0){
      HANDLES->IALP=fopen(tdat->inputAlpha,"r");
    }
  }

  //  printf("T5\n");
  // MITASC format
  if (strcmp(tdat->inputFormat,"MITASC")==0){
    HANDLES->IE=fopen(tdat->inputEnergy,"r");
    HANDLES->IDx=fopen(tdat->inputDipolex,"r");
    HANDLES->IDy=fopen(tdat->inputDipoley,"r");
    HANDLES->IDz=fopen(tdat->inputDipolez,"r");
  }


  // Output

  // GROBIN
  if (strcmp(tdat->outputFormat,"GROBIN")==0){
    HANDLES->OE=fopen(tdat->outputEnergy,"wb");
    HANDLES->OD=fopen(tdat->outputDipole,"wb");
    if (HANDLES->IALP!=NULL) HANDLES->OALP=fopen(tdat->outputAlpha,"wb");
  }

  // SKIBIN
  if (strcmp(tdat->outputFormat,"SKIBIN")==0){
    HANDLES->OE=fopen(tdat->outputEnergy,"wb");
    HANDLES->OD=fopen(tdat->outputDipole,"wb");
    HANDLES->OA=fopen(tdat->outputAnh,"wb");
    HANDLES->OO=fopen(tdat->outputOver,"wb");
  }

  // GROASC/SPECTRON format
  if (strcmp(tdat->outputFormat,"GROASC")==0 || strcmp(tdat->outputFormat,"SPECTRON")==0 || strcmp(tdat->outputFormat,"MITTXT")==0){
    HANDLES->OE=fopen(tdat->outputEnergy,"w");
    HANDLES->OD=fopen(tdat->outputDipole,"w");
    //    if (HANDLES->OD==NULL) printf("q\n");
    if (strcmp(tdat->outputFormat,"GROASC")==0){
      if (HANDLES->IALP!=NULL) HANDLES->OALP=fopen(tdat->outputAlpha,"w");
    }
  }

  // MITASC
  if (strcmp(tdat->outputFormat,"MITASC")==0){
    HANDLES->OE=fopen(tdat->outputEnergy,"w");
    HANDLES->ODx=fopen(tdat->outputDipolex,"w");
    HANDLES->ODy=fopen(tdat->outputDipoley,"w");
    HANDLES->ODz=fopen(tdat->outputDipolez,"w");
  }
  //  printf("T4\n");
  if (HANDLES->IE==NULL){
    printf("Problem opening energy input file\n");
    exit(0);
  }
  if (HANDLES->OE==NULL){
    printf("Problem opening energy output file\n");
    exit(0);
  }
  //  printf("T3\n");
  if (strcmp(tdat->inputFormat,"GROBIN")==0 || strcmp(tdat->inputFormat,"GROASC")==0 || strcmp(tdat->inputFormat,"SPECTRON")==0 || strcmp(tdat->inputFormat,"MITTXT")==0|| strcmp(tdat->inputFormat,"SKIBIN")==0){  
    if (HANDLES->ID==NULL){
      printf("Problem opening dipole input file\n");
      exit(0);
    }
  } else {
    if (HANDLES->IDx==NULL || HANDLES->IDy==NULL || HANDLES->IDz==NULL){
      printf("Problem opening dipole input file\n");
      exit(0);
    }
  }
  //  printf("T2\n");
  if (strcmp(tdat->outputFormat,"GROBIN")==0 || strcmp(tdat->outputFormat,"GROASC")==0 || strcmp(tdat->outputFormat,"SPECTRON")==0 || strcmp(tdat->outputFormat,"MITTXT")==0|| strcmp(tdat->outputFormat,"SKIBIN")==0){
    if (HANDLES->OD==NULL){
      printf("Problem opening dipole output file\n");
      exit(0);
    }
  } else {
    if (HANDLES->ODx==NULL || HANDLES->ODy==NULL || HANDLES->ODz==NULL){
      printf("Problem opening dipole output file\n");
      exit(0);      
    }
  }
  if (HANDLES->IALP==NULL){
    printf("Transition polarizability file not opened!\n");
    //    exit(0);
  }
  //  printf("T1\n");
  return;
}

void closeInOutput(t_files *files,t_trans *tdat){
  fclose(files->IE);
  fclose(files->OE);
  if (strcmp(tdat->inputFormat,"GROBIN")==0 || strcmp(tdat->inputFormat,"GROASC")==0 || strcmp(tdat->inputFormat,"SPECTRON")==0 || strcmp(tdat->inputFormat,"MITTXT")==0){
    fclose(files->ID);
  } else {
    fclose(files->IDx),fclose(files->IDy),fclose(files->IDz);
  }
  if (strcmp(tdat->outputFormat,"GROBIN")==0 || strcmp(tdat->outputFormat,"GROASC")==0 || strcmp(tdat->inputFormat,"SPECTRON")==0 || strcmp(tdat->inputFormat,"MITTXT")==0 ){
    fclose(files->OD);
  } else {
    fclose(files->ODx),fclose(files->ODy),fclose(files->ODz);
  }
  if (strcmp(tdat->outputFormat,"GROBIN")==0){
    fclose(files->ID),fclose(files->OD);
    fclose(files->IA),fclose(files->OA);
    fclose(files->IO),fclose(files->OO);
  }
  if (files->IALP!=NULL){
    if (strcmp(tdat->inputFormat,"GROBIN")==0){
      fclose(files->IALP);
      fclose(files->OALP);
    }
  }
  return;
}

void initializeHam(t_trans *tdat,t_ham *ham){
  int Nf;
  Nf=tdat->singles*(tdat->singles+1)/2;
  ham->He=(float *)calloc(tdat->singles*(tdat->singles+1)/2,sizeof(float));
  ham->Hf=(float *)calloc(Nf*(Nf+1)/2,sizeof(float));
  ham->mu_ge=(float *)calloc(tdat->singles*3,sizeof(float));
  ham->mu_ef=(float *)calloc(tdat->singles*Nf*3,sizeof(float));
  ham->Anh=(float *)calloc(tdat->singles,sizeof(float));
  ham->Over=(float *)calloc(tdat->singles*3,sizeof(float));
  ham->alpha=(float *)calloc(tdat->singles*3,sizeof(float));
  return;
}

void freeHam(t_ham *ham){
  free(ham->He);
  free(ham->Hf);
  free(ham->mu_ge);
  free(ham->mu_ef);
  free(ham->Anh);
  free(ham->Over);
  free(ham->alpha);
}

/*int Sindex(int a,int b,int N){
  int ind;
  if (a>b){
    ind=a+N*b-(b*(b+1)/2);
  } else {
    ind=b+N*a-(a*(a+1)/2);
  }
  return ind;
}
// Index routine
//int Sindex(int a,int b,int N){
//  int ind;
//    if (a>b){
//    ind=a+N*b-(b*(b+1)/2);
//  } else {
//    ind=b+N*a-(a*(a+1)/2);
//    }
//  ind=a+N*b;
//  return ind;
//}
*/

// Read shift input
int keyWordModify(char *keyWord,char *Buffer,int *ivalue,size_t LabelLength,t_modify *modify,FILE *inputFile,int N,t_trans *tdat){
  char *pValue;
  char dummy[256];
  char Buf[256];
  int i;
  char *pStatus;
  int debug=1;

  if (!strncmp(&Buffer[0],&keyWord[0],6)){
    printf("%s:\n",keyWord);
 
    *ivalue=1;
    
    // Read modification information
    // Select
    pStatus = fgets(&Buf[0],sizeof(Buf),inputFile);
    printf("%s\n",Buf);

    if (!strncmp(&Buf[0],"Select",6)){
      sscanf(Buf,"%s %d",dummy,&modify->singles);
      modify->doubles=0;
      if (strcmp(tdat->skipDoubles,"Doubles")){
	modify->doubles=(modify->singles*(modify->singles+1))/2;
      }
      modify->select=(int *)calloc(modify->singles,sizeof(int));
      modify->label=(int *)calloc(modify->singles,sizeof(int));
      modify->shift=(float *)calloc(modify->singles,sizeof(float));	
      for (i=0;i<modify->singles;i++){
	if (modify->singles!=N){
	  fscanf(inputFile,"%d ",&modify->select[i]);
	} else {
	  modify->select[i]=i;
	  //	  printf("MS %d %d\n",i,modify->select[i]);
	}
      }
    } else {
      printf("Select keyword not found after Modify!\n");
      exit(-1);
    }
    // Label
    pStatus = fgets(&Buf[0],sizeof(Buf),inputFile);
    if (!strncmp(&Buf[0],"Label",5)){
      for (i=0;i<modify->singles;i++){
	fscanf(inputFile,"%d ",&modify->label[i]);
      }
    } else {
      printf("Label keyword not found in Modify!\n");
      exit(-1);
    }
   
    // Shift
    pStatus = fgets(&Buf[0],sizeof(Buf),inputFile);
    if (!strncmp(&Buf[0],"Shift",5)){
      for (i=0;i<modify->singles;i++){
	fscanf(inputFile,"%f ",&modify->shift[i]);
      }
    } else {
      printf("Shift keyword not found in Modify!\n");
      exit(-1);
    }
    return 1;
  }
  
  return 0;
}

void modifyHam(t_trans *tdat,t_modify *modify,t_ham *ham){
  float *He,*mu_ge;
  int i,j,x;
  int dim;
  int index;
  float H;
  int swap;
  int debug=1;
  float *Anh,*Over,*alpha;

  dim=(modify->singles*(modify->singles+1))/2;
  He=(float *)calloc(dim,sizeof(float));
  mu_ge=(float *)calloc(modify->singles*3,sizeof(float));
  alpha=(float *)calloc(modify->singles*3,sizeof(float));
  Anh=(float *)calloc(modify->singles,sizeof(float));
  Over=(float *)calloc(modify->singles*3,sizeof(float));

  // Build new Hamiltonian
  index=0;
  //  printf("Q\n");
  for (i=0;i<modify->singles;i++){
    for (j=i;j<modify->singles;j++){
      //      printf("T %d %d %d %d %d\n",i,j,modify->select[i],modify->select[j],Sindex(modify->select[i],modify->select[j],tdat->singles));
      //      printf("B %d %d %f\n",i,j,ham->He[Sindex(modify->select[i],modify->select[j],tdat->singles)]);
      H=ham->He[Sindex(modify->select[i],modify->select[j],tdat->singles)];
      // Shift and Label
      if (i==j){
	//	printf("P %d %f %d %d \n",i,H,modify->label[i],modify->select[i]);
	H+=modify->shift[i];
	if (modify->label[i]==1){
	  H-=41; // C13 label
	} else if (modify->label[i]==2){
	  H-=60; // O18 label
	}
	//	printf("D %d %f\n",i,H);
      }
      //      printf("W\n");
      He[index]=H;
      index++;
    }
    //    printf("X\n");
    Anh[i]=ham->Anh[modify->select[i]];
    // Dipoles
    for (x=0;x<3;x++){
      mu_ge[modify->singles*x+i]=ham->mu_ge[tdat->singles*x+modify->select[i]];
      Over[modify->singles*x+i]=ham->Over[tdat->singles*x+modify->select[i]];
      alpha[modify->singles*x+i]=ham->alpha[tdat->singles*x+modify->select[i]];
    }
  }

  //  printf("R\n");
  // Redefine Hamiltonian
  //  free(ham->He);
  //  free(ham->Hf);
  //  free(ham->mu_ge);
  //  free(ham->mu_ef);
  swap=tdat->singles;
  tdat->singles=modify->singles;
  modify->singles=swap;
  swap=tdat->doubles;
  tdat->doubles=modify->doubles;
  modify->doubles=swap;
  //  ham->He=(float *)calloc(tdat->singles*(tdat->singles+1)/2,sizeof(float));
  //  ham->Hf=(float *)calloc(tdat->doubles*(tdat->doubles+1)/2,sizeof(float));
  //  ham->mu_ge=(float *)calloc(tdat->singles*3,sizeof(float));
  //  ham->mu_ef=(float *)calloc(tdat->singles*tdat->doubles*3,sizeof(float));
  // Copy Hamiltonians
 
  index=0;
  for (i=0;i<tdat->singles;i++){
    for (j=i;j<tdat->singles;j++){
      ham->He[index]=He[index];

      //      printf("E %d %d %f\n",i,j,ham->He[index]);
      
      index++;
    }
    ham->Anh[i]=Anh[i];
    for (x=0;x<3;x++){
      ham->mu_ge[x*tdat->singles+i]=mu_ge[x*tdat->singles+i];
      ham->Over[x*tdat->singles+i]=Over[x*tdat->singles+i];
      ham->alpha[x*tdat->singles+i]=alpha[x*tdat->singles+i];
    }
  }
  
  free(He);
  free(mu_ge);
  free(Anh);
  free(Over);
  free(alpha);
  
}

void revMod(t_trans *tdat,t_modify *modify){
  int swap;
  swap=tdat->singles;
  tdat->singles=modify->singles;
  modify->singles=swap;
  swap=tdat->doubles;
  tdat->doubles=modify->doubles;
  modify->doubles=swap;
}
