#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "types.h"
#include "nrutil.h"

#define round(x) ((x)>0?(long)(x+0.5):(long)(x-0.5))

/* This program takes the Fourier transform of kI and kII processes */
/* and adds up the results to get the absorptive 2D spectra. */
/* The process is identical for 2DIR and 2DUVvis type spectra. */
int main(int argc,char *argv[]){
  int fft;
  fftw_complex *fftIn,*fftOut;
  fftw_plan fftPlan;
  int i,j,k,l;
  float deltat;
  FILE *input;
  FILE *output;
  FILE *axisFH;
  float ti[4],rr,ir;
  char inputFName[256],timeFName[256],frequencyFName[256];
  char format[256];
  char Dummy[256];
  FILE *fileHandle;
  char *pStatus,*pValue;
  size_t LabelLength;
  char Buffer[256];
  int tv1,tv2,tv3;
  t_time t;
  t_rwa w;
  int index;
  float fixedtime;
  float static c_v=2.99792458e-5;/* Speed of light in cm/fs */
  float shift1,shift3;
  int form;
  float w1,w3;
  float sw;
  int pol;
  // Add
  char kIFName[256],kIIFName[256];
  char twoDFName[256];
  char pPFName[256];
  int p1I,p3I,p1II,p3II;
  float **kI,**kI1,**kI3;
  float **kII,**kII1,**kII3;
  float **IR,**IR1,**IR3;
  float **kIi;
  float **kIIi;
  float **IRi;
  float pumpProbe;
  float u,q;
  int p2D1,p2D3;
  float dw1,dw3;
  float dum;
  float rwma1,rwma3,rwmi1,rwmi3;
  float homo,inhomo;

  fixedtime=-1;

  // Read input
  if (argc!=2){
    printf("Specify input file name on command line!\n");
    printf("Program terminated!\n");
    exit(1);
  } else {
    strcpy(&inputFName[0], argv[1]);
    printf("Using input file '%s'.\n",inputFName);
  }
  
  // Open input file
  fileHandle=fopen(inputFName, "r");
  if (fileHandle == NULL) {
    printf("File not found!\n");
    exit(1);
  }

  form=0; // Default format is Dislin
  do {
    pStatus = fgets(&Buffer[0],sizeof(Buffer),fileHandle);
    if (pStatus == NULL) {
      break;
    }
    
    // Compute LabelLength
    LabelLength = strcspn(&Buffer[0], " ");

    // FFT keyword
    if (!strncmp(&Buffer[0],"FFT",LabelLength)){
      printf("FFT: ");
      pValue = &Buffer[LabelLength];
      while (*pValue == ' '){
	pValue++;
      }
      
      fft=atoi(pValue);
      printf("%d\n",fft);
      continue;
    }

    // Timestep keyword
    if (!strncmp(&Buffer[0],"Timestep",LabelLength)){
      printf("Timestep: ");
      pValue = &Buffer[LabelLength];
      while (*pValue == ' '){
        pValue++;
      }

      deltat=(float) atof(pValue);
      printf("%f fs.\n",deltat);
      continue;
    }

    // Homogen keyword
    if (!strncmp(&Buffer[0],"Homogen",LabelLength)){
      printf("Homogeneous Lifetime: ");
      pValue = &Buffer[LabelLength];
      while (*pValue == ' '){
        pValue++;
      }

      homo=(float) atof(pValue);
      printf("%f fs.\n",homo);
      continue;
    }

    // Inhomogen keyword
    if (!strncmp(&Buffer[0],"Inhomogen",LabelLength)){
      printf("Inhomogeneous Lifetime: ");
      pValue = &Buffer[LabelLength];
      while (*pValue == ' '){
        pValue++;
      }

      inhomo=(float) atof(pValue);
      printf("%f fs.\n",inhomo);
      continue;
    }

    // Format
    if (!strncmp(&Buffer[0],"Format",LabelLength)){
      printf("Format:");
      
      sscanf(Buffer,"%s %s",format,format);
      if (!strcmp(&format[0],"Dislin")) form=0;
      if (!strcmp(&format[0],"Matlab")) form=1;
      if (!strcmp(&format[0],"Gnuplot")) form=2;
      if (form==0) printf("Dislin\n");
      if (form==1) printf("Matlab\n");
      if (form==2) printf("Gnuplot\n");
      continue;
    }

    // Timevariables keyword
    tv1=1,tv2=3;
    tv3=2;
    /*    if (!strncmp(&Buffer[0],"Timevariables",LabelLength)){
      printf("Timevariables:");
      
      sscanf(Buffer,"%s %d %d",Dummy,&tv1,&tv2);

      printf(" t%d t%d\n",tv1,tv2);
      // tv1+tv2+tv3=1+2+3=6
      if (tv1==tv2){
	printf("The two time variables must be different!\n");
	exit(0);
      }
      if (tv1<1 || tv1>3 || tv2<1 || tv2>3){
	printf("The time variables must equal 1, 2 or 3!\n");
	exit(0);
      }
      tv3=6-tv1-tv2;
      printf("Fixed time: t%d\n",tv3); 
      continue;
      }*/

    // Fixedtime keyword
    /*    if (!strncmp(&Buffer[0],"Fixedtime",LabelLength)){
      printf("Fixedtime:");
      
      sscanf(Buffer,"%s %f",Dummy,&fixedtime);

      printf(" %f fs\n",fixedtime);
      continue;
      }*/

    // RunTimes keyword
    if (!strncmp(&Buffer[0],"RunTimes",LabelLength)){
      printf("RunTimes:");

      sscanf(Buffer,"RunTimes %d %d %d",&(t.tmax1),&(t.tmax2),&(t.tmax3));

      printf(" %f %f %f (in fs)\n",(t.tmax1-1)*(deltat),(t.tmax2)*(deltat),(t.tmax3-1)*(deltat));
      continue;
    }

    // TimeIncrement keyword
    /*    if (!strncmp(&Buffer[0],"TimeIncrement",LabelLength)){
      printf("TimeIncrement:");

      sscanf(Buffer,"TimeIncrement %d %d %d",&(t.dt[1]),&(t.dt[2]),&(t.dt[3]));

      printf(" %f %f %f (in fs)\n",(t.dt[1])*(deltat),(t.dt[2])*(deltat),(t.dt[3])*(deltat));
      continue;
      }*/
    t.dt[1]=t.dt[2]=t.dt[3]=1;

    // MinTimes keyword
    /*    if (!strncmp(&Buffer[0],"MinTimes",LabelLength)){
      printf("MinTimes:");

      sscanf(Buffer,"MinTimes %d %d %d",&(t.tmin1),&(t.tmin2),&(t.tmin3));

      printf(" %f %f %f (in fs)\n",t.tmin1*(deltat),t.tmin2*(deltat),t.tmin3*(deltat));
      continue;
      }*/
    // MinFrequencies keyword
    if (!strncmp(&Buffer[0],"MinFrequencies",LabelLength)){
      printf("MinFrequencies:");

      sscanf(Buffer,"MinFrequencies %f %f %f",&(w.min1),&(w.min2),&(w.min3));
      printf(" w1: %f  w2: %f  w3: %f\n",w.min1,w.min2,w.min3);
      continue;
    }
    // MaxFrequencies keyword
    if (!strncmp(&Buffer[0],"MaxFrequencies",LabelLength)){
      printf("MaxFrequencies:");

      sscanf(Buffer,"MaxFrequencies %f %f %f",&(w.max1),&(w.max2),&(w.max3));
      printf(" w1: %f  w2: %f  w3: %f\n",w.max1,w.max2,w.max3);
      continue;
    }
    
    
  } while (1==1);
  
  if ((w.min1<0) || (w.min2<0) || (w.min3<0) || (w.max1<0) || (w.max2<0) || (w.max3<0)){
    printf("All field frequencies must be positive!\n");
    exit(0);
  } 
  if ((w.min1>w.max1) || (w.min2>w.max2) || (w.min3>w.max3)){
    printf("The max frequency must be larger than the min frequency!\n");
    exit(0);
  }

  /*
  if (tv1==0 || tv2==0){
    printf("Warning!\n");
    printf("The Timevariable keyword was not specified!\n");
    tv1=1,tv2=3;
    printf("Using default. T2 - waiting time\n");
    tv3=6-tv1-tv2;
    }*/

  shift1=0.5*(w.min1+w.max1),shift3=-0.5*(w.min3+w.max3);

  printf("Initializing FFTW\n");
  if (fft==0){
    printf("Fatal error.\n");
    printf("FFT keyword not specified!\n");
    exit(0);
  }
  // Prepare 2DFFT
  fftIn = fftw_malloc(sizeof(fftw_complex) * (fft*fft));
  fftOut = fftw_malloc(sizeof(fftw_complex) * (fft*fft));
  fftPlan = fftw_plan_dft_2d(fft,fft,fftIn,fftOut,FFTW_FORWARD,FFTW_ESTIMATE);
  printf("-----------------\n");

  // Loop over response functions kI
  for (pol=0;pol<3;pol++){
    if (pol==0) sprintf(timeFName,"RparI.dat"),sprintf(frequencyFName,"Rwpar.I.dat");
    if (pol==1) sprintf(timeFName,"RperI.dat"),sprintf(frequencyFName,"Rwper.I.dat");
    if (pol==2) sprintf(timeFName,"RcroI.dat"),sprintf(frequencyFName,"Rwcro.I.dat");
    
    // Clear array
    for (i=0;i<fft;i++){
      for (j=0;j<fft;j++){
	index=j+fft*i;
	fftIn[index][0]=0;
	fftIn[index][1]=0;
      }
    }

    // Read response function
    input=fopen(timeFName,"r");
    if (input==NULL){
      printf("The response function file %s was not found!\n",timeFName);
//      exit(0); Changed April 2021 to simply skip to the next polarization
//      if the response function is missing for this one. TLC
      printf("Warning: Skipping this response function polarization!\n"); 
	continue;
    }
    
    printf("Reading response function (kI) file!\n");
    for (i=0;i<t.tmax1;i+=t.dt[1]){
      j=t.tmax2;
      for (k=0;k<t.tmax3;k+=t.dt[3]){
	fscanf(input,"%f %f %f %e %e",&ti[1],&ti[2],&ti[3],&rr,&ir);
	// Test is the fixed time is right
	if (round(ti[tv3]/deltat)==round(fixedtime/deltat) || fixedtime<-.5){
	  index=round(ti[tv2]/(deltat*t.dt[tv2])+fft*ti[tv1]/(deltat*t.dt[tv1]));
	  /* Apply appodization */
	  if (homo>0.0){
	     rr=rr*exp(-(ti[1]+ti[3])/2/homo);
             ir=ir*exp(-(ti[1]+ti[3])/2/homo);
	  }
          if (inhomo>0.0){
             rr=rr*exp(-(ti[1]-ti[3])*(ti[1]-ti[3])/2/inhomo/inhomo);
             ir=ir*exp(-(ti[1]-ti[3])*(ti[1]-ti[3])/2/inhomo/inhomo);
          }
	  if (index<fft*fft){
	    fftIn[index][0]=rr;
	    fftIn[index][1]=ir;
	  }
	}
      }
    }
    fclose(input);
    printf("Response function (kI) file read!\n");
    
    // Do fft
    fftw_execute(fftPlan);
    
    // Fix units
    
    // Write response function
    output=fopen(frequencyFName,"w");

    for (i=fft/2;i<fft;i++){
      for (j=fft/2;j<fft;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=(-fft+i)/deltat/c_v/fft-shift1;
	if (i<fft/2) w1=i/deltat/c_v/fft-shift1;
	if (j<fft/2) w3=j/deltat/c_v/fft-shift3;
	if (j>=fft/2) w3=(-fft+j)/deltat/c_v/fft-shift3;
	//	printf("W1 %f W3 %f %f %f\n",w1,w3,shift1,shift3);
	if (w1>-w.max1 && w3<w.max3){
	  if (w1<-w.min1 && w3>w.min3){
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }
      //      if (form==2) fprintf(output,"\n");
      for (j=0;j<fft/2;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=(-fft+i)/deltat/c_v/fft-shift1;
	if (i<fft/2) w1=i/deltat/c_v/fft-shift1;
	if (j<fft/2) w3=j/deltat/c_v/fft-shift3;
	if (j>=fft/2) w3=(-fft+j)/deltat/c_v/fft-shift3;
	if (w1>-w.max1 && w3<w.max3){
	  if (w1<-w.min1 && w3>w.min3){
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }
      if ((form==1 || form==2) && w1>-w.max1 && w1<-w.min1){
	fprintf(output,"\n");
      }

    }
    for (i=0;i<fft/2;i++){
      for (j=fft/2;j<fft;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=(-fft+i)/deltat/c_v/fft-shift1;
	if (i<fft/2) w1=i/deltat/c_v/fft-shift1;
	if (j<fft/2) w3=j/deltat/c_v/fft-shift3;
	if (j>=fft/2) w3=(-fft+j)/deltat/c_v/fft-shift3;
	if (w1>-w.max1 && w3<w.max3){
	  if (w1<-w.min1 && w3>w.min3){
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }
      //      if (form==1) fprintf(output,"\n");
      for (j=0;j<fft/2;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=(-fft+i)/deltat/c_v/fft-shift1;
	if (i<fft/2) w1=i/deltat/c_v/fft-shift1;
	if (j<fft/2) w3=j/deltat/c_v/fft-shift3;
	if (j>=fft/2) w3=(-fft+j)/deltat/c_v/fft-shift3;
	if (w1>-w.max1 && w3<w.max3){
	  if (w1<-w.min1 && w3>w.min3){
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }
      if ((form==1 || form==2) && w1>-w.max1 && w1<-w.min1) fprintf(output,"\n");
    }
    fclose(output);
  }



  // Loop over response functions kII
  shift3=-shift3;
  for (pol=0;pol<3;pol++){
    if (pol==0) sprintf(timeFName,"RparII.dat"),sprintf(frequencyFName,"Rwpar.II.dat");
    if (pol==1) sprintf(timeFName,"RperII.dat"),sprintf(frequencyFName,"Rwper.II.dat");
    if (pol==2) sprintf(timeFName,"RcroII.dat"),sprintf(frequencyFName,"Rwcro.II.dat");

    // Clear array
    for (i=0;i<fft;i++){
      for (j=0;j<fft;j++){
	index=j+fft*i;
	fftIn[index][0]=0;
	fftIn[index][1]=0;
      }
    }
    
    // Read response function kII
    input=fopen(timeFName,"r");
    if (input==NULL){
      printf("The response function file %s was not found!\n",timeFName);
//      exit(0); April 2021 Skips if the polarization is missing.
      printf("Warning: Skipping this response function polarization!\n");
      continue; 
    }
    
    printf("Reading response function file (kII)!\n");
    // Loop starts at zero, bug fixed 16/3-2012 TLC (found by LW)
    for (i=0;i<t.tmax1;i+=t.dt[1]){
      j=t.tmax2;
      for (k=0;k<t.tmax3;k+=t.dt[3]){
	fscanf(input,"%f %f %f %e %e",&ti[1],&ti[2],&ti[3],&rr,&ir);
	// Test is the fixed time is right
	if (round(ti[tv3]/deltat)==round(fixedtime/deltat) || fixedtime<-.5){
	  if (ti[tv1]==0){
	    index=round((ti[tv2])/(deltat*t.dt[tv2])+fft*(ti[tv1]/(deltat*t.dt[tv1])));
	  } else {
	    index=round((ti[tv2])/(deltat*t.dt[tv2])+fft*((ti[tv1])/(deltat*t.dt[tv1])));
	  }
	  /* Apply appodization */
          if (homo>0.0){
             rr=rr*exp(-(ti[1]+ti[3])/2/homo);
             ir=ir*exp(-(ti[1]+ti[3])/2/homo);
          }
          if (inhomo>0.0){
             rr=rr*exp(-(ti[1]+ti[3])*(ti[1]+ti[3])/2/inhomo/inhomo);
             ir=ir*exp(-(ti[1]+ti[3])*(ti[1]+ti[3])/2/inhomo/inhomo);
          }
	  if (index<fft*fft){
	    fftIn[index][0]+=rr;
	    fftIn[index][1]+=ir;
	  }	    
	}
      }
    }
    fclose(input);
    printf("Response function kII file read!\n"); 
    
    // Do fft
    fftw_execute(fftPlan);
    //      printf("%f %f\n",1/deltat/c_v/fft,1/deltat/c_v);
    sw=1/deltat/c_v;
    // Fix units
    
    rwma1=shift1,rwmi1=shift1,rwma3=shift3,rwmi3=shift3;
    // Write response function
    output=fopen(frequencyFName,"w");
    if (form==1) axisFH=fopen("waxis.dat","w");
    for (i=fft/2;i<fft;i++){
      for (j=fft/2;j<fft;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=i/deltat/c_v/fft+shift1-sw;
	if (i<fft/2) w1=(fft+i)/deltat/c_v/fft+shift1-sw;
	if (j<fft/2) w3=(fft+j)/deltat/c_v/fft+shift3-sw;
	if (j>=fft/2) w3=j/deltat/c_v/fft+shift3-sw;
	
	if (w1<w.max1 && w3<w.max3){
	  if (w1>w.min1 && w3>w.min3){
	    if (w1<rwmi1) rwmi1=w1;
	    if (w1>rwma1) rwma1=w1;
	    if (w3<rwmi3) rwmi3=w3;
	    if (w3>rwma3) rwma3=w3;
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }

      //      if (form==1 && w1<w.max1 && w1>w.min1) fprintf(output,"\n");
      for (j=0;j<fft/2;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=i/deltat/c_v/fft+shift1-sw;
	if (i<fft/2) w1=(fft+i)/deltat/c_v/fft+shift1-sw;
	if (j<fft/2) w3=(fft+j)/deltat/c_v/fft+shift3-sw;
	if (j>=fft/2) w3=j/deltat/c_v/fft+shift3-sw;
	if (w1<w.max1 && w3<w.max3){
	  if (w1>w.min1 && w3>w.min3){
	    if (w1<rwmi1) rwmi1=w1;
	    if (w1>rwma1) rwma1=w1;
	    if (w3<rwmi3) rwmi3=w3;
	    if (w3>rwma3) rwma3=w3;
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }

      if ((form==1 || form==2) && w1<w.max1 && w1>w.min1){
	fprintf(output,"\n");
	if (form==1) fprintf(axisFH,"%f\n",w1);
      }
    }
    for (i=0;i<fft/2;i++){
      for (j=fft/2;j<fft;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=i/deltat/c_v/fft+shift1-sw;
	if (i<fft/2) w1=(fft+i)/deltat/c_v/fft+shift1-sw;
	if (j<fft/2) w3=(fft+j)/deltat/c_v/fft+shift3-sw;
	if (j>=fft/2) w3=j/deltat/c_v/fft+shift3-sw;
	if (w1<w.max1 && w3<w.max3){
	  if (w1>w.min1 && w3>w.min3){
	    if (w1<rwmi1) rwmi1=w1;
	    if (w1>rwma1) rwma1=w1;
	    if (w3<rwmi3) rwmi3=w3;
	    if (w3>rwma3) rwma3=w3;
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }

      //      if (form==1 && w1<w.max1 && w1>w.min1) fprintf(output,"\n");
      for (j=0;j<fft/2;j++){
	index=j+fft*i;
	if (i>=fft/2) w1=i/deltat/c_v/fft+shift1-sw;
	if (i<fft/2) w1=(fft+i)/deltat/c_v/fft+shift1-sw;
	if (j<fft/2) w3=(fft+j)/deltat/c_v/fft+shift3-sw;
	if (j>=fft/2) w3=j/deltat/c_v/fft+shift3-sw;
	if (w1<w.max1 && w3<w.max3){
	  if (w1>w.min1 && w3>w.min3){
	    if (w1<rwmi1) rwmi1=w1;
	    if (w1>rwma1) rwma1=w1;
	    if (w3<rwmi3) rwmi3=w3;
	    if (w3>rwma3) rwma3=w3;
	    if (form==0 || form==2) fprintf(output,"%f %f %e %e\n",w1,w3,
				 fftOut[index][0],fftOut[index][1]);
	    if (form==1) fprintf(output,"%e ",fftOut[index][1]);
	  }
	}
      }

      if ((form==1 || form==2) && w1<w.max1 && w1>w.min1){
	fprintf(output,"\n");
	if (form==1) fprintf(axisFH,"%f\n",w1);
      }
    }
    //    printf("y\n");
    if (form==1) fclose(axisFH);
    fclose(output);
    //    printf("z\n");
  }

  // Exit if printing in matlab format (user should add S1 and S2 in matlab);
  if (form==1) exit(0);

  // Add S1 and S2 for dislin use
  dw1=1/deltat/c_v/fft;
  dw3=1/deltat/c_v/fft;
  //  rwma1=1/deltat/c_v/2-shift1;
  //  rwmi1=-1/deltat/c_v/2-shift1;
  //  rwma3=1/deltat/c_v/2-shift3;
  //  rwmi3=-1/deltat/c_v/2-shift3;
  // Find dimensions
  //  p1I=p1II=p2D1=round((w.max1-w.min1)/dw1-0.5);
  //  p3I=p3II=p2D3=round((w.max3-w.min3)/dw3-0.5);
  p1I=p1II=p2D1=round((rwma1-rwmi1)/dw1+1);
  p3I=p3II=p2D3=round((rwma3-rwmi3)/dw3+1);
  //  printf("A1 %d A2 %d %f %f %f %f %f %f\n",p1I,p3I,(w.max1-w.min1)/dw1,(w.max3-w.min3)/dw3,rwma1,rwmi1,rwma3,rwmi3);
  //  p2D1--,p2D3--;
  printf("Adding rephasing and nonrephasing contributions!\n");

  // Define arrays
  kI=matrix(1,p1I,1,p3I);
  kI1=matrix(1,p1I,1,p3I);
  kI3=matrix(1,p1I,1,p3I);
  kII=matrix(1,p1II,1,p3II);
  kII1=matrix(1,p1II,1,p3II);
  kII3=matrix(1,p1II,1,p3II);
  IR=matrix(1,p2D1,1,p2D3);
  IR1=matrix(1,p2D1,1,p2D3);
  IR3=matrix(1,p2D1,1,p2D3);
  kIi=matrix(1,p1I,1,p3I);
  kIIi=matrix(1,p1II,1,p3II);
  IRi=matrix(1,p2D1,1,p2D3);

  for (pol=0;pol<3;pol++){
    if (pol==0) sprintf(kIFName,"Rwpar.I.dat"),sprintf(kIIFName,"Rwpar.II.dat"),sprintf(twoDFName,"2D.par.dat"),sprintf(pPFName,"PP.par.dat");
    if (pol==1) sprintf(kIFName,"Rwper.I.dat"),sprintf(kIIFName,"Rwper.II.dat"),sprintf(twoDFName,"2D.per.dat"),sprintf(pPFName,"PP.per.dat");
    if (pol==2) sprintf(kIFName,"Rwcro.I.dat"),sprintf(kIIFName,"Rwcro.II.dat"),sprintf(twoDFName,"2D.cro.dat"),sprintf(pPFName,"PP.cro.dat");											
    // Read files
    // Read response function kI
    input=fopen(kIFName,"r");
    if (input==NULL){
      printf("The response function file %s was not found!\n",kIFName);
      //exit(0);
      printf("Warning: Skipping this response function polarization!\n");
      continue;
    }
    printf("Reading response function file (kI)!\n");
    for (i=1;i<=p1I;i++){
      for (j=1;j<=p3I;j++){
	fscanf(input,"%f %f %e %e",&kI1[i][j],&kI3[i][j],&kIi[i][j],&kI[i][j]);
	kI1[i][j]=-kI1[i][j];
      }
    }
    fclose(input);
    
    // Read response function kII
    input=fopen(kIIFName,"r");
    if (input==NULL){
      printf("The response function file %s was not found!\n",kIIFName);
      //exit(0);
      printf("Warning: Skipping this response function polarization!\n");
    }
    printf("Reading response function file (kII)!\n");
    for (i=1;i<=p1II;i++){
      for (j=1;j<=p3II;j++){
	fscanf(input,"%f %f %e %e",&kII1[i][j],&kII3[i][j],&kIIi[i][j],&kII[i][j]);
      }
    }
    fclose(input);
    
    // Find grid frequencies
    // Loop over grid points
    for (k=1;k<=p2D1;k++){
      for (l=1;l<=p2D3;l++){
	//	IR1[k][l]=w.min1+(w.max1-w.min1)*k/(p2D1);
	//	IR3[k][l]=w.min3+(w.max3-w.min3)*l/(p2D3);
	IR1[k][l]=rwmi1+(rwma1-rwmi1)*k/(p2D1);
	IR3[k][l]=rwmi3+(rwma3-rwmi3)*l/(p2D3);
      }
    }
    
    // Add contributions from kI
    for (i=1;i<=p1I;i++){
      for (j=1;j<=p3I;j++){
	IR[i][j]=kI[p1I-i+1][j]+kII[i][j];
	IRi[i][j]=kIi[p1I-i+1][j]+kIIi[i][j];
      }
    }
    
    // Save 2DIR spectrum
    output=fopen(twoDFName,"w");
    // Loop over grid points
    for (k=1;k<=p2D1;k++){
      for (l=1;l<=p2D3;l++){
	fprintf(output,"%f %f %e %e\n",kII1[k][l],kII3[k][l],IRi[k][l],IR[k][l]);
      }
      if (form==2) fprintf(output,"\n");
    }
    fclose(output);
    
    // Calculate broad pump - pump probe spectrum
    
    output=fopen(pPFName,"w");
    
    if (output==NULL){
      printf("Name for broadband pump probe not given.\n");
      printf("Skipping pump probe calculation.\n");
    } else {
      fprintf(output,"### Broadband pump probe spectrum\n");  
      // Calculate pump probe spectrum
      for (l=1;l<=p2D3;l++){
	pumpProbe=0;
	for (k=1;k<=p2D1;k++){
	  pumpProbe+=IR[k][l];
	}
	fprintf(output,"%f %e\n",IR3[1][l],pumpProbe);
	//    fprintf(output,"%f %f\n",IR3[k][l],pumpProbe);
      }
      fclose(output);
    }
  }
}
