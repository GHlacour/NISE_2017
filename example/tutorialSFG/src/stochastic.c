#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "randomlib.h"

int main(int argc,char *argv[]){
  int length;
  float w1,w2,wa,wb,wao,wbo;
  float angle,time,sigma,angle2;
  float deltaw,J,deltat,w0;
  float a,b;
  float pi=3.14159265;
  int i;
  FILE *E_FH,*mu_FH,*alpha_FH;

  length=atoi(argv[1]);
  deltat=atof(argv[2]);
  sigma=atof(argv[3]);
  time=atof(argv[4]);
  angle=atof(argv[5]);
  deltaw=atof(argv[6]);
  w0=atof(argv[7]);
  J=atof(argv[8]);
  angle2=atof(argv[9]);

  a=exp(-deltat/time);
  b=sqrt(1-a*a);

  RandomInitialise(2511,1974);
  wao=RandomGaussian(0,sigma);
  wbo=RandomGaussian(0,sigma);

  E_FH=fopen("Energy.txt","w");
  mu_FH=fopen("Dipole.txt","w");
  alpha_FH=fopen("Alpha.txt","w");

  for(i=0;i<length;i++){
    wa=wao*a+b*RandomGaussian(0,sigma);
    wb=wbo*a+b*RandomGaussian(0,sigma);
    w1=wa;
    w2=wa*cos(angle*pi/180.0)+wb*sin(angle*pi/180.0);
    fprintf(E_FH,"%d %f %f %f\n",i,w1-deltaw/2.0+w0,J,w2+deltaw/2.0+w0);
    fprintf(mu_FH,"%d %f %f %f %f %f %f\n",i,0.0,0.0,0.0,sin(angle2*pi/180.0),1.0,cos(angle2*pi/180.0));
    fprintf(alpha_FH,"%d %f %f %f %f %f %f\n",i,1.0,1.0,1.0,1.0,10.0,10.0);
    wao=wa;
    wbo=wb;
  } 

  return 0;
}
