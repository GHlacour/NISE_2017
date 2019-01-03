#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "polar.h"
float static ithirty=1.0/30;

/* Polarization direction averaging */
/* Retruns the weight of the molecular polarization direction 'x' */
/* for the lab frame polarization direction pol */
float polarweight(int pol,int x){
  float weight;
  /* Parallel polarization: ZZZZ */
  if (pol==0){
    if (x<=2) weight=6*ithirty; //6.0/30;
    if (x>=3) weight=2*ithirty; //2.0/30;
  }
  /* Perpendicular polarization: ZZYY */
  if (pol==1){
    if (x<=2) weight=2*ithirty; //2.0/30;
    if (x>=3 && x<=8 ) weight=4*ithirty; //4.0/30;
    if (x>=9) weight=-1*ithirty; //-1.0/30;
  }
  /* Cross polarization -1/2 (ZYYZ-YZYZ) */
  if (pol==2){
    if (x<=2) weight=0.0; // XXXX
    if (x>=3 && x<=8 ) weight=0.0; // XXYY
    if (x>=9 && x<=14) weight=2.5*ithirty; //2.5/30; // XYXY
    if (x>=15) weight=-2.5*ithirty; //-2.5/30; // XYYX
  }
  return weight;
}

/* Return directions */
void polar(int xp[],int x){
  /* XXXX moleculer frame dipoles */
  if (x==0){
    xp[0]=0,xp[1]=0,xp[2]=0,xp[3]=0;
  }
  if (x==1){
    xp[0]=1,xp[1]=1,xp[2]=1,xp[3]=1;
  }
  if (x==2){
    xp[0]=2,xp[1]=2,xp[2]=2,xp[3]=2;
  }
  /* XXYY molecular frame dipoles */
  if (x==3){
    xp[0]=0,xp[1]=0,xp[2]=1,xp[3]=1;
  }
  if (x==4){
    xp[0]=0,xp[1]=0,xp[2]=2,xp[3]=2;
  }
  if (x==5){
    xp[0]=1,xp[1]=1,xp[2]=0,xp[3]=0;
  }
  if (x==6){
    xp[0]=1,xp[1]=1,xp[2]=2,xp[3]=2;
  }
  if (x==7){
    xp[0]=2,xp[1]=2,xp[2]=0,xp[3]=0;
  }
  if (x==8){
    xp[0]=2,xp[1]=2,xp[2]=1,xp[3]=1;
  }
  /* XYXY molecular frame dipoles */
  if (x==9){
    xp[0]=0,xp[1]=1,xp[2]=0,xp[3]=1;
  }
  if (x==10){
    xp[0]=0,xp[1]=2,xp[2]=0,xp[3]=2;
  }
  if (x==11){
    xp[0]=1,xp[1]=0,xp[2]=1,xp[3]=0;
  }
  if (x==12){
    xp[0]=1,xp[1]=2,xp[2]=1,xp[3]=2;
  }
  if (x==13){
    xp[0]=2,xp[1]=0,xp[2]=2,xp[3]=0;
  }
  if (x==14){
    xp[0]=2,xp[1]=1,xp[2]=2,xp[3]=1;
  }
  /* XYYX molecular frame dipoles */
  if (x==15){
    xp[0]=0,xp[1]=1,xp[2]=1,xp[3]=0;
  }
  if (x==16){
    xp[0]=0,xp[1]=2,xp[2]=2,xp[3]=0;
  }
  if (x==17){
    xp[0]=1,xp[1]=0,xp[2]=0,xp[3]=1;
  }
  if (x==18){
    xp[0]=1,xp[1]=2,xp[2]=2,xp[3]=1;
  }
  if (x==19){
    xp[0]=2,xp[1]=0,xp[2]=0,xp[3]=2;
  }
  if (x==20){
    xp[0]=2,xp[1]=1,xp[2]=1,xp[3]=2;
  }
  return;
}


