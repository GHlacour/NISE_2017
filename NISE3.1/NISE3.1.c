#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "omp.h"
#include "types3.1.h"
#include "NISE3.1subs.h"
#include "readinput.h"
#include "NISE3.1.h"
#include "absorption.h"
#include "c_absorption.h"
#include "luminescence.h"
#include "calc_2DIR.h"
#include "calc_2DES.h"
#include "analyse.h"
#include "calc_CD.h"

/* This is the 2017 version of the NISE program
   It allow calculating linear absorption and 2D(IR) spectra
   This version of the program allow to utilize sparce matrices
   and parallel code

   The code was programmed by Thomas la Cour Jansen, RuG, 2017
*/

// The main routine
int main(int argc, char *argv[])
{
  /* Define log file */
  FILE *log;

  /* Time parameters */
  time_t time_now,time_old,time_0;
  /* Initialize time */
  time(&time_now);
  time(&time_0);

  /* Define control structure */
  t_non *non;

  /* Intro */
  printf("----- ----- ----- ----- ----- -----\n");
  printf("  Running the 23/8-2017 version of\n");
  printf("              NISE3.1\n");
  printf("    by Thomas la Cour Jansen.\n");
  printf("----- ----- ----- ----- ----- -----\n");
  printf("\n");

  // Allocate memory for input information
  non=(t_non *)calloc(1,sizeof(t_non));

  // Read the input
  readInput(argc,argv,non);
  // Do initial control
  control(non);

  // Delegate to different subroutines depending on the technique

  // Call the Hamiltonian Analysis routine
  if(!strcmp(non->technique,"Analyse")){
    analyse(non);
  }

  // Call the Population Transfer routine
  if(!strcmp(non->technique,"Pop")){
  }

  // Call the Exciton Diffusion routine
  if(!strcmp(non->technique,"Dif")){
  }

  // Call the Anisotropy and Rotational Correlation routine
  if(!strcmp(non->technique,"Ani")){
  }

  // Call the Linear Absorption Routine
  if(!strcmp(non->technique,"Absorption")){
    if (!strcmp(non->hamiltonian,"Coupling")){
      c_absorption(non);
    } else {
      absorption(non);
    }
  }

  // Call the Luminescence Routine
  if(!strcmp(non->technique,"Luminescence")){
    luminescence(non);
  }

  // Call the Linear Dichroism Routine
  if(!strcmp(non->technique,"LD")){
  }

  // Call the Circular Dichroism Routine
  if(!strcmp(non->technique,"CD")){
    calc_CD(non);
  }

  // Call the Raman Routine
  if(!strcmp(non->technique,"Raman")){
  }

  // Call the Sum Frequency Generation Routine
  if(!strcmp(non->technique,"SFG")){
  }

  // Call the 2DIR calculation routine
  if(!strcmp(non->technique,"2DIR")||(!strcmp(non->technique,"GB"))||(!strcmp(non->technique,"SE"))||(!strcmp(non->technique,"EA"))||(!strcmp(non->technique,"noEA"))){
    calc_2DIR(non);
  }

  // Call the 2DSFG calculation routine
  if(!strcmp(non->technique,"2DSFG")||(!strcmp(non->technique,"GBSFG"))||(!strcmp(non->technique,"SESFG"))||(!strcmp(non->technique,"EASFG"))||(!strcmp(non->technique,"noEASFG"))){
  }

  // Call the 2DUVvis calculation routine
  if(!strcmp(non->technique,"2DUVvis")||(!strcmp(non->technique,"GBUVvis"))||(!strcmp(non->technique,"SEUVvis"))||(!strcmp(non->technique,"EAUVvis"))||(!strcmp(non->technique,"noEAUVvis"))){
    calc_2DES(non);
  }

  // Call the 2DFD calculation routine
  if(!strcmp(non->technique,"2DFD")){
  }


  // The calculation is finished, lets write end log
  log=fopen("NISE.log","a");
  fprintf(log,"Finished Response!\n");
  fprintf(log,"Writing final timing to log file!\n");
  time_now=log_time(time_now,log);  
  fclose(log);

  free(non);

  printf("------------------------------------\n");
  printf(" NISE program succesfully completed\n");
  printf("------------------------------------\n\n");

  return 0;
}
