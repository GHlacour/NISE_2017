#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "types.h"
#include "types_MPI.h"
#include "NISE_subs.h"
#include "readinput.h"
#include "NISE.h"
#include "absorption.h"
//#include "c_absorption.h"
#include "calc_DOS.h"
#include "luminescence.h"
#include "raman.h"
#include "sfg.h"
#include "calc_2DIR.h"
#include "calc_2DIRraman.h"
#include "calc_2DES.h"
#include "calc_CG_2DES.h"
#include "eq_den.h"
#include "analyse.h"
#include "calc_CD.h"
#include "calc_LD.h"
#include "calc_Diffusion.h"
#include "population.h"
#include "anisotropy.h"
#include "mcfret.h"
#include "propagate.h"
#include "correlate.h"
#include <mpi.h>
#include "omp.h"
#include "1DFFT.h"




/* This is the 2017 version of the NISE program
   It allow calculating linear absorption and 2D(IR) spectra
   This version of the program allow to utilize sparce matrices
   and parallel code

   The code was programmed by Thomas la Cour Jansen, RuG, 2017
*/

// The main routine
int main(int argc, char* argv[]) {
    // Initialize MPI
    int parentRank, subRank, parentSize, subSize;
    int thread_id,cpus;
    MPI_Comm subComm, rootComm;

    int threadingProvided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadingProvided);
    MPI_Comm_rank(MPI_COMM_WORLD, &parentRank);
    MPI_Comm_size(MPI_COMM_WORLD, &parentSize);

    // We require funneled multithreading so if that is not allowed we gracefully fail
    if(threadingProvided < MPI_THREAD_FUNNELED) {
        printf("Error: MPI library does not support threading!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // We split up the processing in smaller chunks, each set of MPI processes will make shared memory
    // for the global state. Then only the master processes within each chunk will communicate among
    // each other, to the main master that will also do all logfile printing and reductions.

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &subComm);
    MPI_Comm_rank(subComm, &subRank);
    MPI_Comm_size(subComm, &subSize);

    // We now make another communicator that contains only all roots of the newly created subComms
    MPI_Comm_split(MPI_COMM_WORLD, subRank == 0 ? 0 : MPI_UNDEFINED, 0, &rootComm);

    // Define data structures to transmit via MPI
    MPI_Datatype t_non_type;
    MPI_Type_create_struct(T_NON_TYPE.length, T_NON_TYPE.blocklengths, T_NON_TYPE.offsets, T_NON_TYPE.types, &t_non_type);
    MPI_Type_commit(&t_non_type);

    // Define variables
    t_non* non = calloc(1, sizeof(t_non));
    time_t timeStart;

    // Do master initialization work
    int initResult = 0;
    if (parentRank == 0) {
        /* Time parameters */
        time(&timeStart);

        /* Intro */
        printf("----- ----- ----- ----- ----- -----\n");
        printf("  Running the 23/8-2017 version of\n");
        printf("              NISE3.1\n");
        printf("    by Thomas la Cour Jansen.\n");
        printf("----- ----- ----- ----- ----- -----\n");
        printf("\n");

        // Read the input
        readInput(argc, argv, non);

        // Do initial check of the configuration
        initResult = control(non);
        initResult = autodetect_singles(non);

        // Create new log file
        FILE* logFile = fopen("NISE.log", "w");
        if (logFile == NULL) {
            printf("Could not open log file! Disk full?\n");
            initResult = 1;
        }
        fprintf(logFile, "Log\n");
        fclose(logFile);
    }

    // Sync the initresult to other processes: if we failed we should terminate
    MPI_Bcast(&initResult, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (initResult != 0) {
        MPI_Finalize();
        exit(initResult);
    }

    // Synchronize the configuration
    MPI_Bcast(non, 1, t_non_type, 0, MPI_COMM_WORLD);

    // Synchronize p_sites array, allocate it first if not root process, if necessary
    int has_psites = (parentRank == 0 && non->psites != NULL);
    MPI_Bcast(&has_psites, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (has_psites) {
        if (parentRank != 0) {
            non->psites = calloc(non->singles, sizeof(int));
        }

        MPI_Bcast(non->psites, non->singles, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Inform user of parallel info */
    if (parentRank==0){
       #pragma omp parallel private(thread_id) 
       {
           thread_id=omp_get_thread_num();
           if (thread_id==0){
	       printf("\n=== Parallel computing information ===");	   
               printf("\nDetected %d openMP threads ",omp_get_num_threads());
               printf("and %d MPI instances.\n\n",parentSize);
	       cpus=omp_get_num_threads()*parentSize;
           }
       }
    }

    // Delegate to different subroutines depending on the technique

    // Call the Hamiltonian Analysis routine
    if (string_in_array(non->technique,(char*[]){"Analyse","Analyze"},2)){
        // Does not support MPI
        if (parentRank == 0)
            analyse(non);
    }

    // Call the Hamiltonian Correlate routine
    if (!strcmp(non->technique, "Correlation")) {
        // Does not support MPI
        if (parentRank == 0){
	    if (cpus>1) not_parallel();
            calc_Correlation(non);
	}
    }

    // Call the Population Transfer routine
    if (string_in_array(non->technique,(char*[]){"Pop","Population"},2)){
        // Does not support MPI
        if (parentRank == 0)
            population(non);
    }

    // Call the Exciton Diffusion routine
    if (string_in_array(non->technique,(char*[]){"Dif","Diffusion"},2)){
        // Does not support MPI
        if (parentRank == 0)
            calc_Diffusion(non);
    }

    // Call the Anisotropy and Rotational Correlation routine
    if (string_in_array(non->technique,(char*[]){"Ani","Anisotropy"},2)){
        // Does not support MPI
        if (parentRank == 0)
            anisotropy(non);
    }

    // Call the Linear Absorption Routine
    if (!strcmp(non->technique, "Absorption")) {
        // Does not support MPI
        if (parentRank == 0) {
		if (cpus>1) not_parallel();
                absorption(non);
        }
    }

    /* Call the Linear DOS Routine */
    if (!strcmp(non->technique, "DOS")) {
        /* Does not support MPI */
        if (parentRank == 0) {
                calc_DOS(non);
        }
    }

    /* Call the MCFRET Routine */
        if (string_in_array(non->technique,(char*[]){"MCFRET",
	   "MCFRET-Autodetect","MCFRET-Absorption","MCFRET-Emission",
	   "MCFRET-Coupling","MCFRET-Rate","MCFRET-Analyse",
	   "MCFRET-Density"},8)){
        /* Does not support MPI */
        if (parentRank == 0) {
                mcfret(non);
        }
    }

    // Call the Luminescence Routine
    if (string_in_array(non->technique,(char*[]){"Luminescence","PL","Fluorescence"},3)){
        // Does not support MPI
        if (parentRank == 0)
            luminescence(non);
    }

    // Call the Linear Dichroism Routine
    if (!strcmp(non->technique, "LD")) {
        // Does not support MPI
        if (parentRank == 0){
            if (cpus>1) not_parallel();
            LD(non);
	}
    }

    // Call the Circular Dichroism Routine
    if (!strcmp(non->technique, "CD")) {
        // Does not support MPI
        if (parentRank == 0)
            calc_CD(non);
    }

    /* Call the Raman Routine */
    if (!strcmp(non->technique, "Raman")) {
        /* Does not support MPI */
        if (parentRank == 0){
            if (cpus>1) not_parallel();
            raman(non);
	      }
     }

    /* Call the Sum Frequency Generation Routine */
    if (!strcmp(non->technique, "SFG")) {
        //Does not support MPI
        if (parentRank == 0){
            if (cpus>1) not_parallel();
            sfg(non);
        }
    }

    // Call the 2DIR calculation routine

    if (string_in_array(non->technique,(char*[]){"2DIR","GBIR","SEIR","EAIR","noEAIR"},5)){
	    // Does support MPI
        calc_2DIR(non,parentRank, parentSize, subRank, subSize, subComm, rootComm);
    }

    // Call the 2DIRraman calculation routine
    if (string_in_array(non->technique,(char*[]){"2DIRraman","2DIRraman1","2DIRraman2","2DIRraman3","2DIRramanI","2DIRramanII"},6)){
        // Does support MPI
        calc_2DIRraman(non,parentRank, parentSize, subRank, subSize, subComm, rootComm);
    }

    // Call the 2DSFG calculation routine
    if (string_in_array(non->technique,(char*[]){"2DSFG","GBSFG","SESFG","EASFG","noEASFG"},5)){
	printf("2DSFG is not yet implemented in NISE2017!");
	exit(0);
    }

    // Call the 2DUVvis calculation routine
    if (string_in_array(non->technique,(char*[]){"2DUVvis","GBUVvis","SEUVvis","EAUVvis","noEAUVvis"},5)){
        // Does support MPI
        calc_2DES(non,parentRank, parentSize, subRank, subSize, subComm, rootComm);
    }

    /* Call the CG_2DES Routine */
    if (!strcmp(non->technique, "CG_2DES") ||  (!strcmp(non->technique, "CG_2DES_doorway")) || 
     (!strcmp(non->technique, "CG_2DES_P_DA")) ||  (!strcmp(non->technique, "CG_2DES_window_GB"))
     ||  (!strcmp(non->technique, "CG_2DES_window_SE")) ||  (!strcmp(non->technique, "CG_2DES_window_EA"))
     ||  (!strcmp(non->technique, "CG_full_2DES_segments")) ||  (!strcmp(non->technique, "combine_CG_2DES"))
     ||  (!strcmp(non->technique, "CG_2DES_waitingtime")) ) {
        /* Does not support MPI */
        if (parentRank == 0)
            calc_CG_2DES(non);
    }

    // Call the 2DFD calculation routine
    if (!strcmp(non->technique, "2DFD")) { }

    // Call the 1DFT calculation routine
    if (!strcmp(non->technique, "1DFFT")) {
        // Does not support MPI
        if (parentRank == 0) {
		if (cpus>1) not_parallel();
                ONE_DFFT(non);
        }
    }
    // Call the lineshape funnction for absorption
    if (!strcmp(non->technique, "Lineshape_FFT")) {
        // Does not support MPI
        if (parentRank == 0) {
		if (cpus>1) not_parallel();
                Lineshape_FFT(non);
        }
    }




    // Do Master wrap-up work
    if (parentRank == 0) {
        // The calculation is finished, lets write end log
        time_t timeEnd;
        time(&timeEnd);
        char* timeText = time_diff(timeStart, timeEnd);
        log_item("Finished response!\nWriting final timing to log file!\n%s", timeText);
        free(timeText);

        printf("------------------------------------\n");
        printf(" NISE program succesfully completed\n");
        printf("------------------------------------\n\n");
    }

    // Clean up
    free(non->psites);
    free(non);

    MPI_Comm_free(&subComm);
    if(rootComm != MPI_COMM_NULL) MPI_Comm_free(&rootComm);
    MPI_Type_free(&t_non_type);
    MPI_Finalize();

    return 0;
}
