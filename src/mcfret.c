#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"
#include "propagate.h"
#include "propagate_segment.h"
#include "read_trajectory.h"
#include "mcfret4.h"

/* Main MCFRET routine only calling and combining the other subroutines */ 
void mcfret(t_non *non){
    int nn2;
    int segments;
    /* Response functions for emission and absorption: real and imaginary part*/
    float *re_Abs,*im_Abs;
    float *re_Emi,*im_Emi;
    float *re_e,*im_e; /* Eigenvalues */
    float *vl,*vr; /* Left and right eigenvectors */
    float *energy_cor; /* Effective correction energy for QC */
    float *J;
    float *E;
    float *rate_matrix;
    float *coherence_matrix;
    float *ave_vecr;

    /* Allocate memory for the response functions */
    nn2=non->singles*non->singles;
    re_Abs=(float *)calloc(nn2*non->tmax1,sizeof(float));
    im_Abs=(float *)calloc(nn2*non->tmax1,sizeof(float));
    re_Emi=(float *)calloc(nn2*non->tmax1,sizeof(float));
    im_Emi=(float *)calloc(nn2*non->tmax1,sizeof(float));
    J=(float *)calloc(nn2,sizeof(float));
    E=(float *)calloc(non->singles,sizeof(float));
    ave_vecr=(float *)calloc(non->singles*non->singles,sizeof(float));

    /* The rate matrix is determined by the integral over t1 for */
    /* Tr [ J * Abs(t1) * J * Emi(t1) ] */

    segments=project_dim(non);
    if (segments<2){
        printf(RED "Too few segments defined for MCFRET calculation!" RESET);
        exit(0);
    }
    rate_matrix=(float *)calloc(segments*segments,sizeof(float));
    coherence_matrix=(float *)calloc(segments*segments,sizeof(float));

    /* Tell the user that we are in the MCFRET Routine */
    if (string_in_array(non->technique,(char*[]){"MCFRET",
        "MCFRET-Autodetect","MCFRET-Absorption","ECFRET-Emission",
        "MCFRET-Coupling","MCFRET-Rate","MCFRET-Analyse",
        "MCFRET-Density","MCFRET-4th-approx", "MCFRET-4th-full" },10)){
        printf("Performing MCFRET calculation.\n");
    }

    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Density"))){
        /* Calculate the average density matrix */
	printf("Starting calculation of the average density matrix.\n");
        average_density_matrix(ave_vecr,non);
        write_matrix_to_file("Average_Density.dat",ave_vecr,non->singles);
    }

    /* Call the absorption routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Absorption"))){
	printf("Starting calculation of the MCFRET absorption matrix.\n");
        mcfret_propagation_segmented(re_Abs,im_Abs,non);
    }
   
    // /* Call the emission routine */
    // if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Emission"))){
    //     printf("Starting calculation of the MCFRET emission matrix.\n");
    // 	/* Read precalculated average density matrix */ 
    // 	if (!strcmp(non->technique, "MCFRET-Emission")){
    //        printf("Using precalculated average density matrix from file Average_Density.dat.\n");
    //	    read_matrix_from_file("Average_Density.dat",ave_vecr,non->singles);
    // 	}
    //    mcfret_response_function(re_Emi,im_Emi,non,1,ave_vecr);
    // }
    
    /* Call the coupling routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Coupling"))){
        printf("Starting calculation of the average inter segment coupling.\n");
        mcfret_coupling(J,non);
    }

    /* Call the rate routine routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Rate"))){
        printf("Starting calculation of the rate response function.\n");
        if ((!strcmp(non->technique, "MCFRET-Rate"))){
            /* Read in absorption, emission and coupling from file if needed */
	    printf("Calculating rate from precalculated absorption and coupling!\n");
	    read_matrix_from_file("CouplingMCFRET.dat",J,non->singles);
	    read_matrix_from_file("Average_Density.dat",ave_vecr,non->singles); 
	    read_response_from_file("TD_absorption_matrix.dat",re_Abs,im_Abs,non->singles,non->tmax1);
	    // read_response_from_file("TD_emission_matrix.dat",re_Emi,im_Emi,non->singles,non->tmax1);
            printf("Completed reading pre-calculated data.\n");
        }
        // mcfret_rate(rate_matrix,coherence_matrix,segments,re_Abs,im_Abs,re_Emi,im_Emi,J,non);
        mcfret_rate_from_abs(rate_matrix,coherence_matrix,segments,re_Abs,im_Abs,ave_vecr,J,non);

        /* Write the calculated ratematrix to file */
        write_matrix_to_file("RateMatrix.dat",rate_matrix,segments);
        /* Write the calculated coherence matrix to file */
        write_matrix_to_file("CoherenceMatrix.dat",coherence_matrix,segments);
    }

    /* Call the 4th order approximation routine */
    if (!strcmp(non->technique, "MCFRET-4th-approx")){
        printf("Starting calculation of the 4th order approximation: traces.\n");
	    printf("Calculating 4th order apprpximation from precalculated coupling and density\n");
	
	    read_matrix_from_file("CouplingMCFRET.dat",J,non->singles);
	    read_matrix_from_file("Average_Density.dat",ave_vecr,non->singles);
        printf("Completed reading pre-calculated data.\n");
        
	    // ave_vecr is average density matrix, J is full average intersegment coupling
   	    compute_all_traces_4th_order(ave_vecr, J, non); 
        
	printf("Done with computing the 4th order traces.\n");
    }

    /* Call the full 4th order routine */
    if (!strcmp(non->technique, "MCFRET-4th-full")){
        printf("Starting calculation of the full 4th order rates.\n");
        printf("Calculating 4th order from precalculated coupling and density\n");

        read_matrix_from_file("CouplingMCFRET.dat",J,non->singles);
        read_matrix_from_file("Average_Density.dat",ave_vecr,non->singles);
        printf("Completed reading pre-calculated data.\n");

        // ave_vecr is average density matrix, J is full average intersegment coupling
        full_4th_order_main(ave_vecr,J,non);

        printf("Done with computing the full 4th order rates.\n");
    }

    /* Call the MCFRET Analyse routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Analyse"))){    
        printf("Starting analysis of the MCFRET rate.\n");
	/* If analysis is done as post processing first read the rate matrix */
        if ((!strcmp(non->technique, "MCFRET-Analyse"))){
            read_matrix_from_file("RateMatrix.dat",rate_matrix,segments);
        }

	/* Define various arrays */
	re_e=(float *)calloc(segments,sizeof(float));
	im_e=(float *)calloc(segments,sizeof(float));
	vl=(float *)calloc(segments*segments,sizeof(float));
	vr=(float *)calloc(segments*segments,sizeof(float));
	energy_cor=(float *)calloc(segments,sizeof(float));
	/* Find Eigenvalues and vectors */
	mcfret_eigen(non,rate_matrix,re_e,im_e,vl,vr,segments,energy_cor);
        /* Calculate the expectation value of the segment energies */
        mcfret_energy(E,non,segments, ave_vecr,energy_cor);
        /* Analyse the rate matrix */
        mcfret_analyse(E,rate_matrix,non,segments);	
        free(re_e),free(im_e),free(vl),free(vr);	
	free(energy_cor);
    }


    free(re_Abs);
    free(im_Abs);
    free(re_Emi);
    free(im_Emi);
    free(J);
    free(E);
    free(rate_matrix);
    free(coherence_matrix);
    return;
}

// Segmented Hamiltonian propagation
/* Calculate Absorption matrix */
void mcfret_propagation_segmented(float *re_S_1,float *im_S_1,t_non *non){
    /* Define variables and arrays */
    /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int x,ti,tj,i,j;
    int t1;
    int elements;
    int cl,Ncl;
    int N_segments;

    /* Hamiltonian of the whole system - all donors and acceptors included */
    float *Hamil_i_e;
    /* Vectors representing time dependent states: real and imaginary part */
    float *vecr, *veci;
    /* Transition dipoles for coupling on the fly */
    float *mu_xyz;
    float shift1; 

    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;
    
    /* File handles */
    FILE *H_traj;
    FILE *C_traj;
    FILE *mu_traj;
    FILE *log;
    FILE *Cfile;
    FILE *absorption_matrix; 

    /* Allocate memory for all the variables */
    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));
    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Here we want to call the routine for checking the trajectory files */ 
    control(non);

    itime=0;

    /* Initialize sample numbers */
    N_samples=determine_samples(non);
    N_segments=project_dim(non);
    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);

    /* Read coupling, this is done if the coupling and transition-dipoles are */
    /* time-independent and only one snapshot is stored */
    read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

    clearvec(re_S_1,non->singles*non->singles*non->tmax1);
    
    int N_i, N, si;
    int *H_indices_si;
    N = non->singles;
    H_indices_si = (int *)calloc(N,sizeof(int));
    
    for (si=0;si<N_segments;si++){
            // find the full matrix indices for this segment
	    N_i = find_H_indices_segment(non->psites, H_indices_si, si, non);
	    printf("Segment size: %d\n",N_i); 
	    /* Allocating memory for the real and imaginary part of the wave function that we need to propagate */
        float *U_re, *U_im;
        float *U_re_snap, *U_im_snap;
	    U_re=(float *)calloc(N_i*N_i,sizeof(float));	
	    U_im=(float *)calloc(N_i*N_i,sizeof(float));
	    U_re_snap=(float *)calloc(N_i*N_i,sizeof(float));	
	    U_im_snap=(float *)calloc(N_i*N_i,sizeof(float));
        float *work_re_si, *work_im_si;
        work_re_si =(float *)calloc(N_i*N_i,sizeof(float));
        work_im_si =(float *)calloc(N_i*N_i,sizeof(float));
 
	    /* The segment si hamiltonian in upper triangle format */
            float *Hamiltonian_segment_triu;
	    Hamiltonian_segment_triu = (float *)calloc((N_i+1)*N_i/2,sizeof(float));

            /* Looping over samples: Each sample represents a different starting point on the Hamiltonian trajectory */
	    for (samples=non->begin;samples<non->end;samples++){
		ti=samples*non->sample;
		if (non->cluster!=-1){
		    if (read_cluster(non,ti,&cl,Cfile)!=1){
			    printf("Cluster trajectory file to short, could not fill buffer!!!\n");
			    printf("ITIME %d\n",ti);
			    exit(1);
		    }
		    /* Configuration belong to cluster */ 
		    if (non->cluster==cl){
			    Ncl++;
		    }
		}
		unitmat(U_re,N_i);
		clearvec(U_im,N_i*N_i);
		
		/* Loop over delay */ 
	        for (t1=0;t1<non->tmax1;t1++){
      		    tj=ti+t1;
		        /* Read Hamiltonian */
		        read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
	
	            // isolate segment Hamiltonian
		        isolate_segment_Hamiltonian_triu(Hamil_i_e, Hamiltonian_segment_triu, H_indices_si, N_i, non);
			
		        /* Update the MCFRET 'absorpion matrix' or propagator */
		        mcfret_response_function_sub_segments(re_S_1, im_S_1,t1,non,U_re,U_im,H_indices_si, N_i);

           	    // segmented 'coupling' propagation scheme
                if (non->propagation==1){
		            propagate_matrix_segments(non,Hamiltonian_segment_triu,U_re,U_im,-1,samples,tj*x, N_i);
                }
                else{
            	    // Full diagonal propagation routine
                    if (t1==0 && si==0 && samples==non->begin){
                        printf("Using full propagation scheme, with MKL & OPENBLAS.\n");
                    }
                    time_evolution_mat_non_sparse(non, Hamiltonian_segment_triu, U_re_snap, U_im_snap, N_i);
            	    propagate_snapshot(U_re_snap, U_im_snap, &U_re, &U_im, &work_re_si, &work_im_si, N_i);
                }
	       } /* We are closing the loop over time delays - t1 times */

	    /* Update NISE log file */ 
	    log=fopen("NISE.log","a");
	    fprintf(log,"Finished sample %d\n",samples);
		  
	    time_now=log_time(time_now,log);
	    fclose(log);
	    }/* Closing the loop over samples */
    	free(U_re);
	    free(U_im);
    	free(U_re_snap);
	    free(U_im_snap);
    	free(work_re_si);
	    free(work_im_si);
        free(Hamiltonian_segment_triu);
    }

    /* The calculation is finished, lets write output */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished Calculating MCFRET segmented propagators!\n");
    fprintf(log,"Writing to file!\n");  
    fclose(log);

    if (non->cluster!=-1){
        printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
        if (samples==0){ /* Avoid dividing by zero */ 
            samples=1;
        }
    }

    /* Normalize response */
    for (t1=0;t1<non->tmax1*non->singles*non->singles;t1++){
        re_S_1[t1]=re_S_1[t1]/samples;
        im_S_1[t1]=im_S_1[t1]/samples;
    }

    fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }

    /* Save time domain response */
    absorption_matrix=fopen("TD_absorption_matrix.dat","w");
    fprintf(absorption_matrix,"Samples %d\n",samples);
    fprintf(absorption_matrix,"Dimension %d\n",non->singles*non->singles*non->tmax1);
    for (t1=0;t1<non->tmax1;t1++){
        fprintf(absorption_matrix,"%f ",t1*non->deltat);
	    for (i=0;i<non->singles;i++){
	        for (j=0;j<non->singles;j++){
	            fprintf(absorption_matrix,"%e %e ",re_S_1[t1*non->singles*non->singles+i*non->singles+j],im_S_1[t1*non->singles*non->singles+i*non->singles+j]);
	        }
	    }
	    fprintf(absorption_matrix,"\n");
    }
    fclose(absorption_matrix);
    
    /*Free the memory*/
    free(Hamil_i_e);
    free(mu_xyz);
    free(H_indices_si);
}


// full Hamiltonian propagation
/* Calculate Absorption/Emission matrix (depending on emission variable 0/1) */
void mcfret_response_function(float *re_S_1,float *im_S_1,t_non *non,int emission,float *ave_vecr){
    /* Define variables and arrays */
    /* Integers */
    int nn2;
    int itime,N_samples;
    int samples;
    int x,ti,tj,i,j;
    int t1;
    int elements;
    int cl,Ncl;
    int segments;

    /* Hamiltonian of the whole system - all donors and acceptors included */
    float *Hamil_i_e;
    /* Vectors representing time dependent states: real and imaginary part */
    float *vecr, *veci;
    /* Transition dipoles for coupling on the fly */
    float *mu_xyz;
    float shift1; 

    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;

    /* File handles */
    FILE *H_traj;
    FILE *C_traj;
    FILE *mu_traj;
    FILE *log;
    FILE *Cfile;
    FILE *absorption_matrix; 

    /* Allocate memory for all the variables */
    /* Allocating memory for the real and imaginary part of the wave function that we need to propagate */
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));	
    veci=(float *)calloc(non->singles*non->singles,sizeof(float));
    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));
    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Here we want to call the routine for checking the trajectory files */ 
    control(non);

    itime=0;

    /* Initialize sample numbers */
    N_samples=determine_samples(non);
    segments=project_dim(non);
    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);

    /* Read coupling, this is done if the coupling and transition-dipoles are */
    /* time-independent and only one snapshot is stored */
    read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);

    clearvec(re_S_1,non->singles*non->singles*non->tmax1);
    /* Looping over samples: Each sample represents a different starting point on the Hamiltonian trajectory */
    for (samples=non->begin;samples<non->end;samples++){
        ti=samples*non->sample;
        if (non->cluster!=-1){
            if (read_cluster(non,ti,&cl,Cfile)!=1){
	            printf("Cluster trajectory file to short, could not fill buffer!!!\n");
	            printf("ITIME %d\n",ti);
	            exit(1);
            }
            /* Configuration belong to cluster */ 
            if (non->cluster==cl){
	            Ncl++;
            }
        }
        if (non->cluster==-1 || non->cluster==cl){
            /* Initialize time-evolution operator */
            if (emission==1){
                /* Read Hamiltonian */
                read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
	          
                /* Remove couplings between segments */
                multi_projection_Hamiltonian(Hamil_i_e,non);

                /* Use the provided density matrix as initial state */
                copyvec(ave_vecr,vecr,non->singles*non->singles);
            } else { 
                unitmat(vecr,non->singles);
                /* write_matrix_to_file("Unit.dat",vecr,non->singles); */
            }
            clearvec(veci,non->singles*non->singles);
        
            /* Loop over delay */ 
            for (t1=0;t1<non->tmax1;t1++){
	        tj=ti+t1;
	        /* Read Hamiltonian */
                read_Hamiltonian(non,Hamil_i_e,H_traj,tj);
	          
                /* Remove couplings between segments */
                multi_projection_Hamiltonian(Hamil_i_e,non);
                
                /* Update the MCFRET Response */
                mcfret_response_function_sub(re_S_1, im_S_1,t1,non,vecr,veci);        
                if (emission==0){         
                    propagate_matrix(non,Hamil_i_e,vecr,veci,-1,samples,t1*x);
                } else {
		            propagate_matrix(non,Hamil_i_e,vecr,veci,1,samples,t1*x);
		}	   
            }/* We are closing the loop over time delays - t1 times */
        } /* We are closing the cluster loop */

        /* Update NISE log file */ 
        log=fopen("NISE.log","a");
        fprintf(log,"Finished sample %d\n",samples);
          
        time_now=log_time(time_now,log);
        fclose(log);
    }/* Closing the loop over samples */
    
    /* The calculation is finished, lets write output */
    log=fopen("NISE.log","a");
    if (emission==1){
        fprintf(log,"Finished Calculating Emission Response Matrix!\n");
    } else {
	fprintf(log,"Finished Calculating Absorption Response Matrix!\n");
    }
    fprintf(log,"Writing to file!\n");  
    fclose(log);

    if (non->cluster!=-1){
        printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
        if (samples==0){ /* Avoid dividing by zero */ 
            samples=1;
        }
    }

    /* Normalize response */
    for (t1=0;t1<non->tmax1*non->singles*non->singles;t1++){
        re_S_1[t1]=re_S_1[t1]/samples;
        im_S_1[t1]=im_S_1[t1]/samples;
    }

    fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }

    /* Save time domain response */
    if (emission==0){ 
        absorption_matrix=fopen("TD_absorption_matrix.dat","w");
    } else {
        absorption_matrix=fopen("TD_emission_matrix.dat","w");
    }
    fprintf(absorption_matrix,"Samples %d\n",samples);
    fprintf(absorption_matrix,"Dimension %d\n",non->singles*non->singles*non->tmax1);
    for (t1=0;t1<non->tmax1;t1++){
        fprintf(absorption_matrix,"%f ",t1*non->deltat);
	    for (i=0;i<non->singles;i++){
	        for (j=0;j<non->singles;j++){
	            fprintf(absorption_matrix,"%e %e ",re_S_1[t1*non->singles*non->singles+i*non->singles+j],im_S_1[t1*non->singles*non->singles+i*non->singles+j]);
	        }
	    }
	    fprintf(absorption_matrix,"\n");
    }
    fclose(absorption_matrix);
    
    /*Free the memory*/
    free(vecr);	
    free(veci);  
    free(Hamil_i_e);
    free(mu_xyz);
}


/* Sub routine for adding up the calculated response in the response function */
void mcfret_response_function_sub_segments(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci, int *H_indices_si,int N_i){
    int i1,i2;
    int N,nn2;
    int Ni,tnn;
    int H_a, H_b, N_ref;
    N=non->singles;
    nn2=N*N;
    tnn=t1*nn2;
    /* Update response matrix */
    for (i1=0;i1<N_i;i1++){
        H_a = H_indices_si[i1];
	    N_ref = N * H_a;
        for (i2=0; i2<N_i; i2++){
            H_b = H_indices_si[i2];

            /* We store response function so we can do matrix multiplication */ 
	        re_S_1[tnn+(N_ref+H_b)]+=cr[i1 * N_i + i2];
            im_S_1[tnn+(N_ref+H_b)]+=ci[i1 * N_i + i2];
        }
    }
    return;
}

/* Sub routine for adding up the calculated response in the response function */
void mcfret_response_function_sub(float *re_S_1,float *im_S_1,int t1,t_non *non,float *cr,float *ci){
    int i,k;
    int N,nn2;
    int Ni,tnn;
    N=non->singles;
    nn2=N*N;
    tnn=t1*nn2;

    /* Update response matrix */
    for (i=0;i<N;i++){
        Ni=N*i;
        for (k=0; k<N; k++){
            /* We store response function so we can do matrix multiplication */
            re_S_1[tnn+(Ni+k)]+=cr[Ni+k]; 
            im_S_1[tnn+(Ni+k)]+=ci[Ni+k];
        }
    }
    return;
}

/* Find the average couplings but only between different segments */
void mcfret_coupling(float *J,t_non *non){
    /* Define variables and arrays */
    /* Integers */
    int N;
    int nn2;
    int N_samples;
    int samples;
    int ti,i,j;
    int cl,Ncl;

    /* Hamiltonian of the whole system - all donors and acceptors included */
    float *Hamil_i_e;
    float *mu_xyz;
    float shift1;
    float invsamp;

    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;

    /* File handles */
    FILE *H_traj;
    FILE *C_traj;
    FILE *mu_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;

    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));

    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Here we want to call the routine for checking the trajectory files */ 
    control(non);

    /* Initialize sample numbers */
    N_samples=determine_samples(non);
    Ncl=0;

    /* Read coupling, this is done if the coupling and transition-dipoles are */
    /* time-independent and only one snapshot is stored */
    read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);
    samples=1;


    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);
    N=non->singles;

    /* Looping over samples: Each sample represents a different starting point on the Hamiltonian trajectory */
    for (samples=non->begin;samples<non->end;samples++){ 
        ti=samples*non->sample;
        if (non->cluster!=-1){
            if (read_cluster(non,ti,&cl,Cfile)!=1){
	            printf("Cluster trajectory file to short, could not fill buffer!!!\n");
	            printf("ITIME %d\n",ti);
	            exit(1);
            }
            /* Configuration belong to cluster */ 
            if (non->cluster==cl){
	            Ncl++;
            }
        }
        if (non->cluster==-1 || non->cluster==cl){
            /* Read Hamiltonian */
            read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
          	    
            /* Extract couplings between segments */
            multi_projection_Coupling(Hamil_i_e,non);
            for (i=0;i<N;i++){
                for (j=i+1;j<N;j++){
                    // J[non->singles*i+j]+=Hamil_i_e[Sindex(i,j,non->singles)];
                    // J[non->singles*j+i]+=Hamil_i_e[Sindex(i,j,non->singles)];
                    J[N*i+j]+=Hamil_i_e[j+i*((N*2)-i-1)/2];
                    J[N*j+i]+=Hamil_i_e[j+i*((N*2)-i-1)/2];
                }
            }
        } /* We are closing the cluster loop */

        /* Update NISE log file */ 
        log=fopen("NISE.log","a");
        fprintf(log,"Finished sample %d\n",samples);
          
        time_now=log_time(time_now,log);
        fclose(log);
    }/* Closing the loop over samples */
    
    /* Divide with total number of samples */
    invsamp=1.0/samples;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            J[N*i+j]=J[N*i+j]*invsamp;
        }
    }
    write_matrix_to_file("CouplingMCFRET.dat",J,non->singles);

    /* The calculation is finished, lets write output */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished Averaging Intersegment Coupling Response!\n");
    fprintf(log,"Writing to file!\n");  
    fclose(log);

    if (non->cluster!=-1){
        printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
        if (samples==0){ /* Avoid dividing by zero */ 
            samples=1;
        }
    }

    fclose(mu_traj),fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }
    free(Hamil_i_e);
    free(mu_xyz);  
    return;

}

/* Find MCFRET segments using an automatic scheme */
void mcfret_autodetect(t_non *non, float treshold){
    printf("Use the analyse technique for auto detection.\n");
    return;
}

/* Calculate actual rate matrix */
void mcfret_rate_from_abs(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs, float *rho_0,float *J,t_non *non){
    
    int nn2,N;
    int si,sj;
    int i,j,k;
    int *ns; /* Segment dimensions */
    int t1;
    float *rate_response, *rate_response_imag, *abs_rate_response;
    float rate;
    float isimple,is13; /* Variables for integrals */
    float *re_Abs_hermi,*im_Abs_hermi;
    float *re_aux_mat,*im_aux_mat;
    float *re_aux_mat2,*im_aux_mat2;
    float *Zeros;
    float twoPi2;
    float trace_reaux, trace_imaux;
    FILE *ratefile, *ratefile_imag;
    FILE *log;

    N=non->singles;
    nn2=non->singles*non->singles;
    twoPi2=twoPi*twoPi;

    rate_response=(float *)calloc(non->tmax,sizeof(float));
    rate_response_imag=(float *)calloc(non->tmax,sizeof(float));
    abs_rate_response=(float *)calloc(non->tmax,sizeof(float));
    re_Abs_hermi=(float *)calloc(nn2,sizeof(float));
    im_Abs_hermi=(float *)calloc(nn2,sizeof(float)); 
    re_aux_mat=(float *)calloc(nn2,sizeof(float));
    im_aux_mat=(float *)calloc(nn2,sizeof(float));
    re_aux_mat2=(float *)calloc(nn2,sizeof(float));
    im_aux_mat2=(float *)calloc(nn2,sizeof(float));
    Zeros=(float *)calloc(nn2,sizeof(float));
  
    ratefile=fopen("RateFile.dat","w");
    ratefile_imag=fopen("RateFile_imag.dat","w");
    /* Do one rate at a time - so first we loop over segments */
    /* Tr [ J * Abs(t1) * J * Emi(t1) ] */
    for (si=0;si<segments;si++){
        for (sj=0;sj<segments;sj++){
            /* Exclude rate between same segments */
            if (sj!=si){
                /* Loop over time delay */
                for (t1=0;t1<non->tmax;t1++){   
                    log=fopen("NISE.log","a");
                    fprintf(log,"Between segments %d and %d, computing response at t = %d \n",si, sj, t1);
                    fclose(log);

                    /* compute the hermitian conjugate of the absorption matrix */
                    hermitian_conjugate(re_Abs+nn2*t1,im_Abs+nn2*t1,re_Abs_hermi,im_Abs_hermi,N,N);

                    /* Matrix multiplication - J Abs_hermi */
                    segment_matrix_mul(J,Zeros,re_Abs_hermi,im_Abs_hermi,
                    re_aux_mat,im_aux_mat,non->psites,segments,si,sj,sj,N);

                    /* Matrix multiplication - (J Abs_hermi) rho_0 */
                    segment_matrix_mul(re_aux_mat,im_aux_mat,rho_0,Zeros,
                    re_aux_mat2,im_aux_mat2,non->psites,segments,si,sj,sj,N);


                    // printf("matrix check  %f\n",matrix_sum(J,N));
                    // /* Matrix multiplication - J Emi */
                    // segment_matrix_mul(J,Zeros,re_Emi+nn2*t1,im_Emi+nn2*t1,
                    // re_aux_mat,im_aux_mat,non->psites,segments,si,sj,sj,N);

                    /* Matrix multiplication - Abs (J Abs_hermi rho_0) */
                    segment_matrix_mul(re_Abs+nn2*t1,im_Abs+nn2*t1,re_aux_mat2,im_aux_mat2,
                    re_aux_mat,im_aux_mat,non->psites,segments,si,si,sj,N);
                    /* Matrix multiplication - J  (Abs J Abs_hermi rho_0) */
                    segment_matrix_mul(J,Zeros,re_aux_mat,im_aux_mat,
                    re_aux_mat2,im_aux_mat2,non->psites,segments,sj,si,sj,N);
                    /* Take the trace */
                    trace_reaux=trace_rate(re_aux_mat2,N);
                    trace_imaux=trace_rate(im_aux_mat2,N);


	                rate_response[t1]=trace_reaux*twoPi2;
	                rate_response_imag[t1] = trace_imaux * twoPi2;
		            abs_rate_response[t1]=sqrt(trace_reaux*trace_reaux
                                    +trace_imaux*trace_imaux)*twoPi2;
                    fprintf(ratefile,"%f %f\n",t1*non->deltat,rate_response[t1]);
                    fprintf(ratefile_imag,"%f %f\n",t1*non->deltat,rate_response_imag[t1]);
                }
                /* Update rate matrix */
	            integrate_rate_response(rate_response,non->tmax,&is13,&isimple);
		        /* We use the Trapezium, which is most accurate in most cases */
                rate=2*isimple*non->deltat*icm2ifs*icm2ifs*1000;
                rate_matrix[si*segments+sj]=rate;
                rate_matrix[sj*segments+sj]-=rate;
	            /* Calculate the rate of coherence decay in ps-1 */
	            integrate_rate_response(abs_rate_response,non->tmax,&is13,&isimple);
	            coherence_matrix[si*segments+sj]=1000*abs_rate_response[0]/isimple/non->deltat;
            }
        }
    }
    fclose(ratefile);
    fclose(ratefile_imag);

    free(rate_response);
    free(rate_response_imag);
    free(abs_rate_response);
    free(re_Abs_hermi);
    free(im_Abs_hermi);
    free(re_aux_mat);
    free(im_aux_mat);
    free(re_aux_mat2);
    free(im_aux_mat2);
    free(Zeros);
    free(ns);
    return;
}

/* Find Hermitian conjugate of square NxN matrix */
void hermitian_conjugate(float *A_re, float *A_im, float *hermi_re, float *hermi_im, int N1, int N2){
    int a,b;
    clearvec(hermi_re,N1*N2);
    clearvec(hermi_im,N1*N2);

    for (a=0;a<N1;a++){
        for (b=0;b<N2;b++){
            hermi_re[a+b*N1] = A_re[b+a*N2];
            hermi_im[a+b*N1] = -A_im[b+a*N2];
        }
    }
    return;
}

/* Calculate actual rate matrix */
void mcfret_rate(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs,
    float *re_Emi,float *im_Emi,float *J,t_non *non){
    
    int nn2,N;
    int si,sj;
    int i,j,k;
    int *ns; /* Segment dimensions */
    int t1;
    float *rate_response, *rate_response_imag, *abs_rate_response;
    float rate;
    float isimple,is13; /* Variables for integrals */
    float *re_aux_mat,*im_aux_mat;
    float *re_aux_mat2,*im_aux_mat2;
    float *Zeros;
    float twoPi2;
    float trace_reaux, trace_imaux;
    FILE *ratefile, *ratefile_imag;
    N=non->singles;
    nn2=non->singles*non->singles;
    twoPi2=twoPi*twoPi;

    rate_response=(float *)calloc(non->tmax,sizeof(float));
    rate_response_imag=(float *)calloc(non->tmax,sizeof(float));
    abs_rate_response=(float *)calloc(non->tmax,sizeof(float));
    re_aux_mat=(float *)calloc(nn2,sizeof(float));
    im_aux_mat=(float *)calloc(nn2,sizeof(float));
    re_aux_mat2=(float *)calloc(nn2,sizeof(float));
    im_aux_mat2=(float *)calloc(nn2,sizeof(float));
    Zeros=(float *)calloc(nn2,sizeof(float));
  
    ratefile=fopen("RateFile.dat","w");
    ratefile_imag=fopen("RateFile_imag.dat","w");
    /* Do one rate at a time - so first we loop over segments */
    /* Tr [ J * Abs(t1) * J * Emi(t1) ] */
    for (si=0;si<segments;si++){
        for (sj=0;sj<segments;sj++){
            /* Exclude rate between same segments */
            if (sj!=si){
                /* Loop over time delay */
                for (t1=0;t1<non->tmax;t1++){
                    /* Matrix multiplication - J Emi */
                    segment_matrix_mul(J,Zeros,re_Emi+nn2*t1,im_Emi+nn2*t1,
                    re_aux_mat,im_aux_mat,non->psites,segments,si,sj,sj,N);
                    /* Matrix multiplication - Abs (J Emi) */
                    segment_matrix_mul(re_Abs+nn2*t1,im_Abs+nn2*t1,re_aux_mat,im_aux_mat,
                    re_aux_mat2,im_aux_mat2,non->psites,segments,si,si,sj,N);
                    /* Matrix multiplication - J (Abs J Emi) */
                    segment_matrix_mul(J,Zeros,re_aux_mat2,im_aux_mat2,
                    re_aux_mat,im_aux_mat,non->psites,segments,sj,si,sj,N);
                    /* Take the trace */
                    trace_reaux=trace_rate(re_aux_mat,N);
                    trace_imaux=trace_rate(im_aux_mat,N);
                    rate_response[t1]=trace_reaux*twoPi2;
	            rate_response_imag[t1] = trace_imaux * twoPi2;
		    abs_rate_response[t1]=sqrt(trace_reaux*trace_reaux
                                    +trace_imaux*trace_imaux)*twoPi2;
                    fprintf(ratefile,"%f %f\n",t1*non->deltat,rate_response[t1]);
                    fprintf(ratefile_imag,"%f %f\n",t1*non->deltat,rate_response_imag[t1]);
                }
                /* Update rate matrix */
	        integrate_rate_response(rate_response,non->tmax,&is13,&isimple);
		/* We use the Trapezium, which is most accurate in most cases */
                rate=2*isimple*non->deltat*icm2ifs*icm2ifs*1000;
                rate_matrix[si*segments+sj]=rate;
                rate_matrix[sj*segments+sj]-=rate;
	        /* Calculate the rate of coherence decay in ps-1 */
	        integrate_rate_response(abs_rate_response,non->tmax,&is13,&isimple);
	        coherence_matrix[si*segments+sj]=1000*abs_rate_response[0]/isimple/non->deltat;
            }
        }
    }
    fclose(ratefile);
    fclose(ratefile_imag);

    free(rate_response);
    free(rate_response_imag);
    free(abs_rate_response);
    free(re_aux_mat);
    free(im_aux_mat);
    free(re_aux_mat2);
    free(im_aux_mat2);
    free(Zeros);
    free(ns);
    return;
}

/* Check if mcfret rates are in the incoherent limit */
void mcfret_validate(t_non *non);

/* Find Eigenvalues and eigenvectors of rate matrix */
void mcfret_eigen(t_non *non,float *rate_matrix,float *re_e,float *im_e,float *vl,float *vr,int segments,float *energy_cor){
    //char jobvl = 'V';  // Compute left eigenvectors
    //char jobvr = 'V';  // Compute right eigenvectors
    //int lwork = segments * segments;  // Work array size
    //float work[lwork];
    float *rate; /* Rate Matrix to be destroyed */
    int *degen; /* Degeneracies of segments */
    int i;
    int imax;
    float fmax;
    float popnorm;
    FILE *Efile;
    float *ivr,*ivl;

    rate=(float *)calloc(segments*segments,sizeof(float));
    degen=(int *)calloc(segments,sizeof(int));
    ivr=(float *)calloc(segments*segments,sizeof(float));
    ivl=(float *)calloc(segments*segments,sizeof(float));

    copyvec(rate_matrix,rate,segments*segments);

    diagonalize_real_nonsym(rate_matrix,re_e,im_e,vl,ivl,vr,ivr,segments);

    /* Call LAPACK function sgeev to compute eigenvalues and eigenvectors */
    //sgeev_(&jobvl, &jobvr, &segments, rate, &segments, re_e, im_e, vl, &segments, vr, &segments, work, &lwork, &info);
    //free(rate);

    /* Check for imaginary eigenvalues and find the one closest to zero */
    imax=0;
    fmax=re_e[0];
    
    for (i=0;i<segments;i++){
	    /* Weak check */
	    if (fabs(im_e[i])>0.1*fabs(re_e[i])){
            printf(RED "An imaginary rate matrix eigenvalue larger than 10%%\n");
	        printf("of the real value found! Averaging over more relaizations\n");
	        printf("is adviseable. Use rate matrix with caution!\n" RESET);
	        exit(0);

	    /* Hard Check */
	    } else if (fabs(im_e[i])>1e-8) {
            printf(YELLOW "Warning! An imaginary eigenvalue of the rate matrix was found.");
            printf("Check validity. Averaging over more disorder realizations\n");
	        printf("may remove imaginary eigenvalues." RESET);
	    }

	    /* Check if it is larger than the previous ones */
	    if (re_e[i]>fmax){
            fmax=re_e[i];
	        imax=i;
	    }
    }
    
    /* Write eigenvalues to file and find normalization for the */
    /* equilibrium population. */
    popnorm=0;
    Efile=fopen("RateMatrixEigenvalues.dat","w");
    fprintf(Efile,"# - Eigenvalue (real and imaginary part in ps-1) \n");
    for (i=0;i<segments;i++){
        fprintf(Efile,"%d %f %f\n",i,re_e[i],im_e[i]);
	popnorm+=vl[i+segments*imax];
    }
    fclose(Efile);

    /* Find number of degeneracies */
    project_degeneracies(non,degen,segments);

    /* Write Segment Equilibrium Populations to file */
    /* and find effective energy correction to give equal populations */
    Efile=fopen("SegmentPopulation.dat","w");
    fprintf(Efile,"# - Equilibrium Population\n");
    for (i=0;i<segments;i++){
        fprintf(Efile,"%d %f\n",i,vl[i+segments*imax]/popnorm);
	/* Skip adjusting quantum correction in the high-temperature limit */
	if (non->temperature<100000){
	    energy_cor[i]=-non->temperature*k_B*logf(vl[i+segments*imax]/popnorm/degen[i]);
    }
    }
    fclose(Efile);

    write_matrix_to_file("LeftVectorRateMatrix.dat",vl,segments);
    write_matrix_to_file("RightVectorRateMatrix.dat",vr,segments);
    free(degen);
    free(ivr);
    free(ivl);
    return;
}

/* Analyse rate matrix and find thermal correction */
void mcfret_analyse(float *E,float *rate_matrix,t_non *non,int segments){
      float *qc_rate_matrix,*qc;
      /* Thermal correction */
      float *tc_rate_matrix;
      float *partition_functions;
      float partition_function_i, partition_function_j, equilibration_rate;
      float C;
      float column;
      int i,j;
      float kBT=non->temperature*k_B; /* Kelvin to cm-1 */                     
  
      /* Allocate memory for the partition functions and rate matrices */
      qc_rate_matrix=(float *)calloc(segments*segments,sizeof(float));
      qc=(float *)calloc(segments*segments,sizeof(float)); 
      tc_rate_matrix=(float *)calloc(segments*segments,sizeof(float));
      partition_functions = (float *)calloc(segments,sizeof(float));        
      
      /* load the segment ensemble avg partition functions */
      read_vector_from_file("Segment_Partition_Functions.dat",partition_functions,segments); 

      /* Find quantum correction factors */                                    
      for (i=0;i<segments;i++){
          partition_function_i = partition_functions[i];
	  column = 0;
          for (j=0;j<segments;j++){                                            
              if (i!=j){
                  /* Quantum correction factor from D.W. Oxtoby. */
                  /* Annu. Rev. Phys. Chem., 32(1):77–101, (1981).*/       
                  C=2/(1+exp((E[i]-E[j])/kBT));                            
                  qc_rate_matrix[i*segments+j]=rate_matrix[i*segments+j]*C;
                  qc_rate_matrix[j*segments+j]-=rate_matrix[i*segments+j]*C;
                  qc[i*segments+j]=C;
		  /* Partition Function Based Thermal Correction */
		  equilibration_rate = rate_matrix[i*segments+j] + rate_matrix[j*segments+i];
		  partition_function_j = partition_functions[j];
		  /* rate from segment i to segment j */
		  tc_rate_matrix[i+segments*j] = equilibration_rate * partition_function_j / (partition_function_i + partition_function_j);
                  column += tc_rate_matrix[i+segments*j];
	      }       
              else{ 
                  qc[i*segments+j]=0;
              }       
          }
	  /* Ensure population conservation in the TC_rate_matrix */
	  tc_rate_matrix[i+segments*i] = -column;
      }
  
      /* Write the quantum corrected rate matrix. */
      write_matrix_to_file("QC_RateMatrix.dat",qc_rate_matrix,segments);
      /* Write the thermal corrected rate matrix. */
      write_matrix_to_file("TC_RateMatrix.dat",tc_rate_matrix,segments);
      /* Write the applied quantum correction factors. */
      write_matrix_to_file("QC.dat",qc,segments);                              
      free(qc_rate_matrix);
      free(tc_rate_matrix);
      free(qc);
      return;                                                                  
  }

/* Find the energy of each segment */
void mcfret_energy(float *E,t_non *non,int segments, float *ave_vecr,float *energy_cor){
    /* Define variables and arrays */
    /* Integers */
    int nn2;
    int N_samples;
    int samples;
    int ti,i,j;
    int cl,Ncl;

    /* Hamiltonian of the whole system - all donors and acceptors included */
    float *Hamil_i_e;
    float *vecr;
    float *mu_xyz;
    float shift1;

    /* Time parameters */
    time_t time_now,time_old,time_0;
    /* Initialize time */
    time(&time_now);
    time(&time_0);
    shift1=(non->max1+non->min1)/2;
    printf("Frequency shift %f.\n",shift1);
    non->shifte=shift1;
    /* File handles */
    FILE *H_traj;
    FILE *C_traj;
    FILE *mu_traj;
    /* FILE *absorption_matrix, *emission_matrix,*/
    FILE *log;
    FILE *Cfile;
    FILE *Efile;

    mu_xyz=(float *)calloc(non->singles*3,sizeof(float));
    Hamil_i_e=(float *)calloc((non->singles+1)*non->singles/2,sizeof(float));
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));

    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Here we want to call the routine for checking the trajectory files */
    control(non);

    /* Initialize sample numbers */
    N_samples=determine_samples(non);
    Ncl=0;

    /* Read coupling, this is done if the coupling and transition-dipoles are */
    /* time-independent and only one snapshot is stored */
    read_coupling(non,C_traj,mu_traj,Hamil_i_e,mu_xyz);
    samples=1;

    log=fopen("NISE.log","a");
    fprintf(log,"Begin sample: %d, End sample: %d.\n",non->begin,non->end);
    fclose(log);


    /* Looping over samples: Each sample represents a different starting point on the Hamiltonian trajectory */
    for (samples=non->begin;samples<non->end;samples++){
	ti=samples*non->sample;
        if (non->cluster!=-1){
            if (read_cluster(non,ti,&cl,Cfile)!=1){
                printf("Cluster trajectory file to short, could not fill buffer!!!\n");
                printf("ITIME %d\n",ti);
                exit(1);
            }
            /* Configuration belong to cluster */ 
            if (non->cluster==cl){
                Ncl++;
            }
        }
        if (non->cluster==-1 || non->cluster==cl){
            /* Read Hamiltonian */
            read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
    
            /* Remove couplings between segments */
            multi_projection_Hamiltonian(Hamil_i_e,non);	    
            /* Find density matrix */
            copyvec(ave_vecr,vecr,non->singles*non->singles);
	    /* H * rho */
            triangular_on_square(Hamil_i_e,vecr,non->singles); 
	    /* Add energy contribution for each segment */
	    /* that is take the trace for each segment */
	    for (i=0;i<non->singles;i++){
	        E[non->psites[i]]+=vecr[i*non->singles+i];    
            }
	    
	    clearvec(vecr,non->singles*non->singles);
        } /* We are closing the cluster loop */

        /* Update NISE log file */
        log=fopen("NISE.log","a");
        fprintf(log,"Segement Energy Finished sample %d\n",samples);

        time_now=log_time(time_now,log);
        fclose(log);
    }/* Closing the loop over samples */

    /* Divide with total number of samples */
    for (i=0;i<segments;i++){
        E[i]=E[i]/N_samples;
    }
    Efile=fopen("SegmentEnergies.dat","w");
    fprintf(Efile,"# Segment number - Average segment energy - Energy correction %d\n",N_samples);
    for (i=0;i<segments;i++){
        fprintf(Efile,"%d %f %f\n",i,E[i]+shift1,energy_cor[i]);
	E[i]=E[i]-energy_cor[i];
    }
    fclose(Efile);

    /* The calculation is finished, lets write output */
    log=fopen("NISE.log","a");
    fprintf(log,"Finished Calculating Segment Energies!\n");
    fprintf(log,"Writing to file!\n");
    fclose(log);

    if (non->cluster!=-1){
        printf("Of %d samples %d belonged to cluster %d.\n",samples,Ncl,non->cluster);
        if (samples==0){ /* Avoid dividing by zero */
            samples=1;
        }
    }
    fclose(mu_traj),fclose(H_traj);
    if (non->cluster!=-1){
        fclose(Cfile);
    }
    free(Hamil_i_e);
    free(mu_xyz);
    free(vecr);
    return;
}

void transform_back_to_site(int N,
                           const float *A_float,
                           const double *c2,
                           double *matrix)
{
    int i, j;
    // Perform: matrix_site = H * matrix_eigen * H.T
    // Convert A -> double (BLAS uses double here)
    double *A = (double*) malloc(N*N*sizeof(double));
    double *temp = (double*) calloc(N*N,sizeof(double));

    for (i = 0; i < N*N; i++)
        A[i] = (double) A_float[i];

    // temp = A * diag(c2)
    // scale each column j by c2[i]
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            temp[i * N + j] = A[i * N + j] * c2[j];
        }
    }

    // matrix =  temp * A^T =  (A * D) * A^T 
    cblas_dgemm(
        CblasRowMajor,
        CblasNoTrans,    
        CblasTrans,
        N, N, N,
        1.0,
        temp, N,
        A, N,
        0.0,
        matrix, N);

    free(A);
    free(temp);
}

/* This function will create a density matrix where every term is weighted with a Boltzmann weight */
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int N_segments, float *partition_functions){
    int index,N;
    float *H,*e;
    double *c2;
    double *matrix;

    FILE *log;

    N=non->singles;

    int a,b,c;
    int si;
    int N_site_si;
    int *H_indices_si;
    H_indices_si = (int *)calloc(N,sizeof(int));
    int seg_ref_idx = 0;

    double kBT=(double) non->temperature*k_B; /* Kelvin to cm-1 */
 
    clearvec(density_matrix,N*N);

    double Boltzmann;
    double segment_partition_func = 0;

    double reference_energy;
    reference_energy = non->min1; // Simply use the first value of the minimum indicated frequencies as reference energy

    for (si=0;si<N_segments;si++){
        float *Hamiltonian_segment_triu;
        N_site_si = find_H_indices_segment(non->psites, H_indices_si, si, non);
        Hamiltonian_segment_triu = (float *)calloc((N_site_si+1)*N_site_si/2,sizeof(float));
        H=(float *)calloc(N_site_si*N_site_si,sizeof(float));
        e=(float *)calloc(N_site_si,sizeof(float));
        c2=(double *)calloc(N_site_si,sizeof(double));
        matrix=(double *)calloc(N_site_si*N_site_si,sizeof(double));

        // clear the partition function for each segment
        segment_partition_func = 0;


        // isolate the segment i with projection routine to obtain smaller matrix
        isolate_segment_Hamiltonian_triu(Hamiltonian_i, Hamiltonian_segment_triu, H_indices_si, N_site_si, non);

        /* Build segment Hamiltonian */
        for (a=0;a<N_site_si;a++){
            H[a+N_site_si*a]=Hamiltonian_segment_triu[a+N_site_si*a-(a*(a+1))/2]; /* Diagonal */
            for (b=a+1;b<N_site_si;b++){
                H[a+N_site_si*b]=Hamiltonian_segment_triu[b+N_site_si*a-(a*(a+1))/2];
                H[b+N_site_si*a]=Hamiltonian_segment_triu[b+N_site_si*a-(a*(a+1))/2];

            }
        }
        /* Find eigenvalues and eigenvectors */
        diagonalizeLPD(H,e,N_site_si);
        if (non->printLevel>1){
            log=fopen("NISE.log","a");
            fprintf(log,"Boltzmann factors for segment %d: ", si);
        }
        /* Exponentiate [U=exp(-H/kBT)] */
        for (a=0;a<N_site_si;a++){
            if (non->temperature==0){
                printf("Temperature is 0, the equilibrium density matrix will be nan, we suggest using a low non-zero temperature instead");
                exit(0);
            }

            Boltzmann = exp(-((double) (e[a]-reference_energy))/kBT);
            /* Apply strict high temperature limit when T>100000 */
            if (non->temperature>100000){
                Boltzmann=1.0;
            }
            c2[a] = Boltzmann;
            segment_partition_func += Boltzmann;
            if (non->printLevel>1){
                fprintf(log," %e ", Boltzmann);
            }
        }
        if (non->printLevel>0){
            fprintf(log,"\n" );
            fclose(log);
        }

        /* Transform back to site basis */ 
        transform_back_to_site(N_site_si, H, c2, matrix);

        /* Re-normalize and write to overarching density matrix */
        
        for (a=0;a<N_site_si;a++){
            for (b=0;b<N_site_si;b++){
                density_matrix[H_indices_si[a]+H_indices_si[b]*N]=(float) (matrix[a+b*N_site_si]/segment_partition_func);
            }
        }
    
        /* Update the ensemble average partition function for each segment*/
        partition_functions[si] += (float) segment_partition_func;
        if (non->printLevel>0){
            printf("Partition function for segment %d =  %e \n",si,segment_partition_func);
            // log=fopen("NISE.log","a");
            // fprintf(log,"Partition function for segment %d =  %e \n",si,segment_partition_func);
            // fclose(log);
        }
    
    
        free(Hamiltonian_segment_triu);
        free(H);
        free(c2);
        free(e);
        free(matrix);
    } //close the segment loop
    free(H_indices_si);


    return;
}

void average_density_matrix(float *ave_den_mat,t_non *non){
/* Define variables and arrays */
   /* Integers */
    int ti;
    int segments,s;
    int samples;
    int ele;
    int my_samples;
    int N,a,b;
    float i_samples;
    /* Vectors representing time dependent states: real and imaginary part */
    float *vecr;
    float *Hamiltonian_i;
    /* Vector representing the ensemble average partition function for each segment */
    float *partition_functions;
    /* File handles */
    FILE *H_traj;
    FILE *mu_traj;
    FILE *Cfile;
    FILE *avg_partition_functions;
    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Allocating memory for the real and imaginary part of the wave function that we need to propagate */
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
    Hamiltonian_i=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    //ave_den_mat=(float *)calloc(non->singles*non->singles,sizeof(float));
    /* Initialize sample numbers */
    segments=project_dim(non);
    N=non->singles;
  
    /* Allocate memory for the partition function vector */
    partition_functions = (float *)calloc(segments,sizeof(float));
    
    clearvec(ave_den_mat,N*N);
    /* Initialize sample numbers */
    my_samples=determine_samples(non);

    if (non->end-non->begin<my_samples){
      my_samples=non->end-non->begin;
    }
// #pragma omp parallel for
    for (samples=non->begin;samples<non->end;samples++){
      ti=samples*non->sample; 
      read_Hamiltonian(non,Hamiltonian_i,H_traj,ti);
      /* Use the thermal equilibrium as initial state */
      density_matrix(vecr,Hamiltonian_i,non,segments,partition_functions);

      for (ele=0; ele<non->singles*non->singles; ele++){
          ave_den_mat[ele] +=vecr[ele]; 
      }
    }
    /* Normalise the average density matrix */
    i_samples=1.0/my_samples;
    for (a=0;a<non->singles;a++){
        for (b=0;b<non->singles;b++){
	    ave_den_mat[non->singles*b+a]*=i_samples;
        }
    }

    /* Normalise the segments' partition funcions */
    for (s=0;s<segments;s++){
        partition_functions[s] *= i_samples;
    }

    /* Write partition function vector to a file */
    avg_partition_functions = fopen("Segment_Partition_Functions.dat","w");
    for (s=0;s<segments;s++){
        fprintf(avg_partition_functions,"%f\n",partition_functions[s]);
    }
    fclose(avg_partition_functions);

    free(partition_functions);
    free(vecr); 
    free(Hamiltonian_i); 
    return;
}

// /* Matrix multiplication for different segments */
// void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,
//     float *rC,float *iC,int *psites,int segments,int si,int sk,int sj,int N){
//     int i,j,k;
//     /* Set initial values of results matrix to zero to be sure */
//     clearvec(rC,N*N);
//     clearvec(iC,N*N);

//     int Npsj, Npsk;
//     int cj, ck;
//     Npsj=0;
//     Npsk=0;

//     float *psj, *psk;
//     psj=(float *)calloc(N,sizeof(float));
//     psk=(float *)calloc(N,sizeof(float));

//     for (i=0;i<N;i++){
//         if (psites[i]==sj){
//             psj[Npsj]=i;
//             Npsj++;
//         }
//         if (psites[i]==sk){
//             psk[Npsk]=i;
//             Npsk++;
//         }
//     }

// #pragma omp parallel for
//     for (i=0;i<N;i++){
//         if (psites[i]==si){
//             for (cj=0;cj<Npsj;cj++){
//                 j=psj[cj];
//                 for (ck=0;ck<Npsk;ck++){
//                     k=psk[ck];
//                     rC[i*N+j]+=rA[i*N+k]*rB[k*N+j]-iA[i*N+k]*iB[k*N+j];
//                     iC[i*N+j]+=rA[i*N+k]*iB[k*N+j]+iA[i*N+k]*rB[k*N+j];
//                 }
//             }
//         }
//     }

//     /*
//     for (i=0;i<N;i++){
//         if (psites[i]==si){
//             for (j=0;j<N;j++){
//                 if (psites[j]==sj){
//                     for (k=0;k<N;k++){
//                         if (psites[k]==sk){
//                             rC[i*N+j]+=rA[i*N+k]*rB[k*N+j]-iA[i*N+k]*iB[k*N+j];
//                             iC[i*N+j]+=rA[i*N+k]*iB[k*N+j]+iA[i*N+k]*rB[k*N+j];
//                         } 
//                     } 
//                 }
//             }
//         }
//     }
//     */
//     free(psj);
//     free(psk);
//     return;
// } 

/*  Faster BLAS based matrix multiplication */
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,
    float *rC,float *iC,int *psites,int segments,int si,int sk,int sj,int N){

    int i,j,k;

    clearvec(rC,N*N);
    clearvec(iC,N*N);

    int Npsi=0, Npsj=0, Npsk=0;

    int *psi = (int*)malloc(N*sizeof(int));
    int *psj = (int*)malloc(N*sizeof(int));
    int *psk = (int*)malloc(N*sizeof(int));

    for (i=0;i<N;i++){
        if (psites[i]==si) psi[Npsi++] = i;
        if (psites[i]==sj) psj[Npsj++] = i;
        if (psites[i]==sk) psk[Npsk++] = i;
    }

    /* Create submatrices based in the segment input */
    float *Ar = (float*)malloc(Npsi*Npsk*sizeof(float));
    float *Ai = (float*)malloc(Npsi*Npsk*sizeof(float));
    float *Br = (float*)malloc(Npsk*Npsj*sizeof(float));
    float *Bi = (float*)malloc(Npsk*Npsj*sizeof(float));

    /* Do the segment projection */
    for (i=0;i<Npsi;i++)
        for (k=0;k<Npsk;k++){
            Ar[i*Npsk+k] = rA[psi[i]*N + psk[k]];
            Ai[i*Npsk+k] = iA[psi[i]*N + psk[k]];
        }

    for (k=0;k<Npsk;k++)
        for (j=0;j<Npsj;j++){
            Br[k*Npsj+j] = rB[psk[k]*N + psj[j]];
            Bi[k*Npsj+j] = iB[psk[k]*N + psj[j]];
        }

    /* Allocate result blocks */
    float *Cr = (float*)calloc(Npsi*Npsj,sizeof(float));
    float *Ci = (float*)calloc(Npsi*Npsj,sizeof(float));

    /* Complex GEMM via 4 real GEMMs */

    // Cr += Ar*Br (with Cr presumed initialised at zero)
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        Npsi,Npsj,Npsk,1.0f,Ar,Npsk,Br,Npsj,1.0f,Cr,Npsj);

    // Cr -= Ai*Bi
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        Npsi,Npsj,Npsk,-1.0f,Ai,Npsk,Bi,Npsj,1.0f,Cr,Npsj);

    // Ci += Ar*Bi
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        Npsi,Npsj,Npsk,1.0f,Ar,Npsk,Bi,Npsj,1.0f,Ci,Npsj);

    // Ci += Ai*Br
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
        Npsi,Npsj,Npsk,1.0f,Ai,Npsk,Br,Npsj,1.0f,Ci,Npsj);

    /* place submatrix back in original N*N matrix */
    for (i=0;i<Npsi;i++)
        for (j=0;j<Npsj;j++){
            rC[psi[i]*N + psj[j]] = Cr[i*Npsj+j];
            iC[psi[i]*N + psj[j]] = Ci[i*Npsj+j];
        }

    free(psi); free(psj); free(psk);
    free(Ar); free(Ai); free(Br); free(Bi);
    free(Cr); free(Ci);
}

/* Find the trace of the matrix */
float trace_rate(float *matrix,int N){
    int i;
    float trace;
    trace=0;
    for (i=0;i<N;i++){
        trace=trace+matrix[N*i+i];
    }
    return trace;
}

/* Read the absorption/emission function from file */
void read_response_from_file(char fname[],float *re_R,float *im_R,int N,int tmax){
    FILE *file_handle;
    int i,j;
    int t1;
    int dummy;
    float dummyf;
    file_handle=fopen(fname,"r");
    if (file_handle == NULL) {
        printf("Error opening the file %s.\n",fname);
        exit(0);
    }
  
    /* Read initial info */
    fscanf(file_handle, "Samples %d\n", &dummy);
    fscanf(file_handle, "Dimension %d\n", &dummy);
  
    for (t1=0;t1<tmax;t1++){
        /* Read time */
        fscanf(file_handle,"%f",&dummyf);
        for (i=0;i<N;i++){
            for (j=0;j<N;j++){
	            fscanf(file_handle,"%f %f",&re_R[t1*N*N+i*N+j],&im_R[t1*N*N+i*N+j]);
            }
        }
    }
    fclose(file_handle);
}

/* Multiply a triangular matrix on a square one and return */
/* The result in the square matrix */
/* "S=T*S" */
void triangular_on_square(float *T,float *S,int N){
    float *inter;
    int a,c,b;
    int index;
    inter=(float *)calloc(N*N,sizeof(float));
    /* Do matrix multiplication */
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            for (c=0;c<N;c++){
                index=Sindex(a,b,N);
                inter[a+c*N]+=T[index]*S[b+c*N]; // TLC b -> c
            }       
        }
    }
    /* Copy result back */
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            S[a+b*N]=inter[a+b*N];
        }
    }
    free(inter);
    return;
}