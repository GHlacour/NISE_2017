#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "mcfret.h"
#include "project.h"
#include "propagate.h"
#include "read_trajectory.h"

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
        "MCFRET-Density"},8)){
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
        mcfret_response_function(re_Abs,im_Abs,non,0,ave_vecr);
    }
   
/* Call the emission routine */
    if (!strcmp(non->technique, "MCFRET") || (!strcmp(non->technique, "MCFRET-Emission"))){
        printf("Starting calculation of the MCFRET emission matrix.\n");
	/* Read precalculated average density matrix */ 
	if (!strcmp(non->technique, "MCFRET-Emission")){
            printf("Using precalculated average density matrix from file Average_Density.dat.\n");
	    read_matrix_from_file("Average_Density.dat",ave_vecr,non->singles);
	}
        mcfret_response_function(re_Emi,im_Emi,non,1,ave_vecr);
    }
    
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
	    printf("Calculating rate from precalculated absorption, emission\n");
	    printf("and coupling!\n");
	    read_matrix_from_file("CouplingMCFRET.dat",J,non->singles);
	    read_response_from_file("TD_absorption_matrix.dat",re_Abs,im_Abs,non->singles,non->tmax1);
	    read_response_from_file("TD_emission_matrix.dat",re_Emi,im_Emi,non->singles,non->tmax1);
            printf("Completed reading pre-calculated data.\n");
        }
        mcfret_rate(rate_matrix,coherence_matrix,segments,re_Abs,im_Abs,re_Emi,im_Emi,J,non);

        /* Write the calculated ratematrix to file */
        write_matrix_to_file("RateMatrix.dat",rate_matrix,segments);
        /* Write the calculated coherence matrix to file */
        write_matrix_to_file("CoherenceMatrix.dat",coherence_matrix,segments);
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
void mcfret_rate(float *rate_matrix,float *coherence_matrix,int segments,float *re_Abs,float *im_Abs,
    float *re_Emi,float *im_Emi,float *J,t_non *non){
    
    int nn2,N;
    int si,sj;
    int i,j,k;
    int *ns; /* Segment dimensions */
    int t1;
    float *rate_response,*abs_rate_response;
    float rate;
    float isimple,is13; /* Variables for integrals */
    float *re_aux_mat,*im_aux_mat;
    float *re_aux_mat2,*im_aux_mat2;
    float *Zeros;
    float twoPi2;
    float trace_reaux, trace_imaux;
    FILE *ratefile;
    N=non->singles;
    nn2=non->singles*non->singles;
    twoPi2=twoPi*twoPi;

    rate_response=(float *)calloc(non->tmax,sizeof(float));
    abs_rate_response=(float *)calloc(non->tmax,sizeof(float));
    re_aux_mat=(float *)calloc(nn2,sizeof(float));
    im_aux_mat=(float *)calloc(nn2,sizeof(float));
    re_aux_mat2=(float *)calloc(nn2,sizeof(float));
    im_aux_mat2=(float *)calloc(nn2,sizeof(float));
    Zeros=(float *)calloc(nn2,sizeof(float));
  
    ratefile=fopen("RateFile.dat","w");
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
	            abs_rate_response[t1]=sqrt(trace_reaux*trace_reaux
                                    +trace_imaux*trace_imaux)*twoPi2;
                    fprintf(ratefile,"%f %f\n",t1*non->deltat,rate_response[t1]);
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

    free(rate_response);
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
      float C;
      int i,j;
      float kBT=non->temperature*k_B; /* Kelvin to cm-1 */                     
  
      qc_rate_matrix=(float *)calloc(segments*segments,sizeof(float));
      qc=(float *)calloc(segments*segments,sizeof(float));                     
      /* Find quantum correction factors */                                    
      for (i=0;i<segments;i++){
          for (j=0;j<segments;j++){                                            
              if (i!=j){
                      /* Quantum correction factor from D.W. Oxtoby. */
                      /* Annu. Rev. Phys. Chem., 32(1):77–101, (1981).*/       
                      C=2/(1+exp((E[i]-E[j])/kBT));                            
                      qc_rate_matrix[i*segments+j]=rate_matrix[i*segments+j]*C;
                      qc_rate_matrix[j*segments+j]-=rate_matrix[i*segments+j]*C;
                  qc[i*segments+j]=C;
              }       
              else{   
                  qc[i*segments+j]=0;
              }       
          }
      }
  
      /* Write the quantum corrected rate matrix. */
      write_matrix_to_file("QC_RateMatrix.dat",qc_rate_matrix,segments);
      /* Write the applied quantum correction factors. */
      write_matrix_to_file("QC.dat",qc,segments);                              
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

/* This function will create a density matrix where every term is weighted with a Boltzmann weight */
void density_matrix(float *density_matrix, float *Hamiltonian_i,t_non *non,int segments){
    int index,N;
    float *H,*e;
    double *c2;
    double *cnr;
    double *matrix;

    N=non->singles;
    H=(float *)calloc(N*N,sizeof(float));
    e=(float *)calloc(N,sizeof(float));
    c2=(double *)calloc(N,sizeof(double));
    cnr=(double *)calloc(N*N,sizeof(double));
    matrix=(double *)calloc(N*N,sizeof(double));

    int a,b,c;
    double kBT=(double) non->temperature*k_B; /* Kelvin to cm-1 */
    double *Q,iQ;
 
    clearvec(density_matrix,N*N);

    Q=(double *)calloc(segments,sizeof(double));  

    /* Build Hamiltonian */
    for (a=0;a<N;a++){
        H[a+N*a]=Hamiltonian_i[a+N*a-(a*(a+1))/2]; /* Diagonal */
        for (b=a+1;b<N;b++){
            H[a+N*b]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
            H[b+N*a]=Hamiltonian_i[b+N*a-(a*(a+1))/2];
        }
    }
    /* Find eigenvalues and eigenvectors */
    diagonalizeLPD(H,e,N); 
 
    /* Exponentiate [U=exp(-H/kBT)] */
    for (a=0;a<N;a++){
        if (non->temperature==0){
            printf("Temperature is 0, the equilirbium density matrix will be nan,we suggestion to use a low non-zero temperature instead");
            exit(0);
        }

        c2[a]=exp(-((double) (e[a]-e[N-1]))/kBT);
        /* Apply strict high temperature limit when T>100000 */
        if (non->temperature>100000){
	        c2[a]=1.0;
        }
    }

    /* Transform back to site basis */ 
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            cnr[b+a*N]+=((double) H[b+a*N])*c2[b];
        }
    }  
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            for (c=0;c<N;c++){
                matrix[a+c*N]+=H[b+a*N]*cnr[b+c*N];
            }
        }
    }
  
    /* Find the partition function for each segment */
    for (a=0;a<N;a++){
        Q[non->psites[a]]+=matrix[a+a*N];
    }
    /* Re-normalize */
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
    	    density_matrix[a+b*N]=(float) (matrix[a+b*N]/Q[non->psites[a]]);
        }
    }      

    free(H);
    free(c2);
    free(e);
    free(cnr);
    free(Q);
    return;
}

void average_density_matrix(float *ave_den_mat,t_non *non){
/* Define variables and arrays */
   /* Integers */
    int ti;
    int segments;
    int samples;
    int ele;
    int my_samples;
    int N,a,b;
    float i_samples;
    /* Vectors representing time dependent states: real and imaginary part */
    float *vecr;
    float *Hamiltonian_i;
    /* File handles */
    FILE *H_traj;
    FILE *mu_traj;
    FILE *Cfile;
    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Allocating memory for the real and imaginary part of the wave function that we need to propagate */
    vecr=(float *)calloc(non->singles*non->singles,sizeof(float));
    Hamiltonian_i=(float *)calloc(non->singles*(non->singles+1)/2,sizeof(float));
    //ave_den_mat=(float *)calloc(non->singles*non->singles,sizeof(float));
    /* Initialize sample numbers */
    segments=project_dim(non);
    N=non->singles;
  
    clearvec(ave_den_mat,N*N);
    /* Initialize sample numbers */
    my_samples=determine_samples(non);

    if (non->end-non->begin<my_samples){
      my_samples=non->end-non->begin;
    }
    for (samples=non->begin;samples<non->end;samples++){
      ti=samples*non->sample; 
      read_Hamiltonian(non,Hamiltonian_i,H_traj,ti);
      /* Use the thermal equilibrium as initial state */
      density_matrix(vecr,Hamiltonian_i,non,segments);

      for (ele=0; ele<non->singles*non->singles; ele++){
          ave_den_mat[ele] +=vecr[ele]; 
      }
    }
    /* Zero the coupling between different segments for the averaged density *
     * matrix and normalize */
    i_samples=1.0/my_samples;
    for (a=0;a<non->singles;a++){
        for (b=0;b<non->singles;b++){
	    ave_den_mat[non->singles*b+a]*=i_samples;
            if (non->psites[a] != non->psites[b]){
               ave_den_mat[non->singles*a+b]=0.0;
               /* ave_den_mat[non->singles*b+a]=0.0; */
            } 
        }
    }
    free(vecr); 
    free(Hamiltonian_i); 
    return;
}

/* Matrix multiplication for different segments */
void segment_matrix_mul(float *rA,float *iA,float *rB,float *iB,
    float *rC,float *iC,int *psites,int segments,int si,int sk,int sj,int N){
    int i,j,k;
    /* Set initial values of results matrix to zero to be sure */
    clearvec(rC,N*N);
    clearvec(iC,N*N);

    int Npsj, Npsk;
    int cj, ck;
    Npsj=0;
    Npsk=0;

    float *psj, *psk;
    psj=(float *)calloc(N,sizeof(float));
    psk=(float *)calloc(N,sizeof(float));

    for (i=0;i<N;i++){
        if (psites[i]==sj){
            psj[Npsj]=i;
            Npsj++;
        }
        if (psites[i]==sk){
            psk[Npsk]=i;
            Npsk++;
        }
    }

#pragma omp parallel for
    for (i=0;i<N;i++){
        if (psites[i]==si){
            for (cj=0;cj<Npsj;cj++){
                j=psj[cj];
                for (ck=0;ck<Npsk;ck++){
                    k=psk[ck];
                    rC[i*N+j]+=rA[i*N+k]*rB[k*N+j]-iA[i*N+k]*iB[k*N+j];
                    iC[i*N+j]+=rA[i*N+k]*iB[k*N+j]+iA[i*N+k]*rB[k*N+j];
                }
            }
        }
    }

    /*
    for (i=0;i<N;i++){
        if (psites[i]==si){
            for (j=0;j<N;j++){
                if (psites[j]==sj){
                    for (k=0;k<N;k++){
                        if (psites[k]==sk){
                            rC[i*N+j]+=rA[i*N+k]*rB[k*N+j]-iA[i*N+k]*iB[k*N+j];
                            iC[i*N+j]+=rA[i*N+k]*iB[k*N+j]+iA[i*N+k]*rB[k*N+j];
                        } 
                    } 
                }
            }
        }
    }
    */
    free(psj);
    free(psk);
    return;
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

/* Integrate the rate response */
void integrate_rate_response(float *rate_response,int T,float *is13,float *isimple){
    int i;
    float simple; /* Variable for trapezium integral */
    float simp13; /* Variable for Simpsons 1/3 rule integral */
    simple=0;
    simp13=0;
    for (i=0;i<T;i++){
        if (i==0){
	    simple+=rate_response[i]/2;
	    simp13+=rate_response[i]/3;
	} else if (i%2==0){
	    simple+=rate_response[i];
            simp13+=2*rate_response[i]/3;
        } else {
	     simple+=rate_response[i];
            simp13+=4*rate_response[i]/3;
        }
    }

    /* Check for difference between integration methods */
    if (fabs(simple-simp13)/fabs(simple)>0.05){
        printf("\n");
        printf(YELLOW "Warning the timesteps may be to large for integration!\n" RESET);
        printf(YELLOW "Simple integral value: %f\n Simpson 1/3: %f\n",simple,simp13);
        printf(YELLOW "This difference is larger than 5%%.\n" RESET);
	printf(YELLOW "The trapezium rule value is used.\n\n" RESET);
    }

    /* Check for difference between initial and final value */
    if (fabs(rate_response[T-1])*50>rate_response[0]){
	    printf("\n");
            printf(YELLOW "Final value of rate response is %f %%\n",fabs(rate_response[T-1])*100/fabs(rate_response[0]));
	    printf("of the initial value. You may avearge over too\n");
	    printf("few samples (decrease the value of Samplerate) or\n");
	    printf("your chosen coherence time of %d steps, may\n",T);
	    printf("be too short for the coherence to decay.\n." RESET);
	    printf("\n");
    }

    /* Store results in variables for return */
    *isimple=simple;
    *is13=simp13;
}

/* Write a square matrix to a text file */
void write_matrix_to_file(char fname[],float *matrix,int N){
    FILE *file_handle;
    int i,j;
    file_handle=fopen(fname,"w");
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            fprintf(file_handle,"%10.14e ",matrix[i*N+j]);
        }
        fprintf(file_handle,"\n");
    }
    fclose(file_handle);
}

/* Read a square matrix from a text file */
void read_matrix_from_file(char fname[],float *matrix,int N){
    FILE *file_handle;
    int i,j;
    file_handle=fopen(fname,"r");
    if (file_handle == NULL) {
        printf("Error opening the file %s.\n",fname);
        exit(0);
    }
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            fscanf(file_handle,"%f",&matrix[i*N+j]);
        }
    }
    fclose(file_handle);
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
