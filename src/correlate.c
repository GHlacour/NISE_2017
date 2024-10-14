#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include <complex.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "read_trajectory.h"

// Function to calculate the mean of a signal
void subtractMean(float* signal, int N) {
    float sum;
    /* First subtract a reasonable guess from all numbers */
    sum=(signal[0]+signal[N-1])/2.0;
    for (int i = 0; i < N; i++) {
        signal[i]-=sum;
    }

    /* Subtract the proper average */
    sum=0.0;
    for (int i = 0; i < N; i++) {
        sum += signal[i];
    }
    sum=sum/N;
    for (int i = 0; i < N; i++) {
        signal[i]-=sum;
    }
}


// Function to calculate correlation function using FFTW
void calculateCorrelation(float *input1, float *input2, float *output, float *SD, int N) {
    // Create FFTW plans
    fftwf_plan plan1, plan2, plan3;

    // Allocate memory for FFTW input and output arrays
    fftwf_complex *fft_input1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * N);
    fftwf_complex *fft_input2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * N);
    fftwf_complex *fft_output = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * N);

    // Initialize FFTW plans
    plan1 = fftwf_plan_dft_r2c_1d(N, input1, fft_input1, FFTW_ESTIMATE);
    plan2 = fftwf_plan_dft_r2c_1d(N, input2, fft_input2, FFTW_ESTIMATE);
    plan3 = fftwf_plan_dft_c2r_1d(N, fft_output, output, FFTW_ESTIMATE);

    // Execute FFT for input1 and input2
    fftwf_execute(plan1);
    fftwf_execute(plan2);

    // Calculate the complex conjugate of fft_input2
    for (int i = 0; i < N; i++) {
        fft_input2[i][1] = -fft_input2[i][1];
    }

    // Multiply the two complex spectra element-wise
    // And store real part in Spectral Density
    for (int i = 0; i < N; i++) {
        fft_output[i][0] = (fft_input1[i][0] * fft_input2[i][0]) - (fft_input1[i][1] * fft_input2[i][1]);
        fft_output[i][1] = (fft_input1[i][0] * fft_input2[i][1]) + (fft_input1[i][1] * fft_input2[i][0]);
        SD[i]=fft_output[i][0]/N;
    }

    // Execute inverse FFT to get the correlation function
    fftwf_execute(plan3);

    // Normalize the correlation function
    float normalization = 1.0 / N;
    for (int i = 0; i < N; i++) {
        output[i] *= normalization;
    }

    // Free allocated memory and destroy FFTW plans
    fftwf_destroy_plan(plan1);
    fftwf_destroy_plan(plan2);
    fftwf_destroy_plan(plan3);
    fftwf_free(fft_input1);
    fftwf_free(fft_input2);
    fftwf_free(fft_output);
}

/* Control calculation of correlation functions */
void calc_Correlation(t_non *non){
    float *corr_matrix;
    float *c3_matrix;
    float *c4_matrix;
    float *traj1;
    float *traj2;
    float *traj3;
    float *corr;
    float *Hamil_i_e;
    float *SD,*SD_matrix;
    int TT,T,N,nn2;
    int a,b; /* Indices for chromophores */
    int ti;

    /* File handles */
    FILE *H_traj,*mu_traj;
    FILE *C_traj;
    FILE *outone,*log;
    FILE *outSD;
    FILE *Cfile;

    TT=non->length; /* Total timesteps */
    T=non->tmax1; /* Length of correlation functions */
    N=non->singles; /* Number of chromophores */
    nn2=N*(N+1)/2;
    traj1=(float *)calloc(TT,sizeof(float));
    traj2=(float *)calloc(TT,sizeof(float));
    traj3=(float *)calloc(TT,sizeof(float));
    corr=(float *)calloc(TT,sizeof(float));
    SD=(float *)calloc(TT,sizeof(float));
    corr_matrix=(float *)calloc(nn2*T,sizeof(float));
    SD_matrix=(float *)calloc(N*TT,sizeof(float));
    c3_matrix=(float *)calloc(N*T,sizeof(float));
    c4_matrix=(float *)calloc(N*T,sizeof(float));
    Hamil_i_e=(float *)calloc(nn2,sizeof(float));

    /* Open Trajectory files */
    open_files(non,&H_traj,&mu_traj,&Cfile);

    /* Loop over pairs */
    for (a=0;a<N;a++){
	  for (b=a;b<N;b++){
        for (ti=0;ti<TT;ti++){
          /* Read Hamiltonian */
		  read_Hamiltonian(non,Hamil_i_e,H_traj,ti);
		  traj1[ti]=Hamil_i_e[Sindex(a,a,N)];
		  traj2[ti]=Hamil_i_e[Sindex(b,b,N)];
	    }
	    /* Do WK theorm */
	    subtractMean(traj1,TT);
	    subtractMean(traj2,TT);
        calculateCorrelation(traj1,traj2,corr,SD,TT);
	    /* Store in matrix */
	    for (ti=0;ti<T;ti++){
		//printf("%f ",traj1[ti]);
		  corr_matrix[Sindex(a,b,N)*T+ti]=corr[ti]/TT;
        }
        if (a==b) {
          for (ti=0;ti<TT;ti++){
            SD_matrix[a*TT+ti]=SD[ti]/TT;
	      }
        }
	    /* Higher order correlations */
	    if (a==b){
          for (ti=0;ti<TT;ti++){
		    traj3[ti]=traj1[ti]*traj1[ti];
		  }
		  subtractMean(traj3,TT);
		  /* Third order */
		  calculateCorrelation(traj3,traj1,corr,SD,TT);
		  /* Store in matrix */
          for (ti=0;ti<T;ti++){
            c3_matrix[a*T+ti]=corr[ti]/TT;
		  }
		  /* Fourth order */
		  calculateCorrelation(traj3,traj3,corr,SD,TT);
          /* Store in matrix */
          for (ti=0;ti<T;ti++){
            c4_matrix[a*T+ti]=corr[ti]/TT;
          }
        }
	    /* Skip cross correlations if in autocorrelate mode */
	    if (string_in_array(non->technique,
          (char*[]){"Autocorrelation"},1)){
          break;
	    }
	  }
    }

    /* Save to file */
    outone=fopen("CorrelationMatrix.dat","w");
    for (ti=0;ti<T;ti++){
      fprintf(outone,"%f ",ti*non->deltat);
	  /* Loop through pairs */
	  for (a=0;a<N;a++){
        for (b=a;b<N;b++){
          fprintf(outone,"%e ",corr_matrix[Sindex(a,b,N)*T+ti]);
		  /* Skip cross correlations if in autocorrelate mode */
          if (string_in_array(non->technique,
            (char*[]){"Autocorrelation"},1)){
            break;
          }
	    }
	  }
      fprintf(outone,"\n");
    }
    fclose(outone);

    /* Save Spectral Density to file */
    outone=fopen("SpectralDensity.dat","w");
    for (ti=0;ti<TT;ti++){
      fprintf(outone,"%f ",ti/non->deltat/icm2ifs/twoPi/TT);
	  /* Loop through pairs */
	  for (a=0;a<N;a++){
        fprintf(outone,"%e ",SD_matrix[a*TT+ti]);
      }
      fprintf(outone,"\n");
    }
    fclose(outone);

    outone=fopen("Skewness.dat","w");
    for (ti=0;ti<T;ti++){
        fprintf(outone,"%f ",ti*non->deltat);
        /* Loop through pairs */
        for (a=0;a<N;a++){
            fprintf(outone,"%e ",c3_matrix[a*T+ti]);
        }
        fprintf(outone,"\n");
    }
    fclose(outone);

    outone=fopen("Kurtosis.dat","w");
    for (ti=0;ti<T;ti++){
        fprintf(outone,"%f ",ti*non->deltat);
        /* Loop through pairs */
        for (a=0;a<N;a++){
            fprintf(outone,"%e ",c4_matrix[a*T+ti]);
        }
        fprintf(outone,"\n");
    }
    fclose(outone);

    /* Close Trajectory Files */
    fclose(mu_traj),fclose(H_traj);
    free(traj1),free(traj2),free(traj3),free(corr),free(corr_matrix);
    free(SD_matrix);
    free(SD);
    free(c3_matrix),free(c4_matrix);
    free(Hamil_i_e);
}
