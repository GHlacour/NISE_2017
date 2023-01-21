#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "randomlib.h"
#include "util/asprintf.h"

// Subroutines for nonadiabatic code

// Allocate 2D memory blocks
void** calloc2D(size_t nRows, size_t nCols, size_t size, size_t sizeP) {
    void** result = malloc(nRows * sizeP);
    void* data = calloc(nRows * nCols, size);
    for(int i = 0; i < nRows; i++) {
        result[i] = (unsigned char*) data + size * i * nCols;
    }

    return result;
}

// Free 2D memory blocks
void free2D(void** arr) {
    free(arr[0]);
    free(arr);
}

// Copy a vector
void copyvec(float* a, float* b, int N) {
    int i;
    for (i = 0; i < N; i++) b[i] = a[i];
}

// Set all elements of a vector to zero
void clearvec(float* a, int N) {
    int i;
    for (i = 0; i < N; i++) a[i] = 0;
}

// Construct unit matrix
void unitmat(float *a,int N){
    int i,j;
    for (i = 0; i < N; i++){
	for (j = 0; j < N ; j++){
            a[i+N*j]=0;
    	    if (i==j) a[i+N*j]=1;
	}
    }
}

// Multiply a complex diagonal matrix on a complex vector
void vector_on_vector(float *rr,float *ir,float *vr,float *vi,int N){
    int a;
    float re,im;
    for (a=0;a<N;a++){
	re=rr[a]*vr[a]-ir[a]*vi[a];
	im=ir[a]*vr[a]+rr[a]*vi[a];
	vr[a]=re;
	vi[a]=im;
    }
}

// Multiply a real matrix on a complex vector (vr,vi)
void matrix_on_vector(float *c,float *vr,float *vi,int N){
    float *xr;
    float *xi;
    int a,b;
    xr = (float *)calloc(N * N, sizeof(float));
    xi = (float *)calloc(N * N, sizeof(float));
    // Multiply
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            xr[a]+=c[a+b*N]*vr[b];
	    xi[a]+=c[a+b*N]*vi[b];
	}
    }
    // Copy back
    copyvec(xr,vr,N);
    copyvec(xi,vi,N);
    free(xr);
    free(xi);
}

// Multiply transpose of a real matrix on a complex vector (vr,vi)
void trans_matrix_on_vector(float *c,float *vr,float *vi,int N){
    float *xr;
    float *xi;
    int a,b;
    xr = (float *)calloc(N * N, sizeof(float));
    xi = (float *)calloc(N * N, sizeof(float));
    // Multiply
    for (a=0;a<N;a++){
        for (b=0;b<N;b++){
            xr[a]+=c[b+a*N]*vr[b];
            xi[a]+=c[b+a*N]*vi[b];
        }
    }
    // Copy back
    copyvec(xr,vr,N);
    copyvec(xi,vi,N);
    free(xr);
    free(xi);
}

/**
 * Method that logs a message, in which the message can be formatted like printf accepts.
 */
void log_item(char* msgFormat, ...) {
    // Parse parameters
    va_list args;
    va_start(args, msgFormat);

    // Write to log
    FILE* log = fopen("NISE.log", "a");
    if (log == NULL) {
        printf("Could not open log file!");
        exit(1);
    }
    vfprintf(log, msgFormat, args);
    fclose(log);

    va_end(args);
}

// Set time and write to screen
time_t set_time(time_t t0) {
    return log_time(t0, stdout);
}

// Set time and write to log file
time_t log_time(time_t t0, FILE* log) {
    time_t t1;
    time(&t1);

    char* text = time_diff(t0, t1);
    // fprintf(log, text);
    free(text);
    return t1;
}

// Forms a string with the time difference between the given times
char* time_diff(time_t t0, time_t t1) {
    int s = difftime(t1, t0);
    int h = s / 3600;
    s = s % 3600;
    int m = s / 60;
    s = s % 60;

    char* text;
    asprintf(&text, "Time spent: %dh %dmin %ds\n", h, m, s);
    return text;
}

// Forms a string with the times for MPI_Wtime
char* MPI_time(double t0) {
    int ms =t0*1000; // Convert to milliseconds
    int h = ms / 3600000;
    ms = ms % 3600000;
    int m = ms / 60000;
    ms = ms % 60000;
    int s = ms / 1000;
    ms = ms % 1000;

    char* text;
    asprintf(&text, " %dh %dmin %ds %dms\n", h, m, s, ms);
    return text;
}

/* INDEXING FOR ELECTRONIC STATES */
int Eindex(int a, int b, int N) {
    int ind;
    if (a > b) {
        ind = a  + (b - 1) * (N + N - b - 2) / 2;
    }
    else {
        ind = b  + (a - 1) * (N + N - a - 2) / 2;
    }
//    printf("%d %d %d %d\n",a,b,ind,N);
    return ind;
}

/* Read Hamiltonian */
int read_He(t_non* non, float* He, FILE* FH, int pos) {
    int i,j, N, control, t;
    float* H; // Help Hamiltonain
    float* R,*R2; // Distance if needed
    float* mu; // Dipole moment if needed
    float A; // Conversion factor for TDC coupling
    float q1,q2,J; // Charges for EDC coupling
    float Rx,Ry,Rz,dist,idist,idist3; // Variables for TDC
    float m1x,m1y,m1z,m2x,m2y,m2z; // Variables for TDC
    float f1,f2; // Variables for TDC
    FILE *pos_traj;
    FILE *dip_traj;
    float box[3];
    FILE *pbc_traj;

    /* If pbc is needed read from file */
    box[0]=box[1]=box[2]=0.0;
    if (strlen(non->pbcFName)>0){
        pbc_traj=fopen(non->pbcFName,"rb");
        if (pbc_traj==NULL){
            printf(RED "Periodic Boundary Condition file not found!\n" RESET);
            exit(1);
        }
        fseek(pbc_traj, pos * (sizeof(int) + sizeof(float) * (3)),SEEK_SET);
        fread(&t, sizeof(int), 1, pbc_traj);
        fread(box, sizeof(float), 3, pbc_traj);
        fclose(pbc_traj);
    }

    /* Read only diagonal part */
    if ((!strcmp(non->hamiltonian, "Coupling") && pos >= 0) || (!strcmp(non->hamiltonian, "TransitionDipole")) || (!strcmp(non->hamiltonian, "ExtendedDipole")) ) {
        H = (float *)calloc(non->singles, sizeof(float));
        /* Find position */
        fseek(FH, pos * (sizeof(int) + sizeof(float) * (non->singles)),SEEK_SET);
        /* Read time */
        control = fread(&t, sizeof(int), 1, FH); /* control=1; */
        if (control > non->length + non->begin * non->sample) {
            printf("Control character error in Hamiltonian file!\n");
            printf("Control character is '%d'.\n", control);
            printf("Exceeding max value of '%d'.\n", non->length + non->begin * non->sample);
            printf("Check that the numbers of singles and doubles is correct!\n");
            exit(-1);
        }
        /* Read single excitation Hamiltonian */
        fread(H, sizeof(float), non->singles, FH);

        /* Shift center and update full Hamiltonian */
        for (i = 0; i < non->singles; i++) {
            He[i * non->singles + i - (i * (i + 1)) / 2] = H[i] - non->shifte;
        }
        free(H);
    }
    else {
        /* Read Full Hamiltonian */
        if (pos == -1) { pos = 0; }
        N = non->singles * (non->singles + 1) / 2;
        /* Find position */
        fseek(FH, pos * (sizeof(int) + sizeof(float) * (non->singles * (non->singles + 1) / 2 + non->doubles * (non->doubles + 1) / 2)),SEEK_SET);
        /* Read time */
        control = fread(&t, sizeof(int), 1, FH); /* control=1; */
        if (control > non->length + non->begin * non->sample) {
            printf("Control character error in Hamiltonian file!\n");
            printf("Control character is '%d'.\n", control);
            printf("Exceeding max value of '%d'.\n", non->length + non->begin * non->sample);
            printf("Check that the numbers of singles and doubles is correct!\n");
            exit(-1);
        }
        /* Read single excitation Hamiltonian */
        fread(He, sizeof(float), N, FH);
        /* Shift center */
        for (i = 0; i < non->singles; i++) {
            He[i * non->singles + i - (i * (i + 1)) / 2] -= non->shifte;
        }
    }
    /* Find the couplings from the TDC 'on the fly' scheme  */
    A=5034.11861687; /* Convert to cm-1 from Deb**2/Ang**3 */
    if ((!strcmp(non->hamiltonian, "TransitionDipole"))) {
        R = (float *)calloc(3*non->singles, sizeof(float));
        mu = (float *)calloc(3*non->singles, sizeof(float));
	/* Read in positions */
        pos_traj=fopen(non->positionFName,"rb");
	if (pos_traj==NULL){
    	    printf(RED "Atom position file for TDC on the fly not found!\n" RESET);
            exit(1);
        }
        read_mue(non,R,pos_traj,pos,0);
        read_mue(non,R+non->singles,pos_traj,pos,1);
        read_mue(non,R+2*non->singles,pos_traj,pos,2);
        /* Read in dipoles */
        dip_traj=fopen(non->dipoleFName,"rb");
        read_mue(non,mu,dip_traj,pos,0);
        read_mue(non,mu+non->singles,dip_traj,pos,1);
        read_mue(non,mu+2*non->singles,dip_traj,pos,2);
        /* Calculate the couplings according to TDC */
        for (i = 0; i < non->singles; i++) {
            m1x=mu[i];
            m1y=mu[non->singles+i];
            m1z=mu[2*non->singles+i];
            for (j = i+1; j < non->singles; j++) {
                Rx=pbc1(R[i]-R[j],0,box);
                Ry=pbc1(R[non->singles+i]-R[non->singles+j],1,box);
                Rz=pbc1(R[2*non->singles+i]-R[2*non->singles+j],2,box);
                m2x=mu[j];
                m2y=mu[non->singles+j];
                m2z=mu[2*non->singles+j];
                dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                idist=1.0/dist;
                idist3=idist*idist*idist;
                f1=m1x*m2x+m1y*m2y+m1z*m2z;
                f2=-3*(m1x*Rx+m1y*Ry+m1z*Rz)*(m2x*Rx+m2y*Ry+m2z*Rz)*idist*idist;
                He[Sindex(i,j,non->singles)] = A*(f1+f2)*idist3;
	    }
        }
       free(R);
       free(mu);
       fclose(pos_traj);
       fclose(dip_traj);
    }
    /* Find the couplings from the EDC 'on the fly' scheme  */
    A=5034.11861687; /* Convert to cm-1 from Deb**2/Ang**3 */
    if ((!strcmp(non->hamiltonian, "ExtendedDipole"))) {
        R = (float *)calloc(3*non->singles, sizeof(float));
        R2 = (float *)calloc(3*non->singles, sizeof(float));
        mu = (float *)calloc(3*non->singles, sizeof(float));
        /* Read in positions */
        pos_traj=fopen(non->positionFName,"rb");
	if (pos_traj==NULL){
            printf(RED "Atom position file for TDC on the fly not found!\n" RESET);
            exit(1);
        }
        read_mue(non,R,pos_traj,2*pos,0);
        read_mue(non,R+non->singles,pos_traj,2*pos,1);
        read_mue(non,R+2*non->singles,pos_traj,2*pos,2);
        read_mue(non,R2,pos_traj,2*pos+1,0);
        read_mue(non,R2+non->singles,pos_traj,2*pos+1,1);
        read_mue(non,R2+2*non->singles,pos_traj,2*pos+1,2);
        /* Read in dipoles */
        dip_traj=fopen(non->dipoleFName,"rb");
        read_mue(non,mu,dip_traj,pos,0);
        read_mue(non,mu+non->singles,dip_traj,pos,1);
        read_mue(non,mu+2*non->singles,dip_traj,pos,2);
        /* Calculate the couplings according to EDC */
        for (i = 0; i < non->singles; i++) {
            /* Find q1 */
            Rx=pbc1(R2[i]-R[i],0,box);
            Ry=pbc1(R2[non->singles+i]-R[non->singles+i],1,box);
            Rz=pbc1(R2[2*non->singles+i]-R[2*non->singles+i],2,box);
            m1x=mu[i];
            m1y=mu[non->singles+i];
            m1z=mu[2*non->singles+i];
            dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
            q1=sqrt(m1x*m1x+m1y*m1y+m1z*m1z)/dist;
            for (j = i+1; j < non->singles; j++) {
            	/* Find q2 */
            	Rx=pbc1(R2[j]-R[j],0,box);
            	Ry=pbc1(R2[non->singles+j]-R[non->singles+j],1,box);
            	Rz=pbc1(R2[2*non->singles+j]-R[2*non->singles+j],2,box);
            	m2x=mu[j];
            	m2y=mu[non->singles+j];
            	m2z=mu[2*non->singles+j];
            	dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
            	q2=sqrt(m1x*m1x+m1y*m1y+m1z*m1z)/dist;
                /* ++ */
                Rx=pbc1(R[i]-R[j],0,box);
                Ry=pbc1(R[non->singles+i]-R[non->singles+j],1,box);
                Rz=pbc1(R[2*non->singles+i]-R[2*non->singles+j],2,box);
                dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                idist=1.0/dist;
                J=idist;
                /* +- */
                Rx=pbc1(R[i]-R2[j],0,box);
                Ry=pbc1(R[non->singles+i]-R2[non->singles+j],1,box);
                Rz=pbc1(R[2*non->singles+i]-R2[2*non->singles+j],2,box);
                dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                idist=1.0/dist;
                J=J-idist;
                /* -+ */
                Rx=pbc1(R2[i]-R[j],0,box);
                Ry=pbc1(R2[non->singles+i]-R[non->singles+j],1,box);
                Rz=pbc1(R2[2*non->singles+i]-R[2*non->singles+j],2,box);
                dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                idist=1.0/dist;
                J=J-idist;
                /* -- */
                Rx=pbc1(R2[i]-R2[j],0,box);
                Ry=pbc1(R2[non->singles+i]-R2[non->singles+j],1,box);
                Rz=pbc1(R2[2*non->singles+i]-R2[2*non->singles+j],2,box);
                dist=sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                idist=1.0/dist;
                J=J+idist;
                He[Sindex(i,j,non->singles)] = A*q1*q2*J;
            }
        }
       free(R);
       free(R2);
       free(mu);
       fclose(pos_traj);
       fclose(dip_traj);
    }
    // Remove couplings below coupling cut
    for (i = 0; i < non->singles; i++) {
        for (j = i+1; j < non->singles; j++) {
            if (fabs(He[Sindex(i,j,non->singles)])<non->couplingcut){
               He[Sindex(i,j,non->singles)]=0;
            }
        }
    }

    return control;
}

/* Read Diagonal Hamiltonian */
int read_Dia(t_non* non, float* He, FILE* FH, int pos) {
    int i, N, control, t;
    float* H;
    H = (float *)calloc(non->singles, sizeof(float));
    /* N=non->singles*(non->singles+1)/2; */
    /* Find position */
    fseek(FH, pos * (sizeof(int) + sizeof(float) * (non->singles)),SEEK_SET);
    /* Read time */
    control = fread(&t, sizeof(int), 1, FH); /* control=1; */
    if (control > non->length + non->begin * non->sample) {
        printf("Control character error in Hamiltonian file!\n");
        printf("Control character is '%d'.\n", control);
        printf("Exceeding max value of '%d'.\n", non->length + non->begin * non->sample);
        printf("Check that the numbers of singles and doubles is correct!\n");
        exit(-1);
    }
    /* Read single excitation Hamiltonian */
    fread(H, sizeof(float), non->singles, FH);
    /* Shift center and update full Hamiltonian */
    for (i = 0; i < non->singles; i++) {
        He[i * non->singles + i - (i * (i + 1)) / 2] = H[i] - non->shifte;
    }
    free(H);
    return control;
}

/* Read the diagonal anharmonicities */
int read_A(t_non* non, float* Anh, FILE* FH, int pos) {
    int i, N, control, t;
    N = non->singles;
    /* Find Position */
    fseek(FH, pos * (sizeof(int) + sizeof(float) * non->singles),SEEK_SET);
    /* Read time */
    control = fread(&t, sizeof(int), 1, FH); // control=1;
    /* Read single excitation Hamiltonian */
    fread(Anh, sizeof(float), N, FH);
    return control;
}

/* Read Dipole */
int read_mue(t_non* non, float* mue, FILE* FH, int pos, int x) {
    int control;
    int t;
    int N;
    control = 0;
    // Find position
    fseek(FH, pos * (sizeof(int) + sizeof(float) * (3 * non->singles + 3 * non->singles * non->doubles)) + sizeof(float)
          * x * non->singles,SEEK_SET);
    /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    // Read single excitation Dipoles
    fread(mue, sizeof(float), non->singles, FH);
    return control;
}

/* Read Raman */
int read_alpha(t_non* non, float* alpha, FILE* FH, int pos, int x) {
    int control;
    int t;
    //int N;
    control = 0;

    if (non->doubles == 0){
        //Find position in file to read from
        fseek(FH, pos * (sizeof(int) + sizeof(float) * 6 * non->singles) + sizeof(float)
          * x * non->singles,SEEK_SET); //6 due to six elements for Raman tensor
        /* Read time */
        if (fread(&t, sizeof(int), 1, FH)) control = 1;
        // Read single excitation raman tensor values
        fread(alpha, sizeof(float), non->singles, FH);
    }
    else{
        printf("Doubles not equal to zero, Raman 1D does not support double excited states");
        exit(-1);
    }
    return control;
}

/* Read Dipole for Sequence Transtions */
int read_over(t_non* non, float* over, FILE* FH, int pos, int x) {
    int control;
    int t;
    int N;
    control = 0;
    // Find position
    fseek(FH, pos * (sizeof(int) + sizeof(float) * (3 * non->singles)) + sizeof(float) * x * non->singles,SEEK_SET);
    /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    // Read single excitation Dipoles
    fread(over, sizeof(float), non->singles, FH);
    return control;
}

// Read transition dipole
void muread(t_non* non, float* leftnr, int ti, int x, FILE* mu_traj) {
    /* Read mu(ti) */
    if (read_mue(non, leftnr, mu_traj, ti, x) != 1) {
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n", ti, x);
        exit(1);
    }
    return;
}

// Read transition dipole (including option for Coupling storrage)
void mureadE(t_non* non, float* leftnr, int ti, int x, FILE* mu_traj, float* mu, float* pol) {
    if (!strcmp(non->hamiltonian, "Coupling")) {
        copyvec(mu + non->singles * x, leftnr, non->singles);
    }
    else {
        /* Read mu(ti) */
        if (read_mue(non, leftnr, mu_traj, ti, x) != 1) {
            printf("Dipole trajectory file to short, could not fill buffer!!!\n");
            printf("ITIME %d %d\n", ti, x);
            exit(1);
        }
    }
    return;
}

/* Read Cluster file */
int read_cluster(t_non* non, int pos, int* cl, FILE* FH) {
    int control;
    int t;
    int N;
    control = 0;
    // Find position
    fseek(FH, pos * (sizeof(int) + sizeof(int)),SEEK_SET);
    /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    // Read single excitation Dipoles
    fread(cl, sizeof(int), 1, FH);
    //  printf("%d %d %d\n",pos,t,cl);
    return control;
}




// This subroutine generates vectors with a coordinate transformation
// corresponding to a randomly selected isotropic orientation
void generateCS(float* X, float* Y, float* Z) {
    int no;
    // Generate X vector
    X[0] = RandomGaussian(0, 1);
    X[1] = RandomGaussian(0, 1);
    X[2] = RandomGaussian(0, 1);
    // Normalize
    no = 1 / sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
    X[0] = X[0] / no;
    X[1] = X[1] / no;
    X[2] = X[2] / no;
    // Generate Y vector
    Y[0] = RandomGaussian(0, 1);
    Y[1] = RandomGaussian(0, 1);
    Y[2] = RandomGaussian(0, 1);
    // Make it orthogonal to X
    no = X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2];
    Y[0] = Y[0] - no * X[0];
    Y[1] = Y[1] - no * X[1];
    Y[2] = Y[2] - no * X[2];
    // Normalize Y
    no = 1 / sqrt(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2]);
    Y[0] = Y[0] / no;
    Y[1] = Y[1] / no;
    Y[2] = Y[2] / no;
    // Generate Z by taking cross product
    Z[0] = X[1] * Y[2] - X[2] * Y[1];
    Z[1] = -X[0] * Y[2] + X[2] * Y[0];
    Z[2] = X[0] * Y[1] - X[1] * Y[0];
    // Normalize Z
    no = 1 / sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
    Z[0] = Z[0] / no;
    Z[1] = Z[1] / no;
    Z[2] = Z[2] / no;
    return;
}

/* Test at the start if the Hamiltonian and dipole files are sensible */
int control(t_non* non) {
    float *mu_eg, *Hamil_i_e;
    FILE *H_traj, *mu_traj;
    FILE *x_traj;
    int itime, N_samples;
    int samples;
    int nn2;

    nn2 = non->singles * (non->singles + 1) / 2;
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    mu_eg = (float *)calloc(non->singles, sizeof(float));
    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        return 1;
    }
    mu_traj = fopen(non->dipoleFName, "rb");
    if (mu_traj == NULL) {
        printf("Dipole file %s not found!\n", non->dipoleFName);
        return 1;
    }
    N_samples = (non->length - non->tmax1 - 1) / non->sample + 1;
    if (N_samples < 0) {
        printf("Insufficient data to calculate spectrum.\n");
        printf("Please, lower max times or provide longer\n");
        printf("trajectory.\n");
        return 1;
    }

    // Check first element
    // Read Hamiltonian
    if (read_He(non, Hamil_i_e, H_traj, 0) != 1) {
        printf("Failed initial control\n");
        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
        return 1;
    }
    if (read_mue(non, mu_eg, mu_traj, 0, 0) != 1) {
        printf("Failed initial control\n");
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n", 0, 0);
        return 1;
    }
    if (!strcmp(non->hamiltonian, "Coupling")) { }
    else if (!strcmp(non->hamiltonian, "TransitionDipole") || !strcmp(non->hamiltonian, "ExtendedDipole")) {
        x_traj = fopen(non->positionFName, "rb");
        if (x_traj == NULL) {
          printf("Position file not found!\n");
          return 1;
        }
        fclose(x_traj);
    }

    else {
        // Check last element
        if (read_mue(non, mu_eg, mu_traj, non->length - 1, 2) != 1) {
            printf("Dipole trajectory file to short, could not fill buffer!!!\n");
            printf("ITIME %d %d\n", non->length - 1, 2);
            return 1;
        }
        // Read Hamiltonian
        if (read_He(non, Hamil_i_e, H_traj, non->length - 1) != 1) {
            printf(RED "Failed initial control\n");
            printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
            printf("Real file length shorter than specified with Length keyword!\n" RESET);
            return 1;
        }
        // Check Hamiltonian elements
        if (Hamil_i_e[0] + non->shifte > non->max1 || Hamil_i_e[0] + non->shifte < non->min1) {
            printf(RED "Warning: Hamiltonian value %f outside expected range.\n", Hamil_i_e[0] + non->shifte);
            printf("Expected frequency range: %f to %f.\n", non->min1, non->max1);
            printf("Computation will continue, but check is the number above is realistic\n");
            printf("You may have specified wrong number of sites!\n" RESET);
            printf("---------------------------------------------------------------------\n");
        }
    }
    free(mu_eg);
    free(Hamil_i_e);
    fclose(mu_traj), fclose(H_traj);
    return 0;
}

/* Routine for autodetection of the number of singles */
int autodetect_singles(t_non* non){
    float *Hamil_i_e;
    FILE *H_traj;
    int i;
    float f;
    int samples;
    int n,nn2;
    int identified;

    identified=0;
    nn2=non->singles*(non->singles+1)/2;
    Hamil_i_e = (float *)calloc(nn2, sizeof(float));
    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        return 1;
    }
    for (n=1;n<non->singles*10;n++){
      if (!strcmp(non->hamiltonian, "Coupling")) {
         identified=-2;
      }
      if (!strcmp(non->hamiltonian, "TransitionDipole") || !strcmp(non->hamiltonian, "ExtendedDipole")){
         identified=-2;
      }
      if (!strcmp(non->hamiltonian, "Full")) {
	 fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
	 fread(&i,sizeof(int),1,H_traj);
         fseek(H_traj, 1 * (sizeof(int) + sizeof(float) * (n*(n+1)/2)),SEEK_SET);
         fread(&f,sizeof(float),1,H_traj);
	 //printf("%d %d %f\n",n,i,f);
         if (abs(i)<1000 || f==floorf(f)){
            printf("Autodetected potential singles at %d\n",n);
	    if (n==non->singles){
	       identified=-1;
	       break;
	    }
	    if (n!=non->singles){
	       identified=n;
	       break;
	    }
	 }
      }
    }
    if (identified==0){
       printf(RED "Warning: Autodetection of sites failed. You may need to increase Singles\n" RESET);
    }
    if (identified==-1){
       printf("Singles confirmed by auto detection.\n");
    }
    if (identified>0){
      printf(RED "Warning: Singles keyword may be specified incorrectly!\n");
      printf("Autodetection suggested %d singles.\n" RESET,identified);
    }
    fclose(H_traj);
    free(Hamil_i_e);
    return 0;
}

/* Multiply with double exciton dipole mu_ef on single states */
void dipole_double(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }

    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[index] += over[i] * cr[i];
        fi[index] += over[i] * ci[i];
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[index] += dipole[i] * cr[j];
            fi[index] += dipole[i] * ci[j];
            fr[index] += dipole[j] * cr[i];
            fi[index] += dipole[j] * ci[i];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on ground states */
void dipole_double_ground(t_non* non, float* dipole, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }

    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[index] = over[i] ;
        //! no fi since double excitation from ground state,
        //no i,j loop due to double excitation on one state only
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on single states */
void dipole_double_ES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < N; i++) fr[i] = 0, fi[i] = 0;

    for (i = 0; i < non->singles; i++) {
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[index] += dipole[i] * cr[j];
            fi[index] += dipole[i] * ci[j];
            fr[index] += dipole[j] * cr[i];
            fi[index] += dipole[j] * ci[i];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on double states */
void dipole_double_last(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi, float* over) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < non->singles; i++) fr[i] = 0, fi[i] = 0;
    if (non->anharmonicity != 0) {
        for (i = 0; i < non->singles; i++) {
            over[i] = sqrt2 * dipole[i];
        }
    }
    for (i = 0; i < non->singles; i++) {
        index = Sindex(i, i, non->singles);
        fr[i] += over[i] * cr[index];
        fi[i] += over[i] * ci[index];
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[j] += dipole[i] * cr[index];
            fi[j] += dipole[i] * ci[index];
            fr[i] += dipole[j] * cr[index];
            fi[i] += dipole[j] * ci[index];
        }
    }
    return;
}

/* Multiply with double exciton dipole mu_ef on double states */
void dipole_double_last_ES(t_non* non, float* dipole, float* cr, float* ci, float* fr, float* fi) {
    int N;
    int i, j, k, index;
    N = non->singles * (non->singles + 1) / 2;
    for (i = 0; i < non->singles; i++) fr[i] = 0, fi[i] = 0;
    for (i = 0; i < non->singles; i++) {
        for (j = i + 1; j < non->singles; j++) {
            index = Sindex(i, j, non->singles);
            fr[j] += dipole[i] * cr[index];
            fi[j] += dipole[i] * ci[index];
            fr[i] += dipole[j] * cr[index];
            fi[i] += dipole[j] * ci[index];
        }
    }
    return;
}

/* Return the distance between two locations squared */
float distance(float *rf,float *ri,int a,int b,int N,float box){
  float d,r;
  int x;
  d=0;
  for (x=0;x<3;x++){
    r=rf[3*a+x]-ri[3*b+x];
    if (r>box/2) r=r-box;
    if (r<-box/2) r=r+box;
    d+=r*r;
  }
  return d;
}

/* Return the distance between two locations along a direction x */
float distance_x(float *rf,float *ri,int a,int b,int N,float box,int x){
  float r;
  r=rf[3*a+x]-ri[3*b+x];
  if (r>box/2) r=r-box;
  if (r<-box/2) r=r+box;
  return r;
}

/* Return the distance between two locations squared */
float distance3(float *rf,float *ri,int a,int b,int N,float *box){
  float d,r;
  int x;
  d=0;
  for (x=0;x<3;x++){
    r=rf[3*a+x]-ri[3*b+x];
    if (r>box[x]/2) r=r-box[x];
    if (r<-box[x]/2) r=r+box[x];
    d+=r*r;
  }
  return d;
}

/* Return the distance between two locations along a direction x */
float distance3_x(float *rf,float *ri,int a,int b,int N,float *box,int x){
  float r;
  r=rf[3*a+x]-ri[3*b+x];
  if (r>box[x]/2) r=r-box[x];
  if (r<-box[x]/2) r=r+box[x];
  return r;
}

float pbc1(float r, int x, float *box){
  // Correct for pbc if active
  if (box[0]>0.0){
     if (r>box[x]/2) r=r-box[x];
     if (r<-box[x]/2) r=r+box[x];
  }
  return r;
}
