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


// Propagate using standard matrix exponential
void propagate_vec_DIA(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int sign) {
    float f;
    int index, N;
    float *H, *re_U, *im_U, *e;
    float *cnr, *cni;
    float *crr, *cri;
    float re, im;
    int a, b, c;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N * N, sizeof(float));
    cni = (float *)calloc(N * N, sizeof(float));
    crr = (float *)calloc(N * N, sizeof(float));
    cri = (float *)calloc(N * N, sizeof(float));
    // Build Hamiltonian
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
    // Exponentiate [U=exp(-i/h H dt)]
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    // Transform to site basis
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] += H[b + a * N] * re_U[b], cni[b + a * N] += H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                crr[a + c * N] += H[b + a * N] * cnr[b + c * N], cri[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    // The one exciton propagator has been calculated

    //  elements=0;
    for (a = 0; a < N; a++) {
        cnr[a] = 0, cni[a] = 0;
        for (b = 0; b < N; b++) {
            //      if ((crr[a+b*N]*crr[a+b*N]+cri[a+b*N]*cri[a+b*N])>non->thres){
            //        elements++;
            cnr[a] += crr[a + b * N] * cr[b] - cri[a + b * N] * ci[b];
            cni[a] += crr[a + b * N] * ci[b] + cri[a + b * N] * cr[b];
            //      }
        }
    }

    for (a = 0; a < N; a++) {
        cr[a] = cnr[a], ci[a] = cni[a];
    }


    free(cnr), free(cni), free(re_U), free(im_U), free(H), free(e);
    free(crr), free(cri);
    return;
}

// Do the propagation for t2 using the matrix exponent and using a single diagonalization
// for all vectors
void propagate_t2_DIA(t_non *non,float *Hamiltonian_i,float *cr,float *ci,float **vr,float **vi,int sign){
    float f;
    int index, N;
    float *H, *re_U, *im_U, *e;
    float re, im;
    int a, b, c;
    int t1;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
   
    // Build Hamiltonian
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
    // Exponentiate [U=exp(-i/h H dt)]
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    // Do the single t1 independent vector
    // Transfer to eigen basis
    matrix_on_vector(H,cr,ci,N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);
    // Transfer back to site basis
    trans_matrix_on_vector(H,cr,ci,N);

    // Do all the t1 dependent vectors
#pragma omp parallel for shared(non,re_U,im_U,H,vr,vi) schedule(static,1)
    for (t1=0;t1<non->tmax1;t1++){
        // Transfer to eigen basis
        matrix_on_vector(H,vr[t1],vi[t1],N);
        // Multiply with matrix exponent
        vector_on_vector(re_U,im_U,vr[t1],vi[t1],N);
        // Transfer back to site basis
        trans_matrix_on_vector(H,vr[t1],vi[t1],N);
    }

    free(re_U), free(im_U), free(H), free(e);
    return;
}

// Propagate using matrix exponential sparce
int propagate_vec_DIA_S(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int sign) {
    int elements;
    float f;
    int index, N, N2;
    float *H, *re_U, *im_U, *e;
    float *cnr, *cni;
    float *crr, *cri;
    float re, im;
    int a, b, c;

    N = non->singles;
    N2 = N * N;
    f = non->deltat * icm2ifs * twoPi * sign;
    H = (float *)calloc(N2, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N2, sizeof(float));
    cni = (float *)calloc(N2, sizeof(float));
    crr = (float *)calloc(N2, sizeof(float));
    cri = (float *)calloc(N2, sizeof(float));

    // Build Hamiltonian
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }
    diagonalizeLPD(H, e, N);
    // Exponentiate [U=exp(-i/h H dt)]
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    // Transform to site basis
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] = H[b + a * N] * re_U[b], cni[b + a * N] = H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (c = 0; c < N; c++) {
            crr[a + c * N] = 0, cri[a + c * N] = 0;
            for (b = 0; b < N; b++) {
                crr[a + c * N] += H[b + a * N] * cnr[b + c * N], cri[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    // The one exciton propagator has been calculated

    elements = 0;
    for (a = 0; a < N; a++) {
        cnr[a] = 0, cni[a] = 0;
        for (b = 0; b < N; b++) {
            if ((crr[a + b * N] * crr[a + b * N] + cri[a + b * N] * cri[a + b * N]) > non->thres) {
                elements++;
                cnr[a] += crr[a + b * N] * cr[b] - cri[a + b * N] * ci[b];
                cni[a] += crr[a + b * N] * ci[b] + cri[a + b * N] * cr[b];
            }
        }
    }

    for (a = 0; a < N; a++) {
        cr[a] = cnr[a], ci[a] = cni[a];
    }

    free(crr), free(cri);
    free(cnr), free(cni), free(re_U), free(im_U), free(H), free(e);

    return elements;
}

// Propagate using diagonal vs. coupling sparce algorithm
void propagate_vec_coupling_S(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, int sign) {
    float f;
    int index, N;
    float *H1, *H0, *re_U, *im_U;
    int *col, *row;
    float *ocr, *oci;
    int a, b, c;
    //  float *Norm;
    float J;
    float cr1, cr2, ci1, ci2;
    float co, si;
    int i, k, kmax;

    N = non->singles;
    f = non->deltat * icm2ifs * twoPi * sign / m;
    H0 = (float *)calloc(N, sizeof(float));
    H1 = (float *)calloc(N * N, sizeof(float));
    col = (int *)calloc(N * N / 2, sizeof(int));
    row = (int *)calloc(N * N / 2, sizeof(int));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    ocr = (float *)calloc(N, sizeof(float));
    oci = (float *)calloc(N, sizeof(float));
    //  Norm=(float *)calloc(N,sizeof(float));

    // Build Hamiltonians H0 (diagonal) and H1 (coupling)
    k = 0;
    for (a = 0; a < N; a++) {
        H0[a] = Hamiltonian_i[Sindex(a, a, N)]; // Diagonal
        for (b = a + 1; b < N; b++) {
            index = Sindex(a, b, N);
            if (fabs(Hamiltonian_i[index]) > non->couplingcut) {
                H1[k] = Hamiltonian_i[index];
                col[k] = a, row[k] = b;
                k++;
            }
        }
    }
    kmax = k;

    // Exponentiate diagonal [U=exp(-i/2h H0 dt)]
    for (a = 0; a < N; a++) {
        re_U[a] = cos(0.5 * H0[a] * f);
        im_U[a] = -sin(0.5 * H0[a] * f);
    }

    for (i = 0; i < m; i++) {
        // Multiply on vector first time
        for (a = 0; a < N; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = ci[a] * re_U[a] + cr[a] * im_U[a];
        }


        // Account for couplings
        for (k = 0; k < kmax; k++) {
            a = col[k];
            b = row[k];
            J = H1[k];
            J = J * f;
            si = -sin(J);
            co = sqrt(1 - si * si);
            cr1 = co * ocr[a] - si * oci[b];
            ci1 = co * oci[a] + si * ocr[b];
            cr2 = co * ocr[b] - si * oci[a];
            ci2 = co * oci[b] + si * ocr[a];
            ocr[a] = cr1, oci[a] = ci1, ocr[b] = cr2, oci[b] = ci2;
        }

        // Multiply on vector second time
        for (a = 0; a < N; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = oci[a] * re_U[a] + ocr[a] * im_U[a];
        }
    }

    // Move to back in original
    //  for (a=0;a<N;a++){
    //  cr[a]=ocr[a],ci[a]=oci[a];
    //    nm+=cr[a]*cr[a]+ci[a]*ci[a];
    //}
    //  printf("N3 %f\n",nm);

    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}

/* Propagate doubles using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S_doubles(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m, float* Anh) {
    int N = non->singles;
    int N2 = N * (N + 1) / 2;
    const float f = non->deltat * icm2ifs * twoPi / m;
    float si, co;
    float si2, co2;
    float* H0 = calloc(N2, sizeof(float));
    float* H1 = calloc(N * N / 2, sizeof(float));
    int* col = calloc(N * N / 2, sizeof(int));
    int* row = calloc(N * N / 2, sizeof(int));
    float* re_U = calloc(N2, sizeof(float));
    float* im_U = calloc(N2, sizeof(float));
    float* ocr = calloc(N2, sizeof(float));
    float* oci = calloc(N2, sizeof(float));

    /* Build Hamiltonians H0 (diagonal) and H1 (coupling) */
    for (int a = 0; a < N; a++) {
        const int indexa = Sindex(a, a, N);
        for (int b = a; b < N; b++) {
            int index = Sindex(a, b, N);
            H0[index] = Hamiltonian_i[indexa] + Hamiltonian_i[Sindex(b, b, N)]; // Diagonal
            if (a == b) {
                if (non->anharmonicity == 0) {
                    H0[index] -= Anh[a];
                }
                else {
                    H0[index] -= non->anharmonicity;
                }
            }
        }
    }

    /* Build Hamiltonian H1 (coupling) */
    int kmax = 0;
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            int index = b + a * ((N << 1) - a - 1) / 2; // Part of Sindex, but b > a is always true here

            if (fabsf(Hamiltonian_i[index]) > non->couplingcut) {
                H1[kmax] = Hamiltonian_i[index];
                col[kmax] = a, row[kmax] = b;
                kmax++;
            }
        }
    }

    /* Exponentiate diagonal [U=exp(-i/2h H0 dt)] */
    for (int a = 0; a < N2; a++) {
        re_U[a] = cosf(0.5f * H0[a] * f);
        im_U[a] = -sinf(0.5f * H0[a] * f);
    }

    for (int i = 0; i < m; i++) {

        /* Multiply on vector first time */
        for (int a = 0; a < N2; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = cr[a] * im_U[a] + ci[a] * re_U[a];
        }

        /* Account for couplings */
        /* Loop over couplings */
        for (int k = 0; k < kmax; k++) {
            int a = col[k];
            int b = row[k];
            int c = 0;
            float J = H1[k] * f;
            int aNsum = a * ((N << 1) - a - 1) / 2;  // copied part of Sindex
            int bNsum = b * ((N << 1) - b - 1) / 2;  // copied part of Sindex

            /* Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> */

            // c < a,b
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < a; c++) {
                int index1 = a + c * ((N << 1) - c - 1) / 2; // part of Sindex, but always a > c
                int index2 = b + c * ((N << 1) - c - 1) / 2; // part of Sindex, but always b > c
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            // c == a
            si2 = -sinf(J * sqrt2);
            co2 = sqrtf(1 - si2 * si2);
            for (; c == a; c++) {  // yeah, one iteration. but loop for consistency across all 5
                int index1 = a + c * ((N << 1) - c - 1) / 2;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co2 * ocr[index1] - si2 * oci[index2];
                float ci1 = co2 * oci[index1] + si2 * ocr[index2];
                float cr2 = co2 * ocr[index2] - si2 * oci[index1];
                float ci2 = co2 * oci[index2] + si2 * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            // a < c < b
            //si = -sinf(J);
            //co = sqrtf(1 - si * si);
            for (; c < b; c++) { // could be 0 iterations if (b == a + 1)
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            // c == b
            //si = -sinf(J * sqrt2);
            //co = sqrtf(1 - si * si);
            for (; c == b; c++) {
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co2 * ocr[index1] - si2 * oci[index2];
                float ci1 = co2 * oci[index1] + si2 * ocr[index2];
                float cr2 = co2 * ocr[index2] - si2 * oci[index1];
                float ci2 = co2 * oci[index2] + si2 * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            // c > a,b
            //si = -sinf(J);
            //co = sqrtf(1 - si * si);
            for (; c < N; c++) {
                int index1 = c + aNsum;
                int index2 = c + bNsum;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }
            
        }

        /* Multiply on vector second time */
        for (int a = 0; a < N2; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = ocr[a] * im_U[a] + oci[a] * re_U[a];
        }
    }
    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}

/* Propagate doubles using diagonal vs. coupling sparce algorithm */
void propagate_vec_coupling_S_doubles_ES(t_non* non, float* Hamiltonian_i, float* cr, float* ci, int m) {
    int N = non->singles;
    int N2 = N * (N + 1) / 2;
    int N2old= N * ( N + 1 ) / 2;
    const float f = non->deltat * icm2ifs * twoPi / m;
    float* H0 = calloc(N2, sizeof(float));
    float* H1 = calloc(N * N / 2, sizeof(float));
    int* col = calloc(N * N / 2, sizeof(int));
    int* row = calloc(N * N / 2, sizeof(int));
    float* re_U = calloc(N2, sizeof(float));
    float* im_U = calloc(N2, sizeof(float));
    float* ocr = calloc(N2, sizeof(float));
    float* oci = calloc(N2, sizeof(float));
    float co,si;

    /* Build Hamiltonians H0 (diagonal) and H1 (coupling) */
    for (int a = 0; a < N; a++) {
        const int indexa = Sindex(a, a, N);
        for (int b = a; b < N; b++) {
            int index = Sindex(a, b, N);
            H0[index] = Hamiltonian_i[indexa] + Hamiltonian_i[Sindex(b, b, N)]; // Diagonal
        }
    }

    /* Build Hamiltonian H1 (coupling) */
    int kmax = 0;
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            int index = b + a * ((N << 1) - a - 1) / 2; // Part of Sindex, but b > a is always true here

            if (fabsf(Hamiltonian_i[index]) > non->couplingcut) {
                H1[kmax] = Hamiltonian_i[index];
                col[kmax] = a, row[kmax] = b;
                kmax++;
            }
        }
    }

    /* Exponentiate diagonal [U=exp(-i/2h H0 dt)] */
    for (int a = 0; a < N2; a++) {
        re_U[a] = cosf(0.5f * H0[a] * f);
        im_U[a] = -sinf(0.5f * H0[a] * f);
    }

    // Loop over Trotter steps
    for (int i = 0; i < m; i++) {

        // Multiply on vector first time
        for (int a = 0; a < N2; a++) {
            ocr[a] = cr[a] * re_U[a] - ci[a] * im_U[a];
            oci[a] = cr[a] * im_U[a] + ci[a] * re_U[a];
        }
        
        // Account for couplings 
        // Loop over couplings 
        for (int k = 0; k < kmax; k++) {
            int a = col[k];
            int b = row[k];
            float J = H1[k] * f;

            // Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> 
            // TODO speedup
            int c=0;
            int aNsum = a * ((N << 1) - a - 1) / 2;  // copied part of Sindex
            int bNsum = b * ((N << 1) - b - 1) / 2;  // copied part of Sindex

            /* Loop over wave functions <ca|Hab|cb> and <cb|Hba|ca> */

            // c < a,b
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < a; c++) {
                int index1 = a + c * ((N << 1) - c - 1) / 2; // part of Sindex, but always a > c
                int index2 = b + c * ((N << 1) - c - 1) / 2; // part of Sindex, but always b > c
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }
     
                                                                                            // c == a
            for (; c == a; c++) {  // yeah, one iteration. but loop for consistency across all 5
            }
            
            // a < c < b
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < b; c++) { // could be 0 iterations if (b == a + 1)
                int index1 = c + aNsum;
                int index2 = b + c * ((N << 1) - c - 1) / 2;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }

            // c == b
            for (; c == b; c++) {
            }

            // c > a,b
            si = -sinf(J);
            co = sqrtf(1 - si * si);
            for (; c < N; c++) {
                int index1 = c + aNsum;
                int index2 = c + bNsum;
                float cr1 = co * ocr[index1] - si * oci[index2];
                float ci1 = co * oci[index1] + si * ocr[index2];
                float cr2 = co * ocr[index2] - si * oci[index1];
                float ci2 = co * oci[index2] + si * ocr[index1];
                ocr[index1] = cr1, oci[index1] = ci1, ocr[index2] = cr2, oci[index2] = ci2;
            }
        }
	
        // Multiply on vector second time 
        for (int a = 0; a < N2; a++) {
            cr[a] = ocr[a] * re_U[a] - oci[a] * im_U[a];
            ci[a] = ocr[a] * im_U[a] + oci[a] * re_U[a];
        }
    }
    
    free(ocr), free(oci), free(re_U), free(im_U), free(H1), free(H0);
    free(col), free(row);
}



// Diagonalize with LAPACK (destructive version)
void diagonalizeLPD(float* H, float* v, int N) {
    int INFO, lwork;
    float *work, *Hcopy;
    int i, j;
    // Find lwork;
    lwork = -1;
    work = (float *)calloc(1, sizeof(float));
    ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
    lwork = work[0];
    //  printf("LAPACK work dimension %d\n",lwork);
    //  lwork=8*N;
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    Hcopy = (float *)calloc(N * N, sizeof(float));
    // Copy Hamiltonian
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Hcopy[i * N + j] = H[i * N + j];
        }
    }

    // Call LAPACK routine
    ssyev_("V", "U", &N, Hcopy, &N, v, work, &lwork, &INFO);
    //  printf("LAPACK opt. %f %f\n",work[0],work[0]/N);
    if (INFO != 0) {
        printf("Something went wrong trying to diagonalize a matrix...\n");
        exit(0);
    }

    // Move eigenvectors
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            H[i * N + j] = Hcopy[j * N + i]; // Converting from FORTRAN format
        }
    }
    // Free space
    free(Hcopy), free(work);
    return;
}

// Build and diagonalize Hamiltonian
void build_diag_H(float* Hamiltonian_i, float* H, float* e, int N) {
    int a, b, c;
    // Build square Hamiltonian from triagonal matrix
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }

    diagonalizeLPD(H, e, N);
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

// This subroutine nullify all elements of vector for non selected sites
void projection(float* phi, t_non* non) {
    int i;
    for (i = 0; i < non->singles; i++) {
        phi[i] = phi[i] * non->psites[i];
    }
    return;
}

// Test at the start if the Hamiltonian and dipole files are sensible
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

/* Create truncated time-evolution operator */
int time_evolution_mat(t_non* non, float* Hamiltonian_i, float* Ur, float* Ui, int* R, int* C, int m) {
    float f, g;
    float *H, *re_U, *im_U, *e;
    float *cr, *ci;
    float *cnr, *cni;
    int indexA, indexB, N;
    int a, b, c, d;
    int elements;
    N = non->singles;
    f = non->deltat * icm2ifs * twoPi / m;
    H = (float *)calloc(N * N, sizeof(float));
    cr = (float *)calloc(N * N, sizeof(float));
    ci = (float *)calloc(N * N, sizeof(float));
    re_U = (float *)calloc(N, sizeof(float));
    im_U = (float *)calloc(N, sizeof(float));
    e = (float *)calloc(N, sizeof(float));
    cnr = (float *)calloc(N * N, sizeof(float));
    cni = (float *)calloc(N * N, sizeof(float));
    /* Build Hamiltonian */
    for (a = 0; a < N; a++) {
        H[a + N * a] = Hamiltonian_i[a + N * a - (a * (a + 1)) / 2]; // Diagonal
        for (b = a + 1; b < N; b++) {
            H[a + N * b] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
            H[b + N * a] = Hamiltonian_i[b + N * a - (a * (a + 1)) / 2];
        }
    }
    diagonalizeLPD(H, e, N);
    /* Exponentiate [U=exp(-i/h H dt)] */
    for (a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }
    /* Transform to site basis */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] += H[b + a * N] * re_U[b], cni[b + a * N] += H[b + a * N] * im_U[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                cr[a + c * N] += H[b + a * N] * cnr[b + c * N], ci[a + c * N] += H[b + a * N] * cni[b + c * N];
            }
        }
    }
    /* The one exciton propagator has been calculated */
    elements = 0;
    /* Make sparce */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            if ((cr[a + b * N] * cr[a + b * N] + ci[a + b * N] * ci[a + b * N]) > non->thres) {
                Ur[elements] = cr[a + b * N], Ui[elements] = ci[a + b * N];
                R[elements] = a, C[elements] = b;
                elements++;
            }
        }
    }
    free(H), free(cr), free(ci), free(re_U), free(im_U), free(e), free(cnr), free(cni);
    return elements;
}

void propagate_double_sparce(t_non* non, float* Ur, float* Ui, int* R, int* C, float* fr, float* fi, int elements,
                             int m, float* Anh) {
    int a, b, indexA, indexB;
    int N, Nf, s, t, u, v;
    float g, sqrt12;
    float *vr, *vi;
    float Uar, Ubr, Uai, Ubi;
    float *si, *co;
    float f;
    float norm1, norm2;
    float fm;
    int i;

    sqrt12 = 1.0 / sqrt2;
    f = non->deltat * icm2ifs * twoPi;
    fm = f * 0.5 / m;
    N = non->singles;
    Nf = non->singles * (non->singles + 1) / 2;
    vr = (float *)calloc(Nf, sizeof(float));
    vi = (float *)calloc(Nf, sizeof(float));
    co = (float *)calloc(N, sizeof(float));
    si = (float *)calloc(N, sizeof(float));

    if (non->anharmonicity != 0) {
        for (i = 0; i < N; i++) {
            co[i] = cos(fm * non->anharmonicity), si[i] = sin(fm * non->anharmonicity);
        }
    }
    else {
        for (i = 0; i < N; i++) {
            co[i] = cos(fm * Anh[i]), si[i] = sin(fm * Anh[i]);
        }
    }
    /* Repeat m times */
    for (i = 0; i < m; i++) {

        /* Anharmonicity */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = fr[indexA];
                vi[indexA] = fi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = co[a] * vr[indexA] - si[a] * vi[indexA];
            fi[indexA] = co[a] * vi[indexA] + si[a] * vr[indexA];
        }
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = 0;
                vi[indexA] = 0;
            }
        }

        for (a = 0; a < elements; a++) {
            s = R[a], t = C[a], Uar = Ur[a], Uai = Ui[a];
            for (b = 0; b <= a; b++) {
                u = R[b], v = C[b], Ubr = Ur[b], Ubi = Ui[b];
                g = 1;
                indexA = Sindex(s, u, N), indexB = Sindex(t, v, N);
                if (s != u) g = sqrt2;
                if (v != t) g = sqrt2;
                if (s != u && v != t) g = 1;
                vr[indexA] += (Uar * Ubr - Uai * Ubi) * fr[indexB] * g;
                vr[indexA] -= (Uai * Ubr + Uar * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uar * Ubr - Uai * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uai * Ubr + Uar * Ubi) * fr[indexB] * g;
            }
        }
        /* Anharmonicity */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                fr[indexA] = vr[indexA];
                fi[indexA] = vi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = co[a] * vr[indexA] - si[a] * vi[indexA];
            fi[indexA] = co[a] * vi[indexA] + si[a] * vr[indexA];
        }
    }
    free(si);
    free(co);
    free(vr);
    free(vi);
}

void propagate_double_sparce_ES(t_non* non, float* Ur, float* Ui, int* R, int* C, float* fr, float* fi, int elements,int m) {
    int a, b, indexA, indexB;
    int N, Nf, s, t, u, v;
    float g, sqrt12;
    float *vr, *vi;
    float Uar, Ubr, Uai, Ubi;
    float f;
    float norm1, norm2;
    float fm;
    int i;

    sqrt12 = 1.0 / sqrt2;
    f = non->deltat * icm2ifs * twoPi;
    fm = f * 0.5 / m;
    N = non->singles;
    Nf = non->singles * (non->singles + 1) / 2;
    vr = (float *)calloc(Nf, sizeof(float));
    vi = (float *)calloc(Nf, sizeof(float));

        /* Repeat m times */
    for (i = 0; i < m; i++) {

        /* Copy */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = fr[indexA];
                vi[indexA] = fi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = vr[indexA] ;
            fi[indexA] = vi[indexA] ;
        }
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                vr[indexA] = 0;
                vi[indexA] = 0;
            }
        }

        for (a = 0; a < elements; a++) {
            s = R[a], t = C[a], Uar = Ur[a], Uai = Ui[a];
            for (b = 0; b <= a; b++) {
                u = R[b], v = C[b], Ubr = Ur[b], Ubi = Ui[b];
                g = 1;
                indexA = Sindex(s, u, N), indexB = Sindex(t, v, N);
                if (s != u) g = 0; // No coupling to double excited states
                if (v != t) g = 0;
                if (s != u && v != t) g = 1;
                vr[indexA] += (Uar * Ubr - Uai * Ubi) * fr[indexB] * g;
                vr[indexA] -= (Uai * Ubr + Uar * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uar * Ubr - Uai * Ubi) * fi[indexB] * g;
                vi[indexA] += (Uai * Ubr + Uar * Ubi) * fr[indexB] * g;
            }
	            }
        /* Copy */
        for (a = 0; a < N; a++) {
            for (b = 0; b <= a; b++) {
                indexA = Sindex(a, b, N);
                fr[indexA] = vr[indexA];
                fi[indexA] = vi[indexA];
            }
        }
        for (a = 0; a < N; a++) {
            indexA = Sindex(a, a, N);
            fr[indexA] = vr[indexA] ;
            fi[indexA] = vi[indexA] ;
        }
    }
    free(vr);
    free(vi);
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
