#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "read_trajectory.h"
#include "randomlib.h"
#include "util/asprintf.h"

/* Open Files */
void open_files(t_non *non,FILE **H_traj,FILE **mu_traj,FILE **Cfile){
  /* Open Trajectory files */
  *H_traj=fopen(non->energyFName,"rb");
  if (*H_traj==NULL){
    printf("Hamiltonian file not found!\n");
    exit(1);
  }

  *mu_traj=fopen(non->dipoleFName,"rb");
  if (*mu_traj==NULL){
    printf("Dipole file %s not found!\n",non->dipoleFName);
    exit(1);
  }

  /* Open file with cluster information if appicable */
  if (non->cluster!=-1){
    *Cfile=fopen("Cluster.bin","rb");
    if (*Cfile==NULL){
      printf("Cluster option was activated but no Cluster.bin file provided.\n");
      printf("Please, provide cluster file or remove Cluster keyword from\n");
      printf("input file.\n");
      exit(0);
    }
  }
  return;
}

/* Read coupling, this is done if the coupling and transition-dipoles are */
/* time-independent and only one snapshot is stored */
void read_coupling(t_non *non,FILE *C_traj,FILE *mu_traj,float* Hamil_i_e,float *mu_xyz){ 
  int x;
  if (!strcmp(non->hamiltonian,"Coupling")){
    C_traj=fopen(non->couplingFName,"rb");
    if (C_traj==NULL){
      printf("Coupling file not found!\n");
      exit(1);
    }
    if (read_He(non,Hamil_i_e,C_traj,-1)!=1){
      printf("Coupling trajectory file to short, could not fill buffer!!!\n");
      exit(1);
    }
    fclose(C_traj);
    /* Reading in single fixed transition dipole vector matrix */
    for (x=0;x<3;x++){
      if (read_mue(non,mu_xyz+non->singles*x,mu_traj,0,x)!=1){
         printf("Dipole trajectory file to short, could not fill buffer!!!\n");
         printf("ITIME %d %d\n",0,x);
         exit(1);
      }
    }
  }
  return;
}

/* Read Hamiltonian */
void read_Hamiltonian(t_non *non,float *Hamil_i_e,FILE *H_traj,int pos){
    if (!strcmp(non->hamiltonian,"Coupling")){
        if (read_Dia(non,Hamil_i_e,H_traj,pos)!=1){
            printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
            exit(1);
        }
    } else {
        if (read_He(non,Hamil_i_e,H_traj,pos)!=1){
           printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
           exit(1);
        }
    }
    return;
}

/* Read the Dipole vector */
void read_dipole(t_non *non,FILE *mu_traj,float *mu_eg,float *mu_xyz,int x,int pos){
    if (!strcmp(non->hamiltonian,"Coupling")){
        copyvec(mu_xyz+non->singles*x,mu_eg,non->singles);
    } else {
        if (read_mue(non,mu_eg,mu_traj,pos,x)!=1){
           printf("Dipole trajectory file to short, could not fill buffer!!!\n");
           printf("JTIME %d %d\n",pos,x);
           exit(1);
        }
    }
    return;
}

/* Core functions */

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
    /* Remove couplings below coupling cut */
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
    /* Find position */
    fseek(FH, pos * (sizeof(int) + sizeof(float) * (3 * non->singles + 3 * non->singles * non->doubles)) + sizeof(float)
          * x * non->singles,SEEK_SET);
          /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    /* Read single excitation Dipoles */
    fread(mue, sizeof(float), non->singles, FH);

    if (control!=1){
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n",pos,x);
        exit(1);
    }

    return control;
}

/* Read Raman tensor */
int read_alpha(t_non* non, float* alpha, FILE* FH, int pos, int x) {
    int control;
    int t;
    control = 0;

    if (non->doubles == 0){
        /* Find position in file to read from */
        fseek(FH, pos * (sizeof(int) + sizeof(float) * 6 * non->singles) + sizeof(float)
          * x * non->singles,SEEK_SET); //6 due to six elements for Raman tensor
          /* Read time */
        if (fread(&t, sizeof(int), 1, FH)) control = 1;
        /* Read single excitation raman tensor values */
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
    /* Find position */
    fseek(FH, pos * (sizeof(int) + sizeof(float) * (3 * non->singles)) + sizeof(float) * x * non->singles,SEEK_SET);
    /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    /* Read single excitation Dipoles */
    fread(over, sizeof(float), non->singles, FH);
    return control;
}

/* Read transition dipole */
void muread(t_non* non, float* leftnr, int ti, int x, FILE* mu_traj) {
    /* Read mu(ti) */
    if (read_mue(non, leftnr, mu_traj, ti, x) != 1) {
        printf("Dipole trajectory file to short, could not fill buffer!!!\n");
        printf("ITIME %d %d\n", ti, x);
        exit(1);
    }
    return;
}

/* Read transition dipole (including option for Coupling storrage) */
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
    /* Find position */
    fseek(FH, pos * (sizeof(int) + sizeof(int)),SEEK_SET);
    /* Read time */
    if (fread(&t, sizeof(int), 1, FH)) control = 1;
    /* Read single excitation Dipoles */
    fread(cl, sizeof(int), 1, FH);
    return control;
}
