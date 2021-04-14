#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "pop_single_t2.h"
#include <stdarg.h>

// Print population results to file
void pop_print(char* filename, float* pop, t_non* non, int sampleCount) {
    FILE* out = fopen(filename, "w");
    for (int t2 = 0; t2 < non->tmax2; t2++) {
            pop[t2] /= sampleCount;
            // .csv format
            fprintf(out, "%f,%e\n", t2 * non->deltat, pop[t2]);
    }
    fclose(out);
}

void propagate_NISE(t_non *non, float *H, float *e, float *re_U, float *im_U, float *cr, float *ci) {
    float f;
    int index, N;
    float re, im;

    N = non->singles;
    f = non->deltat * icm2ifs * twoPi;

    // Exponentiate [U=exp(-i/h H dt)]
    for (int a = 0; a < N; a++) {
        re_U[a] = cos(e[a] * f);
        im_U[a] = -sin(e[a] * f);
    }

    // Transfer to eigen basis
    matrix_on_vector(H,cr,ci,N);
    // Multiply with matrix exponent
    vector_on_vector(re_U,im_U,cr,ci,N);
    // Transfer back to site basis
    trans_matrix_on_vector(H,cr,ci,N);
}

// void propagate_prezhdo() {

// }

// void propagate_alt() {

// }

void pop_single_t2(t_non* non) {
    // Initialize each process base variables
    int N = non->singles;
    int nn2 = N * (N + 1) / 2;
    int N2 = N * N;
    int ti, tm, samples;
    int sampleCount;
    float *Hamil_i_e;
    float *H_new;
    float *H_old;
    float *re_U;
    float *im_U;
    float *e;
    float *cr_nise;
    float *ci_nise;
    float *cr_prezhdo;
    float *ci_prezhdo;
    float *cr_alt;
    float *ci_alt;
    float *pop_nise;
    float *pop_prezhdo;
    float *pop_alt;

    Hamil_i_e = (float *) calloc(nn2, sizeof(float));
    H_new = (float *) calloc(N2, sizeof(float));
    H_old = (float *) calloc(N2, sizeof(float));
    re_U = (float *) calloc(N, sizeof(float));
    im_U = (float *) calloc(N, sizeof(float));
    e = (float *) calloc(N, sizeof(float));
    cr_nise = (float *) calloc(N, sizeof(float));
    ci_nise = (float *) calloc(N, sizeof(float));
    cr_prezhdo = (float *) calloc(N, sizeof(float));
    ci_prezhdo = (float *) calloc(N, sizeof(float));
    cr_alt = (float *) calloc(N, sizeof(float));
    ci_alt = (float *) calloc(N, sizeof(float));
    pop_nise = (float *) calloc(non->tmax2, sizeof(float));
    pop_prezhdo = (float *) calloc(non->tmax2, sizeof(float));
    pop_alt = (float *) calloc(non->tmax2, sizeof(float));

    sampleCount = non->end - non->begin;
    
    // Set the initial wavefunctions
    clearvec(cr_nise, N);
    clearvec(ci_nise, N);
    clearvec(cr_prezhdo, N);
    clearvec(ci_prezhdo, N);
    clearvec(cr_alt, N);
    clearvec(ci_alt, N);
    cr_nise[0] = cr_prezhdo[0] = cr_alt[0] = 1.0;
    ci_nise[0] = ci_prezhdo[0] = ci_alt[0] = 0.0;
    pop_nise[0] = pop_prezhdo[0] = pop_alt[0] = 1.0;

    // Read the Hamiltonian file
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    // Loop over samples
    for (samples = 0; samples < sampleCount; samples++) {
        ti = samples * non->sample;
        int tj = ti + non->tmax1;

        // Load first Hamiltonian
        if (read_He(non, Hamil_i_e, H_traj, tj) != 1) {
            printf("Hamiltonian trajectory file too short, could not fill buffer!\n");
            exit(1);
        }
        build_diag_H(Hamil_i_e, H_new, e, N);

        // Start integrating the Schrödinger equation
        // NISE:
        for (int t2 = 0; t2 < non->tmax2 - 1; t2++) {
            int tm = tj + t2;
            // Copy old Hamiltonian
            copyvec(H_new, H_old, N2);
            // Load new Hamiltonian
            if (read_He(non, Hamil_i_e, H_traj, tm) != 1) {
                printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                exit(1);
            }
            build_diag_H(Hamil_i_e, H_new, e, N);

            // Run NISE
            propagate_NISE(non, H_new, e, re_U, im_U, cr_nise, ci_nise);
            float pop = cr_nise[0]*cr_nise[0] + ci_nise[0]*ci_nise[0];
            // Sum over bath (unnormalised)
            pop_nise[t2 + 1] += pop;

            // Run Prezhdo

            // Run alt

        }
    }
    char* fn = "pop_t2.csv";
    pop_print(fn, pop_nise, non, sampleCount);

    free(Hamil_i_e), free(H_new), free(H_old);
    free(re_U), free(im_U), free(e);
    free(cr_nise), free(ci_nise);
    free(cr_prezhdo), free(ci_prezhdo);
    free(cr_alt), free(ci_alt);
    free(pop_nise), free(pop_prezhdo), free(pop_alt);
}