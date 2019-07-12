#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "calc_2DES.h"
#include <stdarg.h>

/**
 * Method that logs a message, in which the message can be formatted like printf accepts.
 */
void log_item(char *msgFormat, ...) {
    // Parse parameters
    va_list args;
    va_start(args, msgFormat);

    // Write to log
    FILE* log = fopen("NISE.log", "a");
    if(log == NULL) {
        printf("Could not open log file!");
        exit(1);
    }
    vfprintf(log, msgFormat, args);
    fclose(log);

    va_end(args);
}

void calc_2DES(t_non* non) {
    /* Define arrays! */

    // Aid arrays 
    float *leftrr, *leftri, *leftnr, *leftni;
    float *leftrr_o, *leftri_o, *leftnr_o, *leftni_o;
    float *rightrr, *rightri, *rightnr, *rightni;
    float *rightrr_o, *rightri_o, *rightnr_o, *rightni_o;
    float *t1rr, *t1ri, *t1nr, *t1ni;
    float *Hamil_i_e;
    float *Anh, *over;
    float *mut2, *mut3r, *mut3i, *mut4;
    float *mut3r_o, *mut3i_o;
    float *fr, *fi, *fr_o, *fi_o;
    float *ft1r, *ft1i, *ft1r_o, *ft1i_o;
    float t3nr, t3ni, t3rr, t3ri;
    float rrI, riI, rrII, riII;
    float *Urs, *Uis;
    int *Rs, *Cs;
    float* mu_xyz;
    float* pol = 0; /* Currently dummy vector that can be used to change coordinate system in the future */

    /* Integers */
    int cl, Ncl;

    /* File handles */
    FILE *H_traj, *mu_traj;
    FILE* Cfile;
    FILE* C_traj;
    FILE *A_traj, *mu2_traj;
    FILE* outttwo;

    /* Time parameters */
    time_t timeNow, timeStart;
    /* Initialize time */
    time(&timeNow);
    time(&timeStart);

    float shift1 = (non->max1 + non->min1) / 2;
    printf("Frequency shift %f.\n", shift1);
    non->shifte = shift1;
    non->shiftf = 2 * shift1;

    int maxtt = non->tmax1 * non->tmax3;

    // 2D response function parallel
    float* rrIpar = calloc(maxtt, sizeof(float));
    float* riIpar = calloc(maxtt, sizeof(float));
    float* rrIIpar = calloc(maxtt, sizeof(float));
    float* riIIpar = calloc(maxtt, sizeof(float));
    // 2D response function perpendic.
    float* rrIper = calloc(maxtt, sizeof(float));
    float* riIper = calloc(maxtt, sizeof(float));
    float* rrIIper = calloc(maxtt, sizeof(float));
    float* riIIper = calloc(maxtt, sizeof(float));
    // 2D response function cross
    float* rrIcro = calloc(maxtt, sizeof(float));
    float* riIcro = calloc(maxtt, sizeof(float));
    float* rrIIcro = calloc(maxtt, sizeof(float));
    float* riIIcro = calloc(maxtt, sizeof(float));

    float* lt_gb_se = calloc(non->tmax1 * non->tmax3, sizeof(float));
    float* lt_ea = calloc(non->tmax1 * non->tmax3, sizeof(float));

    for (int t1 = 0; t1 < non->tmax1; t1++) {
        for (int t3 = 0; t3 < non->tmax3; t3++) {
            lt_gb_se[t1 + t3 * non->tmax1] = (float) exp(-(double)(t1 + t3) * non->deltat / (2 * non->lifetime));
            lt_ea[t1 + t3 * non->tmax1] = (float) exp(-(double)(t1 + t3) * non->deltat / (2 * non->lifetime));
        }
    }

    const int nn2 = non->singles * (non->singles + 1) / 2;

    Hamil_i_e = (float *)calloc(nn2, sizeof(float));

    // Reserve memory for 2D calculation
    Anh = (float *)calloc(non->singles, sizeof(float));
    over = (float *)calloc(non->singles, sizeof(float));
    leftrr = (float *)calloc(non->singles, sizeof(float));
    leftri = (float *)calloc(non->singles, sizeof(float));
    leftnr = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    leftni = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    leftrr_o = (float *)calloc(non->singles, sizeof(float));
    leftri_o = (float *)calloc(non->singles, sizeof(float));
    leftnr_o = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    leftni_o = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    rightrr = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    rightri = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    rightnr = (float *)calloc(non->singles, sizeof(float));
    rightni = (float *)calloc(non->singles, sizeof(float));
    rightrr_o = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    rightri_o = (float *)calloc(non->singles * non->tmax1, sizeof(float));
    rightnr_o = (float *)calloc(non->singles, sizeof(float));
    rightni_o = (float *)calloc(non->singles, sizeof(float));
    mut2 = (float *)calloc(non->singles, sizeof(float));
    mut3r = (float *)calloc(non->singles, sizeof(float));
    mut3i = (float *)calloc(non->singles, sizeof(float));
    mut3r_o = (float *)calloc(non->singles, sizeof(float));
    mut3i_o = (float *)calloc(non->singles, sizeof(float));
    mut4 = (float *)calloc(non->singles, sizeof(float));
    t1rr = (float *)calloc(non->tmax1, sizeof(float));
    t1ri = (float *)calloc(non->tmax1, sizeof(float));
    t1nr = (float *)calloc(non->tmax1, sizeof(float));
    t1ni = (float *)calloc(non->tmax1, sizeof(float));
    fr = (float *)calloc(nn2, sizeof(float));
    fi = (float *)calloc(nn2, sizeof(float));
    fr_o = (float *)calloc(nn2, sizeof(float));
    fi_o = (float *)calloc(nn2, sizeof(float));
    ft1r = (float *)calloc(nn2 * non->tmax1, sizeof(float));
    ft1i = (float *)calloc(nn2 * non->tmax1, sizeof(float));
    ft1r_o = (float *)calloc(nn2 * non->tmax1, sizeof(float));
    ft1i_o = (float *)calloc(nn2 * non->tmax1, sizeof(float));
    Urs = (float *)calloc(non->singles * non->singles, sizeof(float));
    Uis = (float *)calloc(non->singles * non->singles, sizeof(float));
    Rs = (int *)calloc(non->singles * non->singles, sizeof(int));
    Cs = (int *)calloc(non->singles * non->singles, sizeof(int));
    mu_xyz = (float *)calloc(non->singles * 3, sizeof(float));

    /* Open Trajectory files */
    H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        printf("Hamiltonian file not found!\n");
        exit(1);
    }

    mu_traj = fopen(non->dipoleFName, "rb");
    if (mu_traj == NULL) {
        printf("Dipole file %s not found!\n", non->dipoleFName);
        exit(1);
    }

    /* Open file with cluster information if applicable */
    if (non->cluster != -1) {
        Cfile = fopen("Cluster.bin", "rb");
        if (Cfile == NULL) {
            printf("Cluster option was activated but no Cluster.bin file provided.\n");
            printf("Please, provide cluster file or remove Cluster keyword from\n");
            printf("input file.\n");
            exit(0);
        }
        Ncl = 0; // Counter for snapshots calculated
    }

    /* Open file for fluctuating anharmonicities and sequence transition dipoles if needed */
    if (non->anharmonicity == 0 && (!strcmp(non->technique, "2DUVvis") || (!strcmp(non->technique, "EAUVvis")) || (!
        strcmp(non->technique, "noEAUVvis")) || (!strcmp(non->technique, "GBUVvis")) || (!strcmp(
        non->technique, "SEUVvis")))) {
        A_traj = fopen(non->anharFName, "rb");
        if (A_traj == NULL) {
            printf("Anharmonicity file %s not found!\n", non->anharFName);
            exit(1);
        }

        mu2_traj = fopen(non->overdipFName, "rb");
        if (mu2_traj == NULL) {
            printf("Overtone dipole file %s not found!\n", non->overdipFName);
            exit(1);
        }
    }

    /* Read coupling */
    if (!strcmp(non->hamiltonian, "Coupling")) {
        C_traj = fopen(non->couplingFName, "rb");
        if (C_traj == NULL) {
            printf("Coupling file not found!\n");
            exit(1);
        }
        if (read_He(non, Hamil_i_e, C_traj, -1) != 1) {
            printf("Coupling trajectory file to short, could not fill buffer!!!\n");
            exit(1);
        }
        fclose(C_traj);
        for (int x = 0; x < 3; x++) {
            if (read_mue(non, mu_xyz + non->singles * x, mu_traj, 0, x) != 1) {
                printf("Dipole trajectory file to short, could not fill buffer!!!\n");
                printf("ITIME %d %d\n", 0, x);
                exit(1);
            }
        }
    }

    /* Open new log file */
    FILE* logFile = fopen("NISE.log", "w");
    if (logFile == NULL) {
        printf("Could not open log file! Disk full?\n");
        exit(1);
    }
    fprintf(logFile, "Log\n");
    fclose(logFile);

    /* Do calculation */
    int totalSampleCount = (non->length - non->tmax1 - non->tmax2 - non->tmax3 - 1) / non->sample + 1;
    if (totalSampleCount > 0) {
        printf("Making %d samples!\n", totalSampleCount);
    }
    else {
        printf("Insufficient data to calculate spectrum.\n");
        printf("Please, lower max times or provide longer\n");
        printf("trajectory.\n");
        exit(1);
    }

    if (non->end == 0) non->end = totalSampleCount;
    int sampleCount = non->end - non->begin;
    if (sampleCount == 0) {
        // Avoid dividing by zero
        sampleCount = 1;
    }

    log_item("Begin sample: %d, End sample: %d.\n", non->begin, non->end);

    /* Loop over samples */
//pragma?
    for (int currentSample = non->begin; currentSample < non->end; currentSample++) {

        /* Log time */
        time_t timeNew;
        time(&timeNew);
        char* timeText = time_diff(timeNow, timeNew);
        log_item("Starting sample %d\n%s", currentSample, timeText);
        free(timeText);
        timeNow = timeNew;

        /* Calculate 2DIR response */
        int tj = currentSample * non->sample + non->tmax1;
        int tk = tj + non->tmax2;

        // Cluster checking if we are using it
        if (non->cluster != -1) {
            if (read_cluster(non, tj, &cl, Cfile) != 1) {
                printf("Cluster trajectory file to short, could not fill buffer!!!\n");
                printf("ITIME %d\n", tj);
                exit(1);
            }

            // If we do not have the current cluster, skip calculation
            if (non->cluster != cl) continue;
            
            Ncl++; // increase counter
        }

        /* Loop over Polarization components */
//        #pragma omp parallel for default(NONE) \
//            private(mut2, t1nr, t1ni, px) \
//            shared(tj, non)

        for (int molPol = 0; molPol < 21; molPol++) {
            int px[4];
            polar(px, molPol);
            mureadE(non, mut2, tj, px[1], mu_traj, mu_xyz, pol);

            /* Ground state bleach (GB) kI and kII */
            for (int t1 = 0; t1 < non->tmax1; t1++) {
                /* Read dipoles at time 0 */
                mureadE(non, leftnr + t1 * non->singles, tj - t1, px[0], mu_traj, mu_xyz, pol);
                clearvec(leftni + t1 * non->singles, non->singles);
                clearvec(leftnr_o + t1 * non->singles, non->singles);
                clearvec(leftni_o + t1 * non->singles, non->singles);

                /* Propagate */
                for (int tm = 0; tm < t1; tm++) {
                    /* Read Hamiltonian */
                    if (read_He(non, Hamil_i_e, H_traj, tj - t1 + tm) != 1) {
                        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                        exit(1);
                    }
                    if (non->propagation == 1) propagate_vec_coupling_S(
                        non, Hamil_i_e, leftnr + t1 * non->singles, leftni + t1 * non->singles, non->ts, 1);
                    if (non->propagation == 0) propagate_vec_DIA_S(non, Hamil_i_e, leftnr + t1 * non->singles,
                                                                   leftni + t1 * non->singles, 1);

                }
            }

//#pragma omp parallel for
            for (int t1 = 0; t1 < non->tmax1; t1++) {
                t1nr[t1] = 0, t1ni[t1] = 0;
                for (int i = 0; i < non->singles; i++) {
                    t1nr[t1] += mut2[i] * leftnr[i + t1 * non->singles];
                    t1ni[t1] += mut2[i] * leftni[i + t1 * non->singles];
                }
            }

            /* Combine with evolution during t3 */
            mureadE(non, mut3r, tk, px[2], mu_traj, mu_xyz, pol);
            clearvec(mut3i, non->singles);
            for (int t3 = 0; t3 < non->tmax3; t3++) {
                int tl = tk + t3;
                mureadE(non, mut4, tl, px[3], mu_traj, mu_xyz, pol);
                t3nr = 0, t3ni = 0;
                for (int i = 0; i < non->singles; i++) {
                    t3nr += mut4[i] * mut3r[i];
                    t3ni += mut4[i] * mut3i[i];
                }
                /* Calculate GB contributions */
                if ((!strcmp(non->technique, "GBUVvis")) || (!strcmp(non->technique, "2DUVvis")) || (!strcmp(
                    non->technique, "noEAUVvis"))) {
//#pragma omp parallel for private(tt,polWeight)
                    for (int t1 = 0; t1 < non->tmax1; t1++) {
                        int tt = non->tmax1 * t3 + t1;
                        float polWeight = polarweight(0, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIpar[tt] -= (t3ni * t1nr[t1] - t3nr * t1ni[t1]) * polWeight;
                        riIpar[tt] -= (t3nr * t1nr[t1] + t3ni * t1ni[t1]) * polWeight;
                        rrIIpar[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIpar[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                        polWeight = polarweight(1, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIper[tt] -= (t3ni * t1nr[t1] - t3nr * t1ni[t1]) * polWeight;
                        riIper[tt] -= (t3nr * t1nr[t1] + t3ni * t1ni[t1]) * polWeight;
                        rrIIper[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIper[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                        polWeight = polarweight(2, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIcro[tt] -= (t3ni * t1nr[t1] - t3nr * t1ni[t1]) * polWeight;
                        riIcro[tt] -= (t3nr * t1nr[t1] + t3ni * t1ni[t1]) * polWeight;
                        rrIIcro[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIcro[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                    }
                }

                /* Read Hamiltonian */
                if (read_He(non, Hamil_i_e, H_traj, tl) != 1) {
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }

                /* Propagate */
                if (non->propagation == 1) propagate_vec_coupling_S(non, Hamil_i_e, mut3r, mut3i, non->ts, 1);
                if (non->propagation == 0) propagate_vec_DIA_S(non, Hamil_i_e, mut3r, mut3i, 1);
            }

            /* Stimulated emission (SE) */
            /* Calculate evolution during t2 */
            mureadE(non, leftrr, tj, px[1], mu_traj, mu_xyz, pol);
            clearvec(leftri, non->singles);
            for (int t2 = 0; t2 < non->tmax2; t2++) {
                int tm = tj + t2;
                if (read_He(non, Hamil_i_e, H_traj, tm) != 1) {
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }

//#pragma omp parallel for
                for (int t1 = -1; t1 < non->tmax1; t1++) {
                    if (t1 != -1) {
                        propagate_vec_DIA(non, Hamil_i_e, leftnr + t1 * non->singles, leftni + t1 * non->singles,
                                          1);
                    }
                    else {
                        propagate_vec_DIA(non, Hamil_i_e, leftrr, leftri, 1);
                    }
                }
            }

            /* Read dipole for third interaction */
            mureadE(non, mut3r, tk, px[2], mu_traj, mu_xyz, pol);

            if ((!strcmp(non->technique, "EAUVvis")) || (!strcmp(non->technique, "2DUVvis"))) {
                if (non->anharmonicity == 0) {
                    read_over(non, over, mu2_traj, tk, px[2]);
                }
                /* T2 propagation ended store vectors needed for EA */
//#pragma omp parallel for
                for (int t1 = -1; t1 < non->tmax1; t1++) {
                    if (t1 == -1) {
                        dipole_double(non, mut3r, leftrr, leftri, fr, fi, over);
                    }
                    else {
                        dipole_double(non, mut3r, leftnr + t1 * non->singles, leftni + t1 * non->singles,
                                      ft1r + t1 * nn2, ft1i + t1 * nn2, over);
                    }
                }

                copyvec(leftnr, rightrr, non->tmax1 * non->singles);
                copyvec(leftni, rightri, non->tmax1 * non->singles);
//#pragma omp parallel for
                for (int i = 0; i < non->tmax1 * non->singles; i++) rightri[i] = -rightri[i];
                copyvec(leftrr, rightnr, non->singles);
                copyvec(leftri, rightni, non->singles);
                for (int i = 0; i < non->singles; i++) rightni[i] = -rightni[i];
            }

            clearvec(mut3i, non->singles);

            /* Calculate right side of rephasing diagram */
//#pragma omp parallel for
            for (int t1 = -1; t1 < non->tmax1; t1++) {
                if (t1 != -1) {
                    t1rr[t1] = 0, t1ri[t1] = 0;
                    for (int i = 0; i < non->singles; i++) {
                        t1rr[t1] += leftnr[i + t1 * non->singles] * mut3r[i];
                        t1ri[t1] -= leftni[i + t1 * non->singles] * mut3r[i];
                    }
                }
                else {
                    /* Calculate right side of nonrephasing diagram */
                    t3nr = 0, t3ni = 0;
                    for (int i = 0; i < non->singles; i++) {
                        t3nr += leftrr[i] * mut3r[i];
                        t3ni -= leftri[i] * mut3r[i];
                    }
                }
            }

            /* Combine with evolution during t3 */
            for (int t3 = 0; t3 < non->tmax3; t3++) {
                int tl = tk + t3;
                mureadE(non, mut4, tl, px[3], mu_traj, mu_xyz, pol);

//#pragma omp parallel for
                for (int t1 = -1; t1 < non->tmax1; t1++) {
                    if (t1 == -1) {
                        /* Calculate left side of rephasing diagram */
                        t3rr = 0, t3ri = 0;
                        for (int i = 0; i < non->singles; i++) {
                            t3rr += mut4[i] * leftrr[i];
                            t3ri += mut4[i] * leftri[i];
                        }
                    }
                    else {
                        /* Calculate left side of nonrephasing diagram */
                        t1nr[t1] = 0, t1ni[t1] = 0;
                        for (int i = 0; i < non->singles; i++) {
                            t1nr[t1] += leftnr[i + t1 * non->singles] * mut4[i];
                            t1ni[t1] += leftni[i + t1 * non->singles] * mut4[i];
                        }
                    }
                }

                /* Calculate Response */
                if ((!strcmp(non->technique, "SEUVvis")) || (!strcmp(non->technique, "2DUVvis")) || (!strcmp(
                    non->technique, "noEAUVvis"))) {
//#pragma omp parallel for private(tt,polWeight)
                    for (int t1 = 0; t1 < non->tmax1; t1++) {
                        int tt = non->tmax1 * t3 + t1;
                        float polWeight = polarweight(0, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIpar[tt] -= (t3rr * t1ri[t1] + t3ri * t1rr[t1]) * polWeight;
                        riIpar[tt] -= (t3rr * t1rr[t1] - t3ri * t1ri[t1]) * polWeight;
                        rrIIpar[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIpar[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                        polWeight = polarweight(1, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIper[tt] -= (t3rr * t1ri[t1] + t3ri * t1rr[t1]) * polWeight;
                        riIper[tt] -= (t3rr * t1rr[t1] - t3ri * t1ri[t1]) * polWeight;
                        rrIIper[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIper[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                        polWeight = polarweight(2, molPol) * lt_gb_se[t1 + t3 * non->tmax1];
                        rrIcro[tt] -= (t3rr * t1ri[t1] + t3ri * t1rr[t1]) * polWeight;
                        riIcro[tt] -= (t3rr * t1rr[t1] - t3ri * t1ri[t1]) * polWeight;
                        rrIIcro[tt] -= (t3ni * t1nr[t1] + t3nr * t1ni[t1]) * polWeight;
                        riIIcro[tt] -= (t3nr * t1nr[t1] - t3ni * t1ni[t1]) * polWeight;
                    }
                }


                /* Do Propagation */
                /* Read Hamiltonian */
                if (read_He(non, Hamil_i_e, H_traj, tl) != 1) {
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }
//#pragma omp parallel for shared(non,Hamil_i_e,leftnr,leftni,leftrr,leftri)
                for (int t1 = -1; t1 < non->tmax1; t1++) {
                    if (t1 == -1) {
                        /* Propagate left side rephasing */
                        if (non->propagation == 1) propagate_vec_coupling_S(
                            non, Hamil_i_e, leftrr, leftri, non->ts, 1);
                        if (non->propagation == 0) propagate_vec_DIA_S(non, Hamil_i_e, leftrr, leftri, 1);
                    }
                    else {
                        /* Propagate left side nonrephasing */
                        if (non->propagation == 0) propagate_vec_DIA_S(
                            non, Hamil_i_e, leftnr + t1 * non->singles, leftni + t1 * non->singles, 1);
                        if (non->propagation == 1) propagate_vec_coupling_S(
                            non, Hamil_i_e, leftnr + t1 * non->singles, leftni + t1 * non->singles, non->ts, 1);
                    }
                }
            }

            if ((!strcmp(non->technique, "EAUVvis")) || (!strcmp(non->technique, "2DUVvis"))) {
                /* Excited state absorption (EA) */
                /* Combine with evolution during t3 */
                for (int t3 = 0; t3 < non->tmax3; t3++) {
                    int tl = tk + t3;
                    /* Read Dipole t4 */
                    mureadE(non, mut4, tl, px[3], mu_traj, mu_xyz, pol);
                    if (non->anharmonicity == 0) {
                        read_over(non, over, mu2_traj, tl, px[3]);
                    }
                    /* Multiply with the last dipole */
//#pragma omp parallel for shared(non,mut4,fr,fi,ft1r,ft1i,nn2,leftnr,leftni,over,leftrr,leftri) private(t1)
                    for (int t1 = -1; t1 < non->tmax1; t1++) {
                        if (t1 == -1) {
                            dipole_double_last(non, mut4, fr, fi, leftrr, leftri, over);
                        }
                        else {
                            dipole_double_last(non, mut4, ft1r + t1 * nn2, ft1i + t1 * nn2,
                                               leftnr + t1 * non->singles, leftni + t1 * non->singles, over);
                        }
                    }

                    /* Calculate EA response */
//#pragma omp parallel for shared(leftrr,leftri,non,rightrr,rightri,molPol,lt_ea) private(rrI,riI,rrII,riII,i,tt,polWeight)
                    for (int t1 = 0; t1 < non->tmax1; t1++) {
                        int tt = non->tmax1 * t3 + t1;
                        rrI = 0, riI = 0, rrII = 0, riII = 0;
                        for (int i = 0; i < non->singles; i++) {
                            rrI += leftri[i] * rightrr[i + t1 * non->singles] + leftrr[i] * rightri[i + t1 * non->
                                singles];
                            riI += leftrr[i] * rightrr[i + t1 * non->singles] - rightri[i + t1 * non->singles] *
                                leftri[i];

                            rrII += rightnr[i] * leftni[i + t1 * non->singles] + rightni[i] * leftnr[i + t1 * non->
                                singles];
                            riII += rightnr[i] * leftnr[i + t1 * non->singles] - rightni[i] * leftni[i + t1 * non->
                                singles];
                        }

                        float polWeight = polarweight(0, molPol) * lt_ea[t1 + t3 * non->tmax1];
                        rrIpar[tt] += rrI * polWeight;
                        riIpar[tt] += riI * polWeight;
                        rrIIpar[tt] += rrII * polWeight;
                        riIIpar[tt] += riII * polWeight;
                        polWeight = polarweight(1, molPol) * lt_ea[t1 + t3 * non->tmax1];
                        rrIper[tt] += rrI * polWeight;
                        riIper[tt] += riI * polWeight;
                        rrIIper[tt] += rrII * polWeight;
                        riIIper[tt] += riII * polWeight;
                        polWeight = polarweight(2, molPol) * lt_ea[t1 + t3 * non->tmax1];
                        rrIcro[tt] += rrI * polWeight;
                        riIcro[tt] += riI * polWeight;
                        rrIIcro[tt] += rrII * polWeight;
                        riIIcro[tt] += riII * polWeight;
                    }

                    /* Read Hamiltonian */
                    if (read_He(non, Hamil_i_e, H_traj, tl) != 1) {
                        printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                        exit(1);
                    }

                    /* Propagate vectors left */
                    if (non->anharmonicity == 0) {
                        read_A(non, Anh, A_traj, tl);
                    }

                    /* Propagate */
                    if (non->propagation == 0) {
                        int elements = time_evolution_mat(non, Hamil_i_e, Urs, Uis, Cs, Rs, non->ts);
                        if (currentSample == non->begin) {
                            if (molPol == 0) {
                                if (t3 == 0) {
                                    printf("Sparce matrix efficiency: %f pct.\n",
                                           (1 - (1.0 * elements / (non->singles * non->singles))) * 100);
                                    printf("Pressent tuncation %f.\n",
                                           non->thres / ((non->deltat * icm2ifs * twoPi / non->ts) * (non->deltat *
                                               icm2ifs * twoPi / non->ts)));
                                    printf("Suggested truncation %f.\n", 0.001);
                                }
                            }
                        }

                        // Key parallel loop 1
                        // Initial step, former t1=-1
                        propagate_double_sparce(
                            non, Urs, Uis, Rs, Cs, fr, fi, elements, non->ts, Anh
                        );

                        #pragma omp parallel for \
                            shared(non, Anh, Urs, Uis, Rs, Cs, ft1r, ft1i) \
                            schedule(static, 1)

                        for(int t1 = 0; t1 < non->tmax1; t1++) {
                            propagate_double_sparce(
                                non, Urs, Uis, Rs, Cs, ft1r + t1 * nn2,
                                ft1i + t1 * nn2, elements, non->ts, Anh
                            );
                        }

                        // Propagate vectors right
                        // Key parallel loop 2
                        // Initial step
                        propagate_vec_DIA_S(non, Hamil_i_e, rightnr, rightni, -1);

                        for(int t1 = 0; t1 < non->tmax1; t1++) {
                            propagate_vec_DIA_S(
                                non, Hamil_i_e, rightrr + t1 * non->singles, rightri + t1 * non->singles, -1
                            );
                        }
                    }
                    else if(non->propagation == 1) {
                        // Key parallel loop 1
                        // Initial step
                        propagate_vec_coupling_S_doubles(
                            non, Hamil_i_e, fr, fi, non->ts, Anh
                        );

                        #pragma omp parallel for \
                            shared(non,Hamil_i_e,Anh,ft1r,ft1i) \
                            schedule(static, 1)

                        for (int t1 = 0; t1 < non->tmax1; t1++) {
                            propagate_vec_coupling_S_doubles(
                                non, Hamil_i_e, ft1r + t1 * nn2,
                                ft1i + t1 * nn2, non->ts, Anh
                            );
                        }

                        // Key parallel loop 2
                        // Initial step
                        propagate_vec_coupling_S(
                            non, Hamil_i_e, rightnr, rightni, non->ts, -1
                        );

                        for (int t1 = 0; t1 < non->tmax1; t1++) {
                            propagate_vec_coupling_S(
                                non, Hamil_i_e, rightrr + t1 * non->singles, rightri + t1 * non->singles, non->ts, -1
                            );
                        }
                    }
                }
            }
        }
    }

    /* The calculation is finished, lets write output */
    log_item("Finished Calculating Response!\nWriting to file\n");

    printf("Samples %d\n", sampleCount);

    if (non->cluster != -1) {
        printf("Of %d samples %d belonged to cluster %d.\n", sampleCount, Ncl, non->cluster);
    }

    /* Close Files */
    fclose(mu_traj), fclose(H_traj);
    if ((!strcmp(non->technique, "2DUVvis")) || (!strcmp(non->technique, "GBUVvis")) || (!
        strcmp(non->technique, "SEUVvis")) || (!strcmp(non->technique, "EAUVvis")) || (!strcmp(
        non->technique, "noEAUVvis"))) {
        if (non->anharmonicity == 0) {
            fclose(mu2_traj), fclose(A_traj);
        }
    }
    if (non->cluster != -1) {
        fclose(Cfile);
    }

    /* Print 2D */
    outttwo = fopen("RparI.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIpar[non->tmax1 * t3 + t1] /= sampleCount;
            riIpar[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIpar[non->tmax1 * t3 + t1], riIpar[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    outttwo = fopen("RparII.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIIpar[non->tmax1 * t3 + t1] /= sampleCount;
            riIIpar[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIIpar[non->tmax1 * t3 + t1], riIIpar[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    outttwo = fopen("RperI.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIper[non->tmax1 * t3 + t1] /= sampleCount;
            riIper[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIper[non->tmax1 * t3 + t1], riIper[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    outttwo = fopen("RperII.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIIper[non->tmax1 * t3 + t1] /= sampleCount;
            riIIper[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIIper[non->tmax1 * t3 + t1], riIIper[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    outttwo = fopen("RcroI.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIcro[non->tmax1 * t3 + t1] /= sampleCount;
            riIcro[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIcro[non->tmax1 * t3 + t1], riIcro[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    outttwo = fopen("RcroII.dat", "w");
    for (int t1 = 0; t1 < non->tmax1; t1 += non->dt1) {
        int t2 = non->tmax2;
        for (int t3 = 0; t3 < non->tmax3; t3 += non->dt3) {
            rrIIcro[non->tmax1 * t3 + t1] /= sampleCount;
            riIIcro[non->tmax1 * t3 + t1] /= sampleCount;
            fprintf(outttwo, "%f %f %f %e %e\n", t1 * non->deltat, t2 * non->deltat, t3 * non->deltat,
                    rrIIcro[non->tmax1 * t3 + t1], riIIcro[non->tmax1 * t3 + t1]);
        }
    }
    fclose(outttwo);

    /* Free memory for 2D calculation */
    free(leftrr), free(leftri), free(leftnr), free(leftni);
    free(leftrr_o), free(leftri_o), free(leftnr_o), free(leftni_o);
    free(rightrr), free(rightri), free(rightnr), free(rightni);
    free(rightrr_o), free(rightri_o), free(rightnr_o), free(rightni_o);
    free(t1rr), free(t1ri), free(t1nr), free(t1ni);
    free(mut2), free(mut3r), free(mut3i), free(mut4);
    free(mut3r_o), free(mut3i_o);
    free(fr), free(fi);
    free(ft1r), free(ft1i);
    free(fr_o), free(fi_o), free(ft1r_o), free(ft1i_o);
    free(Urs), free(Uis), free(Rs), free(Cs);
    if (non->anharmonicity == 0) {
        free(Anh), free(over);
    }
    free(rrIpar), free(riIpar);
    free(rrIIpar), free(riIIpar);
    free(rrIper), free(riIper);
    free(rrIIper), free(riIIper);
    free(rrIcro), free(riIcro);
    free(rrIIcro), free(riIIcro);
    free(lt_gb_se);
    free(lt_ea);
    printf("----------------------------------------\n");
    printf(" 2DIR calculation succesfully completed\n");
    printf("----------------------------------------\n\n");
}
