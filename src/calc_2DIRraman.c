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
#include "polar.h"
#include "calc_2DIRraman.h"
#include <stdarg.h>
#include "mpi.h"
#include "MPI_subs.h"
#include "project.h"

void calc_2DIRraman(t_non* non, int parentRank, int parentSize, int subRank, int subSize, MPI_Comm subComm, MPI_Comm rootComm) {
    // Start by determining the work to be done, make an array of samples/poldir to simulate
    // and distribute those statically over all processors
    int clusterCount = 0, sampleCount = 0, * workset, * worksetSizes = calloc(parentSize, sizeof(int));;
    int counter,counter_pass;
    float counter_current;
    double my_time,my_current_time;
    counter=0,counter_pass=1;
    if(parentRank == 0) {
        // Master process calculates the work items to be performed
        int* fullWorkset;
        calculateWorkset(non, &fullWorkset, &sampleCount, &clusterCount);
        log_item("Begin sample: %d, End sample: %d.\n", non->begin, non->end);

        // Distribute work, each process does as many items as the others (static decomposition)
        MPI_Bcast(&clusterCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&sampleCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int totalWorkItems = 21 * sampleCount;
        int baseWorksetSize = totalWorkItems / parentSize;
        int remainder = totalWorkItems % parentSize;

        // Make array of worksetsizes, since each process might have a different workset size, distribute that array to the other processes
        int* worksetOffsets = calloc(parentSize, sizeof(int));
        for(int i = 0; i < parentSize; i++) {
            worksetSizes[i] = i < remainder ? (baseWorksetSize + 1) * 2 : baseWorksetSize * 2;
            worksetOffsets[i] = i == 0 ? 0 : worksetOffsets[i - 1] + worksetSizes[i];
        }

        MPI_Bcast(worksetSizes, parentSize, MPI_INT, 0, MPI_COMM_WORLD);
        workset = malloc(worksetSizes[0] * sizeof(int));
        MPI_Scatterv(fullWorkset, worksetSizes, worksetOffsets, MPI_INT, workset, worksetSizes[0], MPI_INT, 0, MPI_COMM_WORLD);

        free(fullWorkset);
        free(worksetOffsets);
    } else {
        // Fix non settings
        const int totalSampleCount = (non->length - non->tmax1 - non->tmax2 - non->tmax3 - 1) / non->sample + 1;
        if (non->end == 0) non->end = totalSampleCount;

        // Receive information on workset
        MPI_Bcast(&clusterCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&sampleCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(worksetSizes, parentSize, MPI_INT, 0, MPI_COMM_WORLD);
        workset = malloc(worksetSizes[parentRank] * sizeof(int));
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT, workset, worksetSizes[parentRank], MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Initialize each process base variables
    float* pol = 0; /* Currently dummy vector that can be used to change coordinate system in the future */ // RO
    const int nn2 = non->singles * (non->singles + 1) / 2;
    float* Hamil_i_e = calloc(nn2, sizeof(float));

    // Frequency shifts
    float shift1 = (non->max1 + non->min1) / 2;
    if(parentRank == 0) printf("Rotating Wave Frequency shift %f.\n", shift1);
    non->shifte = shift1;
    non->shiftf = 2 * shift1;

    // Arrays where the result is stored, these will be reduced (summed) at the end!

    // 2D response function parallel
    float** rrIpar = (float**) calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIpar = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** rrIIpar = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIIpar = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    // 2D response function perpendicular
    float** rrIper = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIper = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** rrIIper = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIIper = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    // 2D response function cross
    float** rrIcro = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIcro = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** rrIIcro = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE
    float** riIIcro = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // REDUCE

    // These arrays are initialized here and only read in the loops
    float** lt_gb_se = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // RO
    float** lt_ea = (float**)calloc2D(non->tmax3, non->tmax1, sizeof(float), sizeof(float*)); // RO

    for (int t1 = 0; t1 < non->tmax1; t1++) {
        for (int t3 = 0; t3 < non->tmax3; t3++) {
            lt_gb_se[t3][t1] = (float) exp(-(double)(t1 + t3) * non->deltat / (2 * non->lifetime));
            lt_ea[t3][t1] = (float) exp(-(double)(t1 + t3) * non->deltat / (2 * non->lifetime));
        }
    }

    // Open the various files, these are readonly per process so safe for concurrent use
    FILE* H_traj = fopen(non->energyFName, "rb");
    if (H_traj == NULL) {
        if(parentRank == 0) printf("Hamiltonian file not found!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    FILE* mu_traj = fopen(non->dipoleFName, "rb");
    if (mu_traj == NULL) {
        if (parentRank == 0) printf("Dipole file %s not found!\n", non->dipoleFName);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    FILE* alpha_traj = fopen(non->alphaFName, "rb");
    if (alpha_traj == NULL) {
        if (parentRank == 0) printf("Raman file %s not found!\n", non->alphaFName);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //ToDo: Raman -- anharmonicities

    /* Open file for fluctuating anharmonicities and sequence transition dipoles if needed */
    FILE* A_traj, *mu2_traj;
    if (non->anharmonicity == 0) {
        A_traj = fopen(non->anharFName, "rb");
        if (A_traj == NULL) {
            if (parentRank == 0) printf("Anharmonicity file %s not found!\n", non->anharFName);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        mu2_traj = fopen(non->overdipFName, "rb");
        if (mu2_traj == NULL) {
            if (parentRank == 0) printf("Overtone dipole file %s not found!\n", non->overdipFName);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* Read coupling */
    float* mu_xyz = calloc(3 * non->singles, sizeof(float));
    if (!strcmp(non->hamiltonian, "Coupling")) {
        FILE* C_traj = fopen(non->couplingFName, "rb");
        if (C_traj == NULL) {
            if (parentRank == 0) printf("Coupling file not found!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (read_He(non, Hamil_i_e, C_traj, -1) != 1) {
            if (parentRank == 0) printf("Coupling trajectory file to short, could not fill buffer!!!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fclose(C_traj);

        for (int x = 0; x < 3; x++) {
            if (read_mue(non, mu_xyz + non->singles * x, mu_traj, 0, x) != 1) {
                if (parentRank == 0) printf("Dipole trajectory file to short, could not fill buffer!!!\n");
                if (parentRank == 0) printf("ITIME %d %d\n", 0, x);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // Start clock
    if (parentRank==0){
	my_time=MPI_Wtime();
    }

    // From now on we'll do the calculations
    for (int currentWorkItem = 0; currentWorkItem < worksetSizes[parentRank]; currentWorkItem += 2) {
        int currentSample = workset[currentWorkItem];
        int molPol = workset[currentWorkItem + 1];

        if (currentSample == -1 || molPol == -1) continue;

        /* Calculate 2DIR-Raman response */
        int tj = currentSample * non->sample + non->tmax1;
        int tk = tj + non->tmax2;
        int px[4];
        int alphapol; //variable for Raman polarization from px[2] px[3] combination
        polar(px, molPol);

        // Allocate arrays
        float* Anh = calloc(non->singles, sizeof(float));
        float* over = calloc(non->singles, sizeof(float));

        float* leftrr = calloc(non->singles, sizeof(float));
        float* leftri = calloc(non->singles, sizeof(float));
        float** leftnr = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float** leftni = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float** rightrr = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float** rightri = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float* rightnr = calloc(non->singles, sizeof(float));
        float* rightni = calloc(non->singles, sizeof(float));
        float** leftn3r = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float** leftn3i = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));

        float* mut3r = calloc(non->singles, sizeof(float)); // Compared to 2DIR: now dipole of px[1] at tj
        float* alpha4 = calloc(non->singles, sizeof(float));

        float* fr = calloc(nn2, sizeof(float));
        float* fi = calloc(nn2, sizeof(float));
        float** ft1r = (float**)calloc2D(non->tmax1, nn2, sizeof(float), sizeof(float*));
        float** ft1i = (float**)calloc2D(non->tmax1, nn2, sizeof(float), sizeof(float*));

        /* First excitation kI and kII */
        for (int t1 = 0; t1 < non->tmax1; t1++) {
            /* Read dipoles at time 0, single excitation for kI */
            mureadE(non, leftnr[t1], tj - t1, px[0], mu_traj, mu_xyz, pol);
            clearvec(leftni[t1], non->singles);

            /*Double excitation for kII */
            if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2"))||
            (!strcmp(non->technique, "2DIRraman3")) || (!strcmp(non->technique, "2DIRramanII"))) {
                //! Generates double excited state (ft1r) from single excited state (leftnr) via: mu_eg = sqrt(2) mu_fg
                //! ft1i is 0's due to excitation from the ground state
                dipole_double_ground(non, leftnr[t1],ft1r[t1], ft1i[t1],over);
            }

            /* Propagate */
            for (int tm = 0; tm < t1; tm++) {
                /* Read Hamiltonian */
                if (read_He(non, Hamil_i_e, H_traj, tj - t1 + tm) != 1) {
                    printf("Hamiltonian trajectory file to short, could not fill buffer!!!\n");
                    exit(1);
                }

                /* Propagation of singles for kI*/
                if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman1"))|| 
                    (!strcmp(non->technique, "2DIRramanI"))) {
                    if (non->propagation == 1) {
                        propagate_vec_coupling_S(non, Hamil_i_e, leftnr[t1], leftni[t1], non->ts, 1);
                    } else if (non->propagation == 0) {
                        propagate_vec_DIA_S(non, Hamil_i_e, leftnr[t1], leftni[t1], 1);
                    }
                }

                /* Propagation of doubles for kII */
                if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2")) ||
                    (!strcmp(non->technique, "2DIRraman3"))|| (!strcmp(non->technique, "2DIRramanII"))) {

                    if (non->anharmonicity == 0) {
                        printf("Anharmonicity feature for Raman not implemented yet");
                        exit(1);
                        read_A(non, Anh, A_traj, tj-t1+tm);
                    }

                    if (non->propagation == 1) {
                        propagate_vec_coupling_S_doubles(non, Hamil_i_e, ft1r[t1], ft1i[t1], non->ts, Anh);
                    } else if (non->propagation == 0) {

                        float* Urs = calloc(non->singles * non->singles, sizeof(float));
                        float* Uis = calloc(non->singles * non->singles, sizeof(float));

                        int* Rs = calloc(non->singles * non->singles, sizeof(int));
                        int* Cs = calloc(non->singles * non->singles, sizeof(int));

                        int elements = time_evolution_mat(non, Hamil_i_e, Urs, Uis, Cs, Rs, non->ts);

                        propagate_double_sparce(non, Urs, Uis, Cs, Rs, ft1r[t1], ft1i[t1],elements,non->ts, Anh);

                        free(Urs), free(Uis), free(Rs), free(Cs);
                    }
                }

            }
        }

        /* Read dipole for second interaction */
        mureadE(non, mut3r, tj, px[1], mu_traj, mu_xyz, pol);

        if (non->anharmonicity == 0) {
            read_over(non, over, mu2_traj, tj, px[1]);
        }

        /*Double excitation for kI */
        if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman1")) 
            || (!strcmp(non->technique, "2DIRramanI"))) {
            dipole_double_ground(non, mut3r, fr, fi,over);
            //! Transpose dipoles of t1 from left-side array to right-side array to align with Diagram 1
            memcpy(rightrr[0], leftnr[0], non->tmax1 * non->singles * sizeof(float));
            memcpy(rightri[0], leftni[0], non->tmax1* non->singles * sizeof(float));
            for (int i = 0; i < non->tmax1 * non->singles; i++) rightri[0][i] = -rightri[0][i];
        }

        /* Excitation on right side for Diagram 2 (kII)*/
        if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2"))||
            (!strcmp(non->technique, "2DIRramanII"))) {
            //! Transpose dipoles of t3 to right-side array to align with Diagram 2
            //! rightni is zero since mut3 has no imaginary part, due to excitation from the ground state
            memcpy(rightnr, mut3r, non->singles * sizeof(float));
            for (int i = 0; i < non->singles; i++) rightni[i] = 0;
        }

        float** t1nr = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));
        float** t1ni = (float**)calloc2D(non->tmax1, non->singles, sizeof(float), sizeof(float*));

        /* De-excitation on left side for Diagram 3 (kII)*/
        if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman3"))||
            (!strcmp(non->technique, "2DIRramanII"))) {
            for(int t1 = 0; t1 < non->tmax1; t1++) {
            //! Multiply double excited states (ft1) with double exciton dipole mut3r, new state is t1n
            dipole_double_last(non, mut3r, ft1r[t1], ft1i[t1],t1nr[t1], t1ni[t1],over);
            }
        }

        /* Combine with evolution during t3 */
        for (int t3 = 0; t3 < non->tmax3; t3++) {
            int tl = tj + t3;

            /* Determine polarization of the Raman Tensor*/
            find_raman_pol(px[2],px[3],alphapol);

            /* Read Raman Tensor t3 */
            read_alpha(non, alpha4, alpha_traj, tl, alphapol);

            if (non->anharmonicity == 0) {
                // anharmonicity for raman is different from dipole due to 6 in stead of 3 terms
                printf("Anharmonicity feature for Raman not implemented yet");
                exit(1);
                read_over(non, over, mu2_traj, tl, px[3]);
            }

            /* Multiply with the Raman tensor and Calculate Diagram 1 response */
            if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman1"))
                || (!strcmp(non->technique, "2DIRramanI"))) {
                //! since step at t1 (g -> e) is on right side, no t1 dependence for leftrr/leftri
                dipole_double_last(non, alpha4, fr, fi, leftrr, leftri,over);
                for (int t1 = 0; t1 < non->tmax1; t1++) {
                    float rrI = 0, riI = 0;
                    for (int i = 0; i < non->singles; i++) {
                        rrI += leftri[i] * rightrr[t1][i] + leftrr[i] * rightri[t1][i];
                        riI += leftrr[i] * rightrr[t1][i] - rightri[t1][i] * leftri[i];
                    }

                    float polWeight = polarweight(0, molPol) * lt_ea[t3][t1];
                    rrIpar[t3][t1] += rrI * polWeight;
                    riIpar[t3][t1] += riI * polWeight;
                    polWeight = polarweight(1, molPol) * lt_ea[t3][t1];
                    rrIper[t3][t1] += rrI * polWeight;
                    riIper[t3][t1] += riI * polWeight;
                    polWeight = polarweight(2, molPol) * lt_ea[t3][t1];
                    rrIcro[t3][t1] += rrI * polWeight;
                    riIcro[t3][t1] += riI * polWeight;
                }
            }

            /* Multiply with the Raman tensor and Calculate Diagram 2 response */
            if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2"))
                || (!strcmp(non->technique, "2DIRramanII"))) {
                for (int t1 = 0; t1 < non->tmax1; t1++) {
                    //! since step at t1 (g -> f) is on left side, t1 dependence for leftnr/leftni
                    dipole_double_last(non, alpha4, ft1r[t1], ft1i[t1], leftnr[t1], leftni[t1],over);

                    float rrII = 0, riII = 0;
                    for (int i = 0; i < non->singles; i++) {
                        rrII += rightnr[i] * leftni[t1][i] + rightni[i] * leftnr[t1][i];
                        riII += rightnr[i] * leftnr[t1][i] - rightni[i] * leftni[t1][i];
                    }

                    float polWeight = polarweight(0, molPol) * lt_ea[t3][t1];
                    rrIIpar[t3][t1] += rrII * polWeight;
                    riIIpar[t3][t1] += riII * polWeight;
                    polWeight = polarweight(1, molPol) * lt_ea[t3][t1];
                    rrIIper[t3][t1] += rrII * polWeight;
                    riIIper[t3][t1] += riII * polWeight;
                    polWeight = polarweight(2, molPol) * lt_ea[t3][t1];
                    rrIIcro[t3][t1] += rrII * polWeight;
                    riIIcro[t3][t1] += riII * polWeight;
                }
            }

            /* Multiply with the Raman tensor and Calculate Diagram 3 response */
            if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman3"))
                || (!strcmp(non->technique, "2DIRramanII"))) {
                for (int t1 = 0; t1 < non->tmax1; t1++) {
                    //! since step at t1 (g -> f) is on left side, t1 dependence for leftnr/leftni

                   //! empty leftnr/leftni since it was already filled during t1 for use in (g -> f)
                   //! with dipole_double_ground ft1r/ft1i arrays were calculated with leftnr/leftni,
                   //! ft1r/ft1i arrays were used to calculate t1nr/t1ni
                   for(int i = 0; i<non->singles; i++){
                       leftn3r[t1][i] = 0;
                       leftn3i[t1][i] = 0;
                   }

                   for(int i =0; i <non->singles; i++){
                    leftn3r[t1][i] += alpha4[i] * t1nr[t1][i];
                    leftn3i[t1][i] += alpha4[i] * t1ni[t1][i];
                    }

                    float rrII = 0, riII = 0;

                    int rightnr3 = 1, rightni3 = 0; //! since the right side is in the ground state
                    for (int i = 0; i < non->singles; i++) {
                        rrII += rightnr3 * leftn3i[t1][i] + rightni3 * leftn3r[t1][i];
                        riII += rightnr3 * leftn3r[t1][i] - rightni3 * leftn3i[t1][i];
                    }

                    //! -= is used since this pahtway has opposite sign to the other pathways
                    float polWeight = polarweight(0, molPol) * lt_ea[t3][t1];
                    rrIIpar[t3][t1] -= rrII * polWeight;
                    riIIpar[t3][t1] -= riII * polWeight;
                    polWeight = polarweight(1, molPol) * lt_ea[t3][t1];
                    rrIIper[t3][t1] -= rrII * polWeight;
                    riIIper[t3][t1] -= riII * polWeight;
                    polWeight = polarweight(2, molPol) * lt_ea[t3][t1];
                    rrIIcro[t3][t1] -= rrII * polWeight;
                    riIIcro[t3][t1] -= riII * polWeight;
                }
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

            /* Propagate:
            The doubles on the left-side (fr/fi) and the single excitation on the right-side (rightr) for Diagram 1
            The doubles on the left-side (ft1) and the single excitation on the right-side (rightn) for Diagram 2
            The single excitation on the left-side (t1n) only for Diagram 3 */
            if (non->propagation == 0) {
                propagate_0(non, Hamil_i_e, currentSample, molPol, t3, fr, fi, ft1r, ft1i, t1nr, t1ni, rightnr, rightni, rightrr, rightri, Anh);
            }
            else if(non->propagation == 1) {
                propagate_1(non, Hamil_i_e, fr, fi, ft1r, ft1i, t1nr, t1ni, rightnr, rightni, rightrr, rightri, Anh);
            }
        }


        free(leftrr), free(leftri), free2D((void**) leftnr), free2D((void**) leftni),free2D((void**) leftn3r), free2D((void**) leftn3i);
        free2D((void**) rightrr), free2D((void**) rightri), free(rightnr), free(rightni);
        free2D((void**) t1nr), free2D((void**) t1ni);
        free(mut3r);
        free(alpha4);
        free(Anh), free(over);
        free(fr), free(fi);
        free2D((void**) ft1r), free2D((void**) ft1i);

        counter++;
	if (parentRank==0){
	    counter_current=counter*100.0/sampleCount/21*parentSize;
            if (counter_current>counter_pass){
		if (non->printLevel>0){
                    my_current_time=MPI_Wtime();
                    char* timeText=MPI_time(my_current_time-my_time);
  		    printf("Passed %d pct. of expected calculation time in %s",(int)floor(counter_current),timeText);
                    fflush(stdout);
                    if (non->printLevel==1){
                      counter_pass=floor(counter_current)+10;
                    } else {
                      counter_pass=floor(counter_current)+1;
                    }
		    free(timeText);
                }
	    }
	}
    }

    // Reduce calculation results, in a tiered approach to save network bandwidth
    int reduceArraySize = non->tmax3 * non->tmax1;
    float** reductionArrays[12] = {
        rrIpar, riIpar, rrIIpar, riIIpar, rrIper, riIper, rrIIper, riIIper, rrIcro, riIcro, rrIIcro, riIIcro
    };
    MPI_Request reductions[2][12];

    // Reduce at small scope
    if(subRank == 0) { // parent will do in-place reduce to save memory
        for (int i = 0; i < 12; i++) {
            MPI_Ireduce(MPI_IN_PLACE, reductionArrays[i][0], reduceArraySize, MPI_FLOAT, MPI_SUM, 0, subComm, &(reductions[0][i]));
        }
    } else {
        for (int i = 0; i < 12; i++) {
            MPI_Ireduce(reductionArrays[i][0], NULL, reduceArraySize, MPI_FLOAT, MPI_SUM, 0, subComm, &(reductions[0][i]));
        }
    }

    // Wait for the small-scope reduction to finish
    asyncWaitForMPI(reductions[0], 12, 1, 5000);

    // Reduce at global scope
    if (parentRank == 0) { // parent will do in-place reduce to save memory
        for (int i = 0; i < 12; i++) {
            MPI_Ireduce(MPI_IN_PLACE, reductionArrays[i][0], reduceArraySize, MPI_FLOAT, MPI_SUM, 0, rootComm, &(reductions[1][i]));
        }
    }
    else if(subRank == 0) {
        for (int i = 0; i < 12; i++) {
            MPI_Ireduce(reductionArrays[i][0], NULL, reduceArraySize, MPI_FLOAT, MPI_SUM, 0, rootComm, &(reductions[1][i]));
        }
    }

    // Wait for the global scope reduction to finish
    if(parentRank == 0 || subRank == 0)
        asyncWaitForMPI(reductions[1], 12, 1, 5000);

    /* The calculation is finished, lets write output */
    if (parentRank == 0) {
        log_item("Finished Calculating Response!\nWriting to file\n");

        printf("Samples %d\n", sampleCount);
        if (non->cluster != -1) {
            printf("Of %d samples %d belonged to cluster %d.\n", sampleCount, clusterCount, non->cluster);
        }

        /* Close Files */
        fclose(mu_traj), fclose(H_traj),fclose(alpha_traj);

        if (non->anharmonicity == 0) {
            fclose(mu2_traj), fclose(A_traj);
        }


        /* Print 2D */
            print2D("RparIraman.dat", rrIpar, riIpar, non, sampleCount);
            print2D("RperIraman.dat", rrIper, riIper, non, sampleCount);
            print2D("RcroIraman.dat", rrIcro, riIcro, non, sampleCount);

            print2D("RparIIraman.dat", rrIIpar, riIIpar, non, sampleCount);
            print2D("RperIIraman.dat", rrIIper, riIIper, non, sampleCount);
            print2D("RcroIIraman.dat", rrIIcro, riIIcro, non, sampleCount);

        printf("----------------------------------------\n");
        printf(" 2DIR-Raman calculation succesfully completed\n");
        my_current_time=MPI_Wtime();
        char* timeText=MPI_time(my_current_time-my_time);
        printf(" Total time elapsed %s",timeText);
	free(timeText);
	printf("----------------------------------------\n\n");
    }

    /* Free memory for 2D calculation */
    free2D((void**)rrIpar), free2D((void**)riIpar);
    free2D((void**)rrIIpar), free2D((void**)riIIpar);
    free2D((void**)rrIper), free2D((void**)riIper);
    free2D((void**)rrIIper), free2D((void**)riIIper);
    free2D((void**)rrIcro), free2D((void**)riIcro);
    free2D((void**)rrIIcro), free2D((void**)riIIcro);
    free(mu_xyz);
    free2D((void**)lt_gb_se);
    free2D((void**)lt_ea);
    free(workset); free(worksetSizes);
    free(Hamil_i_e);

    // Barrier to gather all threads and exit simultaneously
    MPI_Request barrierRequest[1];
    MPI_Ibarrier(MPI_COMM_WORLD, &(barrierRequest[0]));
    asyncWaitForMPI(barrierRequest, 1, 1, 5000);
}

/*Propagation function for propagate scheme 0*/
void propagate_0(t_non* non, float* Hamil_i_e, int currentSample, int molPol, int t3, float* fr, float* fi,float** ft1r, float** ft1i, float** t1nr, float** t1ni, float* rightnr, float* rightni, float** rightrr, float** rightri,float* Anh){

    float* Urs = calloc(non->singles * non->singles, sizeof(float));
    float* Uis = calloc(non->singles * non->singles, sizeof(float));

    int* Rs = calloc(non->singles * non->singles, sizeof(int));
    int* Cs = calloc(non->singles * non->singles, sizeof(int));

    // bug? Cs and Rs seems to be exchanged (TLC just multiplying with imaginary number)
    int elements = time_evolution_mat(non, Hamil_i_e, Urs, Uis, Cs, Rs, non->ts);
    if (currentSample == non->begin && molPol == 0 && t3 == 0) {
        printf("Sparse matrix efficiency: %f pct.\n",
                (1 - (1.0 * elements / (non->singles * non->singles))) * 100);
        printf("Present truncation %f.\n",
                non->thres / ((double) non->deltat * icm2ifs * (double) twoPi / non->ts * (non->deltat *
                    icm2ifs * twoPi / non->ts)));
        printf("Suggested truncation %f.\n", 0.001);
    }

    // Key parallel loop 1
    // Initial step, former t1=-1
    propagate_double_sparce(
        non, Urs, Uis, Rs, Cs, fr, fi, elements, non->ts, Anh
    );

    int t1; // MSVC can't deal with C99 declarations inside a for with OpenMP
    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2"))
        || (!strcmp(non->technique, "2DIRramanII"))) {
        #pragma omp parallel for \
            shared(non, Anh, Urs, Uis, Rs, Cs, ft1r, ft1i) \
            schedule(static, 1)

        for(t1 = 0; t1 < non->tmax1; t1++) {
            propagate_double_sparce(
                non, Urs, Uis, Rs, Cs, ft1r[t1],
                ft1i[t1], elements, non->ts, Anh
            );
        }
    }

    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman3"))
        || (!strcmp(non->technique, "2DIRramanII"))) {
        for(t1 = 0; t1 < non->tmax1; t1++) {
            propagate_vec_DIA_S(non, Hamil_i_e, t1nr[t1], t1ni[t1], -1);
        }
    }
    // Propagate vectors right
    // Key parallel loop 2
    // Initial step

    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman1"))
        || (!strcmp(non->technique, "2DIRramanI"))) {
        propagate_vec_DIA_S(non, Hamil_i_e, rightnr, rightni, -1);

        for(t1 = 0; t1 < non->tmax1; t1++) {
            propagate_vec_DIA_S(non, Hamil_i_e, rightrr[t1], rightri[t1], -1);
        }
    }
    free(Urs), free(Uis), free(Rs), free(Cs);
}

/*Propagation function for propagation scheme 1*/
void propagate_1(t_non* non, float* Hamil_i_e, float* fr, float* fi,float** ft1r, float** ft1i,float** t1nr, float** t1ni, float* rightnr, float* rightni, float** rightrr, float** rightri,float* Anh){

    int t1;
    //Propagate left and right side diagram 2, during t3
    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman2"))
        || (!strcmp(non->technique, "2DIRramanII"))) {
        #pragma omp parallel for \
            shared(non,Hamil_i_e,Anh,ft1r,ft1i) \
            schedule(static, 1)

        for (t1 = 0; t1 < non->tmax1; t1++) {
            propagate_vec_coupling_S_doubles(non, Hamil_i_e, ft1r[t1], ft1i[t1], non->ts,Anh);
        }

        propagate_vec_coupling_S(non, Hamil_i_e, rightnr, rightni, non->ts, -1);
    }

    //Propagate left side diagram 3, during t3 (no right side propagation due to ground state)
    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman3"))
        || (!strcmp(non->technique, "2DIRramanII"))) {

        for (t1 = 0; t1 < non->tmax1; t1++) {
            //! sign +1 for left side needed
            propagate_vec_coupling_S(non, Hamil_i_e, t1nr[t1], t1ni[t1], non->ts, 1);
        }
    }

    //Propagate right and left side of diagram 1, during t3
    if ((!strcmp(non->technique, "2DIRraman")) || (!strcmp(non->technique, "2DIRraman1"))
        || (!strcmp(non->technique, "2DIRraman2"))) {
        propagate_vec_coupling_S_doubles(non, Hamil_i_e, fr, fi, non->ts,Anh);

        for (t1 = 0; t1 < non->tmax1; t1++) {
            propagate_vec_coupling_S(non, Hamil_i_e, rightrr[t1], rightri[t1], non->ts, -1);
        }
    }
}