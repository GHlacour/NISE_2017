#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "NISE_subs.h"
#include "projection.h"
#include "randomlib.h"
#include "util/asprintf.h"

// This subroutine nullify all elements of vector for non selected sites
void projection(float* phi, t_non* non) {
    int i;
    for (i = 0; i < non->singles; i++) {
        phi[i] = phi[i] * non->psites[i];
    }
    return;
}
