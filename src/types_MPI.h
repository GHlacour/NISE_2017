#ifndef _TYPES_MPI_
#define _TYPES_MPI_
#include "types.h"
#include <stddef.h>
#include <mpi.h>

// Struct macro for an MPI data type
// Macro to work around the VLA (var. len. arrays) limitations in older compilers
#define CUSTOM_MPI_DATATYPE(LEN) struct { \
    int length;\
    int blocklengths[LEN];\
    MPI_Datatype types[LEN];\
    MPI_Aint offsets[LEN];\
}

typedef CUSTOM_MPI_DATATYPE(62) t_non_datatype;
extern const t_non_datatype T_NON_TYPE;

#endif
