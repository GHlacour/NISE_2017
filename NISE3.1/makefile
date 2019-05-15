# Makefile

### Define the variable MACHINE in your environment
### export MACHINE=$HOSTNAME may be a good option
### If not define the default settings below will be used

### Default settings
ifeq ($(MACHINE), )
CC = gcc
FFTW3=-lfftw3
FFTW3L = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
FFTW3I = -fopenmp
CC_FLAGS = -B -O3 -lm -fopenmp
LAPACK = -llapack -lblas -lpthread
endif

ifeq ($(MACHINE), peri-gpu)
CC = pgcc
FFTW3=-lfftw3
FFTW3L = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
FFTW3I =
CC_FLAGS = -B -O3 -lm -tp sandybridge-64 -acc -ta=nvidia -Minfo=accel -fast -mp -Minfo
LAPACK = -lrt -lpgftnrtl -lm -lpgc -lgcc -llapack -lblas
endif

ifeq ($(MACHINE), peri)
CC = icc
FFTW3=-lfftw3
FFTW3L = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
FFTW3I = -qopenmp
CC_FLAGS = -B -O3 -lm -qopenmp
LAPACK = -llapack -lblas -lpthread
endif

### My settings
ifeq ($(MACHINE), mymachine)
CC=        # Your compiler name
LAPACK=    # Libraries for lapack
LAPACKL=   # Location of liblapack
FFTW3=     # Libraries for fftw3
FFTW3L=    # Location of libfftw3
FFTW3I=    # Location of fftw3.h file
CC_FLAGS=  # Optional compiler flags
endif

### My settings
ifeq ($(MACHINE), osx9)
CC=clang 
LAPACK=-framework Accelerate
LAPACKL=
FFTW3= -lfftw3
FFTW3L=-L/usr/local/lib
FFTW3I=-I/usr/local/include
CC_FLAGS=-O3 -lm 
endif

### My settings
ifeq ($(MACHINE), osx)
CC=gcc
LAPACK=-framework veclib    
LAPACKL=
FFTW3= -lfftw3
FFTW3L= -L/usr/local/lib   
FFTW3I=    
CC_FLAGS= -O3 -lm
endif

### My settings
ifeq ($(MACHINE), baiz)
CC=pgcc
LAPACK=-llapack -lblas -pgf77libs    
LAPACKL=/opt/pgi/linux86-64/10.8/lib/
FFTW3= -lfftw3
FFTW3L= -L/usr/lib   
FFTW3I=    
CC_FLAGS= -O3 -lm
endif

### Compiler settings medusa
ifeq ($(MACHINE), medusa)
CC=/usr/local/intel/cc/9.0/bin/icc
LAPACK=-lmkl_lapack -lmkl_ia32 -lguide -lpthread
LAPACKL=-L/usr/local/intel/mkl/8.0.1/lib/32
FFTW3=-lfftw3
FFTW3L=-L/home/lacour/fftw3/lib
FFTW3I=-I/home/lacour/fftw3/include
CC_FLAGS= -O3 -lm
endif

### Compiler settings new medusa
ifeq ($(MACHINE), newmedusa)
CC=icc
LAPACK=-lmkl_intel_lp64 -lmkl_lapack -lmkl_core -lguide -lmkl_intel_thread -lpthread 
LAPACKL=-L/usr/local/intel/mkl10.2/lib/em64t
FFTW3=-lfftw3f
FFTW3L=-L/home/lacour/fftw312/lib
FFTW3I=-I/home/lacour/fftw312/include
CC_FLAGS= -O3 -lm
endif

### Compiler Settings Hydron
ifeq ($(MACHINE), hydron)
CC = pgcc
LAPACK_DIR=/usr/local/pgi/linux86-64/7.0-7/lib/
CC_FLAGS = -O3 -L/usr/local/lib -lfftw3 -lm
LAPACK = -llapack -lblas -pgf77libs
LAPACKL = -L$(LAPACK_DIR) 
endif

### Compiler Settings HPCIBM
ifeq ($(MACHINE), hpcibm)
CC = /opt/pgi/linux86/6.0/bin/pgcc
CPPFLAGS= -B -I/opt/fftw-3.1.1/include -L/opt/fftw-3.1.1/lib
CC_FLAGS = -B -O3 -L/opt/fftw-3.1.1/lib -lfftw3 -lm
LAPACK = -llapack -lblas -lpgftnrtl
endif

### Compiler Settings Millipede
ifeq ($(MACHINE), millipede)
CC = icc
CPPFLAGS= -B -I/cm/shared/apps/fftw/intel/64/3.2.2/include -L/cm/shared/apps/fftw/intel/64/3.2.2/lib
CC_FLAGS = -B -O3 -L/cm/shared/apps/fftw/intel/64/3.2.2/lib -lfftw3 -lm
LAPACK = -llapack -lblas -lifcore
LAPACKL = -L$(LAPACK_DIR) -L$(BLASDIR)
endif

### Compiler Settings Bohr
ifeq ($(MACHINE), bohr)
CC = icc
FFTW3=-lfftw3
FFTW3L=-L/software/apps/fftw/intel/3.3.3/lib -L/software/apps/intel/compiler/2013/13.0.1/compiler/lib/intel64 -L/software/apps/intel/compiler/2013/13.0.1/mkl/lib/intel64
FFTW3I=-I/software/apps/fftw/intel/3.3.3/include -qopenmp
CC_FLAGS = -B -O3 -lm -qopenmp
LAPACK = -mkl -lmkl_intel_lp64 -lmkl_core -liomp5 -lmkl_intel_thread -lpthread
endif

### Compiler Settings Nieuwpoort
ifeq ($(MACHINE), nieuwpoort)
CC = icc
FFTW3=-lfftw3
FFTW3L=-L$(FFTWDIR) -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64
FFTW3I=-I$(FFTWINCLUDE) -qopenmp
CC_FLAGS = -B -O3 -lm -qopenmp
LAPACK = -mkl -lmkl_intel_lp64 -lmkl_core -liomp5 -lmkl_intel_thread -lpthread
endif

### Compiler Settings newNieuwpoort
ifeq ($(MACHINE), newnieuwpoortgcc)
CC = gcc
FFTW3=-lfftw3
FFTW3L = -L/cm/shared/apps/intel/compilers_and_libraries/linux/lib/intel64 
FFTW3I =
CC_FLAGS = -B -O3 -lm -fopenmp
LAPACK = -llapack -lblas -liomp5 -lpthread
endif


CPPFLAGS= $(FFTW3I)
LIBS=$(LAPACKL) $(LAPACK) $(FFTW3L) $(FFTW3)
#CC_FLAGS=-B -O3 -lm $(LAPACKL) $(LAPACK) $(FFTW3L) $(FFTW3)

all: echo NISE3.1 translate 

echo:
	@echo Compiling on the machine $(MACHINE) with $(CC)

NISE3.1: NISE3.1.o NISE3.1subs.o polar.o 1DFFT.o absorption.o c_absorption.o calc_2DIR.o calc_2DES.o luminescence.o calc_CD.o analyse.o readinput.o randomlib.o
	$(CC) NISE3.1.o NISE3.1subs.o polar.o 1DFFT.o absorption.o c_absorption.o calc_2DIR.o calc_2DES.o luminescence.o calc_CD.o analyse.o readinput.o randomlib.o $(LIBS) -o NISE3.1 $(CC_FLAGS)

NISE3.1.o: NISE3.1.c NISE3.1.h types3.1.h lapack.h

polar.o: polar.c polar.h 

NISE3.1subs.o: NISE3.1subs.c NISE3.1subs.h types3.1.h lapack.h

readinput.o: readinput.c readinput.h types3.1.h

1DFFT.o: 1DFFT.c 1DFFT.h types3.1.h

absorption.o: absorption.c types3.1.h lapack.h NISE3.1subs.h 1DFFT.h

c_absorption.o: c_absorption.c types3.1.h lapack.h NISE3.1subs.h 1DFFT.h

calc_2DIR.o: calc_2DIR.c types3.1.h polar.h lapack.h NISE3.1subs.h

calc_2DES.o: calc_2DES.c types3.1.h polar.h lapack.h NISE3.1subs.h

luminescence.o: luminescence.c absorption.h types3.1.h lapack.h NISE3.1subs.h

calc_CD.o: calc_CD.c calc_CD.h types3.1.h lapack.h NISE3.1subs.h

analyse.o: analyse.c types3.1.h lapack.h NISE3.1subs.h

translate.o: translate.c translate.h

translate: translate.o NISE3.1subs.o readinput.o randomlib.o
	$(CC) translate.o NISE3.1subs.o readinput.o randomlib.o $(LIBS) -o translate $(CC_FLAGS)


nonsubs.o: nonsubs.c nonsubs.h types.h lapack.h

2DFFT: 2DFFT.o nrutil.o
	$(CC) nrutil.o 2DFFT.o $(LIBS) -o 2DFFT $(CC_FLAGS) 

2DFFT.o: 2DFFT.c

nrutil.o: nrutil.h nrutil.c

randomlib.o: randomlib.h randomlib.c

clean:
	\rm *.o
