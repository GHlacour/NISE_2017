#ifndef _TYPES_
#define _TYPES_

// Commonly used constants
float static k_B=0.6950389 ;/* cm-1K-1 */
float static c_v=2.99792458e-5;/* Speed of light in cm/fs */
float static twoPi=2*3.14159265359;
float static icm2ifs=2.99792458e-5;
float static ifs2icm=1.0/2.99792458e-5;
float static sqrt2=1.41421356237;
float static ithirty=1.0/30;

// Structures
typedef struct {
    int tmax1, tmax2, tmax3;
    int tmin1, tmin2, tmin3;
    int dt[4];
    int d1, d2, d3;
} t_time;

typedef struct {
    float min1, max1;
    float min2, max2;
    float min3, max3;
    int k[4];
} t_rwa;

typedef struct {
    int length;
    int sample;
    int begin;
    int end;
} t_steps;

// REMINDER Change the lengths/types for the MPI data in types_MPI.c when changing this model!
typedef struct {
  int tmax1,tmax2,tmax3;
  int tmin1,tmin2,tmin3;
  int dt1,dt2,dt3;
  int d1,d2,d3;
  float min1,max1;
  float min2,max2;
  float min3,max3;
  int k[4];
  int length;
  int sample;
  int fft;
  int begin;
  int end;
  int is;
  int ts;
  int interpol;
  int propagation;
  float lifetime;
  float deltat;
  int buffer,singles,doubles;
  char energyFName[256];
  char dipoleFName[256];
  char alphaFName[256];
  char positionFName[256];
  char anharFName[256];
  char overdipFName[256];
  char technique[256];
  char basis[256];
  char pdbFName[256];
  char hamiltonian[256];
  char couplingFName[256];
  int tmax;
  int cluster;
  float shifte,shiftf;
  int labPol; // 0=parallel, 1=perpendicular
  float dephasing;
  float statstart,statend,statstep;
  int statsteps;
  float thres;
  float couplingcut;
  float temperature;
  float anharmonicity;
  int Npsites;
  int *psites;
} t_non;

#endif // _TYPES_
