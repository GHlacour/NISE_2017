// Just here to not confuse IDEs with the includes, the header file of the compiled
// version in ./src is actually used.

#ifndef _RANDOM_
#define _RANDOM_

void   RandomInitialise(int,int);
double RandomUniform(void);
double RandomGaussian(double,double);
int    RandomInt(int,int);
double RandomDouble(double,double);

#endif // _RANDOM_
