# Makefile

CC      = gcc

CC_FLAGS= -O3 -lm 

all: stochastic stochasticSFG

stochastic: stochastic.o randomlib.o
	$(CC) $(CC_FLAGS) randomlib.o stochastic.o -o stochastic

stochastic.o: stochastic.c stochastic.h

stochasticSFG: stochasticSFG.o randomlib.o
	$(CC) $(CC_FLAGS) randomlib.o stochasticSFG.o -o stochasticSFG

stochasticSFG.o: stochasticSFG.c stochasticSFG.h

randomlib.o: randomlib.h randomlib.c
