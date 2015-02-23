# compile flags
SDL=1
OPTIMIZATION_LEVEL=3

# object files
OBJ=fpu.o

# compilers
CC=gcc

# compiler options
ifeq ($(SDL),1)
	OPTS=-Wall -pthread `sdl-config --cflags`
else
	OPTS=-Wall -pthread
endif

ifeq ($(SDL),1)
	CFLAGS=-O$(OPTIMIZATION_LEVEL) -DSDL=$(SDL) -lm -lrt -lsndfile -ltiff `sdl-config --libs`
else
	CFLAGS=-O$(OPTIMIZATION_LEVEL) -lm -lrt -lsndfile
endif

fpu: 	$(OBJ)
	$(CC) -o $@ $+ $(OPTS) $(CFLAGS) 

fpu.o:	fpu.c
	$(CC) $(OPTS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm *.o fpu

