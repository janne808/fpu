# compile flags
SDL=1
WAV_ENABLE=0
TIFF_ENABLE=0
OPTIMIZATION_LEVEL=3

# object files
OBJ=fpu.o

# compilers
CC=gcc

# SDL options
ifeq ($(SDL),1)
	SDL_OPTS=-DSDL=1 `sdl-config --libs` `sdl-config --cflags`
else
	SDL_OPTS=
endif

# sndfile options
ifeq ($(WAV_ENABLE),1)
	SNDFILE_OPTS=-DWAV_ENABLE=1 -lsndfile
else
	SNDFILE_OPTS=
endif

# TIFFLib options
ifeq ($(TIFF_ENABLE),1)
	TIFFLIB_OPTS=-DTIFF_ENABLE=1 -ltiff
else
	TIFFLIB_OPTS=
endif

# compiler options
OPTS=-Wall -pthread
CFLAGS=-O$(OPTIMIZATION_LEVEL) -lm -lrt

fpu: 	$(OBJ)
	$(CC) -o $@ $+ $(OPTS) $(CFLAGS) $(SDL_OPTS) $(SNDFILE_OPTS) $(TIFFLIB_OPTS)

fpu.o:	fpu.c
	$(CC) $(OPTS) $(CFLAGS) $(SDL_OPTS) $(SNDFILE_OPTS) $(TIFFLIB_OPTS) -c $<

.PHONY: clean
clean:
	rm *.o fpu

