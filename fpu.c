// fpu routines

/*
 *  (C) 2015 Janne Heikkarainen <janne.heikkarainen@student.tut.fi>
 *
 *  All rights reserved.
 *
 *  This file is part of 2-D Fermi-Pasta-Ulam solver.
 *
 *  Fpu is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Fpu is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Fpu.  If not, see <http://www.gnu.org/licenses/>.
 */

// wave writer enable
#define WAV_ENABLE 0

// tiff writer enable
#define TIFF_ENABLE 0

#include <pthread.h>

#if WAV_ENABLE
#include <sndfile.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <SDL.h>

#if TIFF_ENABLE
#include <tiffio.h>
#endif

#include "fpu.h"

// number of posix threads
#define NUM_THREADS 2

// draw every Nth frame
#define FRAMESKIP 2

// graphics scale up factor
#define SCALE 4

// minimum and maximum geometries
#define XMIN -1.0
#define XMAX 1.0

// number of oscillators per dimension
#define NUM 128

// time step factor
#define DT 8.0

// nonlinear factor
#define BETA 5.0

struct thread_data thread_data_array[NUM_THREADS];

#if TIFF_ENABLE
void writeframe(char* path, SDL_Surface *screen){
  TIFF *file;
  Uint8 *p;
  int ii;

  int width=screen->w;
  int height=screen->h;

  file=TIFFOpen(path,"w");

  if(file){
    TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
    TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
    TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
    TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(file, TIFFTAG_EXTRASAMPLES, 0);
    TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, "");

    p = (Uint8 *)screen->pixels;
    for (ii = height - 1; ii >= 0; ii--) {
      if (TIFFWriteScanline(file, p, ii, 0) < 0) {
	TIFFClose(file);
	printf("Error writing TIFF file.\n");
	exit(1);
      }
      p += 3 * width * sizeof(Uint8);
    }
    TIFFClose(file);
  }
}
#endif

/* field integration thread */
void *integration_thread(void *threadarg){
  // loop variables
  int ii;
  int jj;

  // thread slice boundaries
  int lo;
  int hi;

  // state vectors
  double *u0;
  double *u1;
  double *u2;

  // transmission vector
  double *c;

  // nonlinear term variable
  double dd;

  double dt;
  double beta;

  // thread data pointer
  struct thread_data *my_data;

  // thread enumeration variable
  int thread_id;

  // set up pointers
  my_data=(struct thread_data *) threadarg;

  u0=my_data->u0;
  u1=my_data->u1;
  u2=my_data->u2;

  c=my_data->c;

  dt=my_data->dt;
  beta=my_data->beta;

  thread_id=my_data->thread_id;

  // compute thread slice
  lo=thread_id*NUM/NUM_THREADS;
  hi=(thread_id+1)*NUM/NUM_THREADS;

  // don't integrate on lattice edges
  if(lo==0)
    lo=1;
  if(hi==NUM)
    hi=NUM-1;

  // update state
  for(jj=1; jj<NUM-1; jj++){
    for(ii=lo; ii<hi; ii++){
      // compute time difference
      u2[jj*NUM+ii]=2*u1[jj*NUM+ii]-u0[jj*NUM+ii];

      // compute spatial difference
      u2[jj*NUM+ii]+=dt*(u1[jj*NUM+ii+1]+u1[jj*NUM+ii-1]-4*u1[jj*NUM+ii]+
			 u1[(jj-1)*NUM+ii]+u1[(jj+1)*NUM+ii]);

      // compute nonlinear terms
      dd=u1[jj*NUM+ii+1]-u1[jj*NUM+ii];
      u2[jj*NUM+ii]+=beta*(dd*dd*dd);
      dd=u1[jj*NUM+ii-1]-u1[jj*NUM+ii];
      u2[jj*NUM+ii]+=beta*(dd*dd*dd);
      dd=u1[(jj+1)*NUM+ii]-u1[jj*NUM+ii];
      u2[jj*NUM+ii]+=beta*(dd*dd*dd);
      dd=u1[(jj-1)*NUM+ii]-u1[jj*NUM+ii];
      u2[jj*NUM+ii]+=beta*(dd*dd*dd);

      // transmission factor
      u2[jj*NUM+ii]*=c[jj*NUM+ii];
    }
  }

  pthread_exit(NULL);
}

/* compute time difference in seconds and nanoseconds */
void timediff(struct timespec start, struct timespec end, struct timespec *out){
  /* compute time difference */
  if(end.tv_nsec<start.tv_nsec){
    out->tv_nsec=end.tv_nsec-start.tv_nsec+1000000000;
    out->tv_sec=end.tv_sec-start.tv_sec-1;
  }
  else{
    out->tv_nsec=end.tv_nsec-start.tv_nsec;
    out->tv_sec=end.tv_sec-start.tv_sec;
  }
}

int main(int argc, char *argv[])
{
  // loop variables
  int nn;
  int ii;
  int jj;
  int tt;
  int ii2;
  int jj2;

  int steps;

  // lattice point coordinate vectors
  double *x;
  double *y;

  double xx;
  double yy;

  // initial condition variables
  double d;
  double rho=40;

  // state vectors
  double *u0;
  double *u1;
  double *u2;
  double *u3;

  double *c;

  // time step variable
  double dt=1.0/(DT*DT);

  // nonlinear factor variable
  double beta=BETA;

  // SDL variables
  SDL_Surface *screen;
  SDL_Event event;
  Uint8 *pixels;
  int done;

  // posix thread variables
  int thread_rc;
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  void *thread_status;
  
  // time measurement variables
  struct timespec int_time;
  struct timespec time1, time2;

#if TIFF_ENABLE
  // tiff writer variables
  char filename[128];
  int tiff_frame;
#endif

#if WAV_ENABLE
  // wav write variables 
  const char* outfilename="fpu0001.wav";
  SF_INFO sfinfo;
  SNDFILE * outfile;
  double buffer[FRAMESKIP];

  // set up sndfile wav write
  sfinfo.channels = 1;
  sfinfo.samplerate = 44100;
  sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    
  // open to file
  outfile=(SNDFILE *)sf_open(outfilename, SFM_WRITE, &sfinfo);
  if (!outfile){
    printf("Error: Can't open output file.\n");
    exit(1);
  }
#endif

  // open a SDL window
  SDL_Init(SDL_INIT_VIDEO);
  screen = SDL_SetVideoMode(SCALE*NUM, SCALE*NUM, 24, SDL_SWSURFACE);
  SDL_WM_SetCaption("FPU", "FPU");

  // allocate lattice coordinate vectors
  x=(double *)malloc(NUM*NUM*sizeof(double));
  if(!x){
    printf("Out of memory: x not allocated.\n");
    exit(1);
  }

  y=(double *)malloc(NUM*NUM*sizeof(double));
  if(!x){
    printf("Out of memory: y not allocated.\n");
    exit(1);
  }

  // allocate state vectors
  u0=(double *)malloc(NUM*NUM*sizeof(double));
  if(!u0){
    printf("Out of memory: u0 not allocated.\n");
    exit(1);
  }

  u1=(double *)malloc(NUM*NUM*sizeof(double));
  if(!u1){
    printf("Out of memory: u1 not allocated.\n");
    exit(1);
  }

  u2=(double *)malloc(NUM*NUM*sizeof(double));
  if(!u2){
    printf("Out of memory: u2 not allocated.\n");
    exit(1);
  }

  c=(double *)malloc(NUM*NUM*sizeof(double));
  if(!u2){
    printf("Out of memory: c not allocated.\n");
    exit(1);
  }

  // build up lattice coordinate vectors
  for(jj=0; jj<NUM; jj++){
    for(ii=0; ii<NUM; ii++){
      x[jj*NUM+ii]=XMIN+(XMAX-XMIN)*(double)((double)ii/(NUM-1));
      y[jj*NUM+ii]=XMIN+(XMAX-XMIN)*(double)((double)jj/(NUM-1));
    }
  }

  // initialize transmisson vector
  for(jj=0; jj<NUM; jj++){
    for(ii=0; ii<NUM; ii++){
      d=sqrt(x[jj*NUM+ii]*x[jj*NUM+ii]+y[jj*NUM+ii]*y[jj*NUM+ii]);
      if(d<1.0)
	c[jj*NUM+ii]=1.0;
      else
	c[jj*NUM+ii]=0.0;	
    }
  }

  // build up initial condition
  for(jj=0; jj<NUM; jj++){
    for(ii=0; ii<NUM; ii++){
      xx=x[jj*NUM+ii];
      yy=y[jj*NUM+ii];
      u0[jj*NUM+ii]=exp(-rho*(xx*xx+yy*yy));
      u1[jj*NUM+ii]=u0[jj*NUM+ii];
    }
  }

  // start timer
  clock_gettime(CLOCK_MONOTONIC, &time1);

  // integrate time steps
  done=0;
  tt=0;
  steps=0;

#if TIFF_FRAME
  tiff_frame=0;
#endif

  while(!done){
    // check for SDL event
    if (SDL_PollEvent(&event)) {
      switch (event.type) {
        // close button clicked
      case SDL_QUIT:
	done=1;
	break;
      }
    }

    // create integration threads
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(nn=0;nn<NUM_THREADS;nn++){
      thread_data_array[nn].thread_id=nn;
      thread_data_array[nn].u0=u0;
      thread_data_array[nn].u1=u1;
      thread_data_array[nn].u2=u2;
      thread_data_array[nn].c=c;
      thread_data_array[nn].dt=dt;
      thread_data_array[nn].beta=beta;
      thread_rc=pthread_create(&threads[nn], &attr, integration_thread, (void *) &thread_data_array[nn]);
    
      if(thread_rc){
	printf("ERROR: pthread_create() returned %d.\n", thread_rc);
	exit(-1);
      }
    }
  
    /* join threads */
    pthread_attr_destroy(&attr);
    for(nn=0;nn<NUM_THREADS;nn++){
      thread_rc=pthread_join(threads[nn], &thread_status);
      
      if(thread_rc){
	printf("ERROR: pthread_join() returned %d.\n", thread_rc);
	exit(-1);
      }
    }

    // update time slices
    u3=u0;
    u0=u1;
    u1=u2;
    u2=u3;

#if WAV_ENABLE
    // update wave buffer
    buffer[steps]=u0[NUM/2*NUM+NUM/2];
#endif

    // draw state
    if(steps==FRAMESKIP-1){
      // timer stop
      clock_gettime(CLOCK_MONOTONIC, &time2);

      pixels=(Uint8 *)screen->pixels;
      SDL_LockSurface(screen);
      for(jj=0; jj<NUM; jj++){
	for(ii=0; ii<NUM; ii++){
	  // get field value
	  d=256.0*(4.0*u2[jj*NUM+ii]+0.5)*c[jj*NUM+ii];

	  // clip value
	  if(d>255.0)
	    d=255.0;
	  if(d<0.0)
	    d=0.0;
	  
	  // update framebuffer
	  for(jj2=0; jj2<SCALE; jj2++){
	    for(ii2=0; ii2<SCALE; ii2++){
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+0]=(Uint8)d;
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+1]=(Uint8)d;
	      pixels[3*(SCALE*jj+jj2)*screen->w+3*(SCALE*ii+ii2)+2]=(Uint8)d;
	    }
	  }
	}
      }
      SDL_UnlockSurface(screen);
      SDL_Flip(screen);

      timediff(time1, time2, &int_time);

      // print state
      printf("tt: %d beta: %f steps/s: %d\n", tt, beta, (int)(steps/((double)(int_time.tv_sec)+(double)(int_time.tv_nsec)*1.0E-9)));

#if WAV_ENABLE      
      // write wave to disk
      sf_write_double(outfile, &buffer[0], FRAMESKIP) ;
#endif

#if TIFF_ENABLE
      sprintf(filename, "/home/janne808/testrun/%08d.tif", tiff_frame++);
      writeframe(filename, screen);
#endif

      // timer start
      clock_gettime(CLOCK_MONOTONIC, &time1);

      // reset steps variable
      steps=-1;
    }

    // update time step variable
    tt++;
    steps++;
  }

#if WAV_ENABLE
  // force write to disk and close
  sf_write_sync(outfile);
  sf_close(outfile);
#endif

  // free lattice coordinate vectors
  free(y);
  free(x);

  // free state vectors
  free(u2);
  free(u1);
  free(u0);

  // clean up SDL
  SDL_Quit();

  return 0;
}

