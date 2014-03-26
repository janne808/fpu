// fpu headers
//
#ifndef __FPUH__
#define __FPUH__

struct thread_data{
  int thread_id;

  double *u0;
  double *u1;
  double *u2;

  double *c;

  double dt;
  double beta;
};

#endif
