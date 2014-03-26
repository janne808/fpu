// fpu headers

/*
 *  (C) 2014 Janne Heikkarainen <janne.heikkarainen@tut.fi>
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
