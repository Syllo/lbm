/*
 * Copyright (c) 2020 Philippe Helluy <helluy@unistra.fr>
 *                    Maxime Schmitt <maxime.schmitt@manchester.ac.uk>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef D2Q9__H_
#define D2Q9__H_

#include <stdbool.h>
#include <stdio.h>

// simple class for solving the 2d navier-stokes equation
// in a rectangle by the d2q9 Lattice Boltzmann Method (LBM)
// flow around a cylinder

// relaxation parameter
#define _RELAX 1.9

#include <stdbool.h>
#include <stddef.h>

#define VLA_3D_definition(size1, size2, size3, name, ptr)                      \
  double(*restrict name)[size2][size3] = ptr;
#define VLA_3D_size(size1, size2, size3) sizeof(double[size1][size2][size3])

typedef struct d2q9 {

  // number of points in each direction
  size_t nx, ny;

  // sizes
  double Lx, Ly;

  // grid step
  double dx;

  // nine shift velocities
  int vel[9][2];

  // scaling speed parameter
  double smax;

  // one kinetic function per velocity
  void *f;
  void *fnext;

  // density rho
  // velocity vector vx, vy
  // three macro variables: rho, rho * vx and rho * vy
  void *w;

  // current and final time
  double tnow, tmax;
  // time step
  double dt;

  // cylinder radius
  double Rc;
  // cylinder distance from the left boundary
  double Dx;

} d2q9;

// returns true if a boundary condition has to be
// imposed at point x and y
bool mask(double x, double y, d2q9 *lbm);

// imposed data at the boundaries
// x=0 or y=0 -> vx= 1 and vy = 0 rho = 1
// internal disk vx=vy=0 and rho = 1
void imposed_data(double x, double y, double t, double *w, d2q9 *lbm);

// void struct
extern d2q9 d2q9_null;

// constructor
void d2q9_init(double Lx, size_t nx, size_t ny, double Rc, double Dx,
               d2q9 *lbm);

// compute w[3] from f[9]
inline void kin_to_fluid(const double *restrict f, double *restrict w,
                         d2q9 *lbm) {

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  double c = lbm->smax;

  for (size_t i = 0; i < 9; i++) {
    w[0] = w[0] + f[i];
    w[1] = w[1] + c * lbm->vel[i][0] * f[i];
    w[2] = w[2] + c * lbm->vel[i][1] * f[i];
  }
}

// compute f[9] ("maxwellian state) from w[3]
void fluid_to_kin(const double *w, double *f, d2q9 *lbm);

// perform the shift process
void d2q9_shift(d2q9 *lbm);

// perform the relax process
void d2q9_relax(d2q9 *lbm);

// apply boundary conditions
void d2q9_boundary(d2q9 *lbm);

// performs one time step of the LBM algorithm
void d2q9_step(d2q9 *lbm);
// and all the steps
void d2q9_solve(d2q9 *lbm, double tmax, bool verbose);

// display the results iv (iv=0 -> rho,
// iv=1 -> rho * u, iv=2 -> rho * v, iv=3 -> sqrt(u*u+v*v))
void d2q9_dump(FILE *out, d2q9 *lbm, int var);
// dumps the data as a 2D CSV (but space separated) matrix
void d2q9_dump_matrix(FILE *out, d2q9 *lbm, int var);

#endif // D2Q9__H_