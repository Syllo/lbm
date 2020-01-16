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

#include "d2q9.h"
#include "time_measurement.h"
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool mask(double x, double y, d2q9 *lbm) {

  double xmil = lbm->Dx;

  double Lx = lbm->Lx;
  double Ly = lbm->Ly;

  double ymil = Ly / 2.;

  if (x * x < lbm->dx * lbm->dx || y * y < lbm->dx * lbm->dx ||
      (x - xmil) * (x - xmil) + (y - ymil) * (y - ymil) < lbm->Rc * lbm->Rc)
    return true;
  else
    return false;
}

void imposed_data(double x, double y, double t, double *w, d2q9 *lbm) {
  (void)t;

  double rho, u, v;
  rho = 1;
  u = 0.03;
  v = 0.00001; // to get instability

  double xmil = lbm->Dx;

  double Lx = lbm->Lx;
  double Ly = lbm->Ly;

  double ymil = Ly / 2;

  if ((x - xmil) * (x - xmil) + (y - ymil) * (y - ymil) < lbm->Rc * lbm->Rc) {
    u = 0;
  }

  w[0] = rho;
  w[1] = rho * u;
  w[2] = rho * v;
}

void d2q9_init(double Lx, size_t nx, size_t ny, double Rc, double Dx,
               d2q9 *lbm) {

  lbm->Lx = Lx;
  lbm->nx = nx;
  lbm->ny = ny;
  lbm->Rc = Rc;
  lbm->Dx = Dx;

  lbm->Ly = Lx * (double)ny / (double)nx;

  lbm->dx = Lx / (double)nx;

  lbm->vel[0][0] = 0;
  lbm->vel[0][1] = 0;

  lbm->vel[1][0] = 1;
  lbm->vel[1][1] = 0;

  lbm->vel[2][0] = 0;
  lbm->vel[2][1] = 1;

  lbm->vel[3][0] = -1;
  lbm->vel[3][1] = 0;

  lbm->vel[4][0] = 0;
  lbm->vel[4][1] = -1;

  lbm->vel[5][0] = 1;
  lbm->vel[5][1] = 1;

  lbm->vel[6][0] = -1;
  lbm->vel[6][1] = 1;

  lbm->vel[7][0] = -1;
  lbm->vel[7][1] = -1;

  lbm->vel[8][0] = 1;
  lbm->vel[8][1] = -1;

  lbm->smax = 1;

  lbm->f = calloc(1, VLA_3D_size(9, nx, ny));
  lbm->fnext = calloc(1, VLA_3D_size(9, nx, ny));
  lbm->w = calloc(1, VLA_3D_size(3, nx, ny));

  VLA_3D_definition(9, nx, ny, f_, lbm->f);
  VLA_3D_definition(3, nx, ny, w_, lbm->w);

  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      double w[3];
      double x = (double)i * lbm->dx;
      double y = (double)j * lbm->dx;
      double t = 0;
      imposed_data(x, y, t, w, lbm);
      w_[0][i][j] = w[0];
      w_[1][i][j] = w[1];
      w_[2][i][j] = w[2];

      double f[9];
      // for(int k = 0; k < 3; k++) w[k] = lbm->w[k][i + j * nx];
      fluid_to_kin(w, f, lbm);
      for (size_t k = 0; k < 9; k++)
        f_[k][i][j] = f[k];
    }
  }
}

void d2q9_shift(d2q9 *lbm) {
  VLA_3D_definition(9, lbm->nx, lbm->ny, f_, lbm->f);
  VLA_3D_definition(9, lbm->nx, lbm->ny, fnext_, lbm->fnext);
  for (size_t k = 0; k < 9; k++) {
#pragma omp for schedule(static) nowait
    for (size_t i = 0; i < lbm->nx; i++) {
      size_t i2 = (lbm->nx + i - (size_t)lbm->vel[k][0]) % lbm->nx;
      for (size_t j = 0; j < lbm->ny; j++) {
        size_t j2 = (lbm->ny + j - (size_t)lbm->vel[k][1]) % lbm->ny;
        // printf("j=%d j2=%d\n",j,j2);
        // assert(j2 >=0 && j2 < lbm->ny);
        fnext_[k][i][j] = f_[k][i2][j2];
      }
    }
  }
}

void d2q9_step(d2q9 *lbm) {
  d2q9_shift(lbm);
  d2q9_relax(lbm);
  d2q9_boundary(lbm);
}

void d2q9_relax(d2q9 *lbm) {

  VLA_3D_definition(9, lbm->nx, lbm->ny, f_, lbm->f);
  VLA_3D_definition(9, lbm->nx, lbm->ny, fnext_, lbm->fnext);
  VLA_3D_definition(3, lbm->nx, lbm->ny, w_, lbm->w);
#pragma omp for schedule(static) nowait
  for (size_t i = 0; i < lbm->nx; i++) {
    for (size_t j = 0; j < lbm->ny; j++) {
      double f[9];
      double feq[9];
      for (size_t k = 0; k < 9; k++) {
        f[k] = fnext_[k][i][j];
      }
      double w[3];
      kin_to_fluid(f, w, lbm);
      for (size_t k = 0; k < 3; k++) {
        w_[k][i][j] = w[k];
        // printf("u=%f\n",w[1]/w[0]);
      }
      fluid_to_kin(w, feq, lbm);
      for (size_t k = 0; k < 9; k++) {
        f_[k][i][j] = _RELAX * feq[k] + (1 - _RELAX) * fnext_[k][i][j];
      }
    }
  }
}

void d2q9_boundary(d2q9 *lbm) {

  VLA_3D_definition(9, lbm->nx, lbm->ny, f_, lbm->f);
#pragma omp for schedule(static) nowait
  for (size_t i = 0; i < lbm->nx; i++) {
    for (size_t j = 0; j < lbm->ny; j++) {
      double x = (double)i * lbm->dx;
      double y = (double)j * lbm->dx;
      if (mask(x, y, lbm)) {
        double wb[3];
        imposed_data(x, y, lbm->tnow, wb, lbm);
        double fb[9];
        fluid_to_kin(wb, fb, lbm);
        for (size_t k = 0; k < 9; k++) {
          f_[k][i][j] = _RELAX * fb[k] + (1 - _RELAX) * f_[k][i][j];
        }
      }
    }
  }
}

void d2q9_solve(d2q9 *lbm, double tmax, bool verbose) {
  lbm->tmax = tmax;
  double local_tnow = lbm->tnow;
  const double dt = lbm->dx / lbm->smax;

  const double num_iter_d = ceil((tmax - local_tnow) / dt);
  double print_interval_d;
  if (num_iter_d >= 10.) {
    double divide = 1.;
    do {
      print_interval_d = ceil(num_iter_d / divide);
      divide = divide + 1.;
    } while (num_iter_d / print_interval_d < 10.);
  } else {
    print_interval_d = 1.;
  }
  const size_t print_interval = (size_t)print_interval_d;
  const double percent_increment = 100. / (num_iter_d / print_interval_d);
  const size_t inter_print = print_interval - 1;
#pragma omp parallel default(none) shared(lbm, tmax, dt, inter_print, verbose, \
                                          percent_increment, print_interval)   \
    firstprivate(local_tnow)
  {
    time_measure tstart_chunk;
    get_current_time(&tstart_chunk);
    size_t iter_count = 0;
    double percentage = percent_increment;

    while (local_tnow < tmax) {
      d2q9_step(lbm);
      local_tnow += dt;
#pragma omp single nowait
      lbm->tnow += dt;

      iter_count = iter_count == inter_print ? 0 : iter_count + 1;
      if (verbose && iter_count == 0) {
#pragma omp master
        {
          time_measure tend_chunk;
          get_current_time(&tend_chunk);
          double difference = measuring_difftime(tstart_chunk, tend_chunk);
          printf("%.0f%% -- t=%e dt=%e tend=%e (%zu iter in %.3fs)\n",
                 percentage, lbm->tnow, dt, lbm->tmax, print_interval,
                 difference);
          percentage += percent_increment;
          tstart_chunk = tend_chunk;
        }
      }
#pragma omp barrier
    }
  }
}

void fluid_to_kin(const double w[restrict 3], double f[restrict 9], d2q9 *lbm) {
  static const double c2 = 1. / 3.;
  double dotvel = 0, vel2 = 0, l2 = 0, l4 = 0, c4 = 0;

  l2 = lbm->smax * lbm->smax;
  double l2_ov_c2 = l2 / c2;

  l4 = l2 * l2;
  c4 = c2 * c2;

  vel2 = (w[1] * w[1] + w[2] * w[2]) / (w[0] * w[0]);
  dotvel = sqrt(l2) * (lbm->vel[0][0] * w[1] + lbm->vel[0][1] * w[2]) / w[0];

  f[0] = (4. / 9.) * w[0] *
         (1.0 + (l2_ov_c2)*dotvel + l4 / (2. * c4) * dotvel * dotvel -
          l2 / (2. * c2) * vel2);

  for (size_t i = 1; i < 5; i++) {
    dotvel = sqrt(l2) * (lbm->vel[i][0] * w[1] + lbm->vel[i][1] * w[2]) / w[0];
    f[i] = (1. / 9.) * w[0] *
           (1.0 + (l2_ov_c2)*dotvel + l4 / (2. * c4) * dotvel * dotvel -
            l2 / (2. * c2) * vel2);
  }
  for (size_t i = 5; i < 9; i++) {
    dotvel = sqrt(l2) * (lbm->vel[i][0] * w[1] + lbm->vel[i][1] * w[2]) / w[0];
    f[i] = (1. / 36.) * w[0] *
           (1.0 + (l2_ov_c2)*dotvel + l4 / (2. * c4) * dotvel * dotvel -
            l2 / (2. * c2) * vel2);
  }
}

extern inline void kin_to_fluid(const double *restrict f, double *restrict w,
                                d2q9 *lbm);

void d2q9_dump_matrix(FILE *out, d2q9 *lbm, int var) {
  VLA_3D_definition(3, lbm->nx, lbm->ny, w_, lbm->w);

  fprintf(out, "%zu ", lbm->nx);
  for (size_t i = 0; i < lbm->nx; i++) {
    fprintf(out, "%.3f ", (double)i * lbm->dx);
  }
  fprintf(out, "\n");
  for (size_t j = 0; j < lbm->ny; j++) {
    fprintf(out, "%.3f ", (double)j * lbm->dx);
    for (size_t i = 0; i < lbm->nx; i++) {
      double w[3];
      for (size_t k = 0; k < 3; k++)
        w[k] = w_[k][i][j];
      double val;
      if (var < 3) {
        val = w[var];
      } else {
        double r = w[0];
        double u = w[1] / r;
        double v = w[2] / r;
        val = sqrt(u * u + v * v);
      }
      fprintf(out, "%.3f ", val * 100.);
    }
    fprintf(out, "\n");
  }
}

void d2q9_dump(FILE *out, d2q9 *lbm, int var) {
  VLA_3D_definition(3, lbm->nx, lbm->ny, w_, lbm->w);

  double posY = 0.;
  for (size_t j = 0; j < lbm->ny; j++, posY += lbm->dx) {
    double posX = 0.;
    for (size_t i = 0; i < lbm->nx; i++, posX += lbm->dx) {
      double val;
      if (var < 3) {
        val = w_[var][i][j];
      } else {
        double r = w_[0][i][j];
        double u = w_[1][i][j] / r;
        double v = w_[2][i][j] / r;
        val = sqrt(u * u + v * v);
      }
      fprintf(out, "%e %e %e\n", posX, posY, val);
    }
  }
}