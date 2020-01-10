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
#include <getopt.h>
#include <stdbool.h>
#include <stdlib.h>

static struct option opt_options[] = {
    {"domain-length", required_argument, 0, 'l'},
    {"grid-x", required_argument, 0, 'x'},
    {"grid-y", required_argument, 0, 'y'},
    {"output", optional_argument, 0, 'o'},
    {"tmax", required_argument, 0, 't'},
    {"cylinder-radius", required_argument, 0, 'r'},
    {"cylinder-right-shift", required_argument, 0, 's'},
    {"help", no_argument, 0, 'h'},
    {"quiet", no_argument, 0, 'q'},
    {0, 0, 0, 0}};

static const char options[] = ":l:x:y:o:t:r:s:hq";

static const char help_string[] =
    "Options:"
    "\n  -l --domain-length        : Domain length (e.g., 10.)"
    "\n  -x --grid-x               : Number of grid points"
    "\n  -y --grid-y               : Number of grid points"
    "\n  -o --output               : Select the output file name"
    "\n  -t --tmax                 : Stop the simulation when the time is "
    "reached"
    "\n  -r --cylinder-radius      : Radius of the cylinder"
    "\n  -s --cylinder-right-shift : Radius of the cylinder"
    "\n  -h --help                 : Print this help";

// domain length
#define _LX 10.

// cylinder radius
#define _RC 0.3

// cylinder distance from the left boundary
#define _DX 3.

// number of grid points in each direction
#define _NX 512
#define _NY 128
#define _TMAX 600

int main(int argc, char **argv) {

  d2q9 lbm = {0};

  double Lx = _LX;
  size_t nx = _NX;
  size_t ny = _NY;
  double Rc = _RC;
  double Dx = _DX;
  char *output_filename = NULL;
  double tmax = _TMAX;
  bool verbose = true;

  while (true) {
    int sscanf_return;
    int optchar = getopt_long(argc, argv, options, opt_options, NULL);
    if (optchar == -1)
      break;
    switch (optchar) {
    case 'q':
      verbose = false;
      break;
    case 'o':
      if (optarg != NULL && optarg[0] != '-' && optarg[0] != '\0') {
        output_filename = optarg;
      } else {
        output_filename = "gridData.dat";
        if (optarg != NULL)
          optind--;
      }
      break;
    case 't':
      sscanf_return = sscanf(optarg, "%lf", &tmax);
      if (sscanf_return == EOF || sscanf_return == 0 || tmax <= 0.) {
        fprintf(stderr,
                "Please enter a positive floating point value for the "
                "simulation time instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        tmax = _TMAX;
      }
      break;
    case 'l':
      sscanf_return = sscanf(optarg, "%lf", &Lx);
      if (sscanf_return == EOF || sscanf_return == 0 || Lx <= 0.) {
        fprintf(stderr,
                "Please enter a positive floating point value for the "
                "domain length instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        Lx = _LX;
      }
      break;
    case 'x':
      sscanf_return = sscanf(optarg, "%zu", &nx);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the "
                "X number of grid points instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        nx = _NX;
      }
      break;
    case 'y':
      sscanf_return = sscanf(optarg, "%zu", &ny);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the "
                "Y number of grid points instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        ny = _NY;
      }
      break;
    case 'r':
      sscanf_return = sscanf(optarg, "%lf", &Rc);
      if (sscanf_return == EOF || sscanf_return == 0 || Rc <= 0.) {
        fprintf(stderr,
                "Please enter a positive floating point value for the "
                "cylinder radius instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        Rc = _RC;
      }
      break;
    case 's':
      sscanf_return = sscanf(optarg, "%lf", &Dx);
      if (sscanf_return == EOF || sscanf_return == 0 || Rc <= 0.) {
        fprintf(stderr,
                "Please enter a positive floating point value for the "
                "cylinder radius instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        Dx = _DX;
      }
      break;
    case 'h':
      printf("Usage: %s <options>\n%s\n", argv[0], help_string);
      return EXIT_SUCCESS;
    case ':':
      fprintf(stderr, "Option %c requires an argument\n", optopt);
      exit(EXIT_FAILURE);
      break;
    default:
      fprintf(stderr, "Unrecognized option %c\n", optopt);
      exit(EXIT_FAILURE);
      break;
    }
  }

  d2q9_init(Lx, nx, ny, Rc, Dx, &lbm);

  printf("Domain size (%f,%f) [%zu x %zu points] cylinder (%f,%f) radius %f\n",
         lbm.Lx, lbm.Ly, lbm.nx, lbm.ny, lbm.Dx, lbm.Ly / 2., lbm.Rc);

  time_measure startTime, endTime;
  get_current_time(&startTime);
  d2q9_solve(&lbm, tmax, verbose);
  get_current_time(&endTime);
  fprintf(stdout, "Kernel time %.4fs\n",
          measuring_difftime(startTime, endTime));
  if (output_filename) {
    FILE *out = fopen(output_filename, "w");
    d2q9_dump(out, &lbm, 3);
    fclose(out);
    printf("Data written to %s\n", output_filename);
  }

  return EXIT_SUCCESS;
}