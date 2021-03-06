/****************************************************************
 *
 * random.cpp:
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include <cstdlib>

#include "random.h"

using namespace POTFIT_NS;

Random::Random(POTFIT *ptf) :
  Pointers(ptf),
  seed(0)
{}

Random::~Random() {}

void Random::set_seed(int new_seed) {
  seed = new_seed;

  /* properly initialize random number generator */
#define R_SIZE 624
#define RAND_MAX 2147483647
  uint32_t *array;
  array = (uint32_t *) malloc(R_SIZE * sizeof(uint32_t));
  srand(seed);
  for (int i = 0; i < R_SIZE; i++)
    array[i] = rand();

  dsfmt_init_by_array(&dsfmt, array, R_SIZE);
  for (int i = 0; i < 10e5; i++)
    eqdist();
  free(array);
#undef R_SIZE
#undef RAND_MAX
}

double Random::normdist(void) {
  static int have = 0;
  static double nd2;
  double x1, x2, sqr, cnst;

  if (!(have)) {
    do {
      x1 = 2.0 * eqdist() - 1.0;
      x2 = 2.0 * eqdist() - 1.0;
      sqr = x1 * x1 + x2 * x2;
    } while (!(sqr <= 1.0 && sqr > 0));
    /* Box Muller Transformation */
    cnst = sqrt(-2.0 * log(sqr) / sqr);
    nd2 = x2 * cnst;
    have = 1;
    return x1 * cnst;
  } else {
    have = 0;
    return nd2;
  }
}
