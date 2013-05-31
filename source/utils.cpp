/****************************************************************
 *
 * utils.cpp:
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

#include <cmath>
#include <cstring>

#include "potential.h"
#include "templates.h"
#include "utils.h"

using namespace POTFIT_NS;

Utils::Utils(POTFIT *ptf) : Pointers(ptf) {
  return;
}

Utils::~Utils() {
  return;
}

char *Utils::tolowercase(char *str) {
  int differ = 'A'-'a';
  char ch;
  int ii = strlen(str);

  for (int i=0; i <ii; i++)
  {
    strncpy(&ch,str+i,1);
    if (ch>='A' && ch<='Z')
    {
      ch = ch-differ;
      memcpy(str+i,&ch,1);
    }
  }
  return str;
}

double Utils::vect_dist(vector a, vector b) {
  return sqrt(square(b.x-a.x)+square(b.y-a.y)+square(b.z-a.z));
}

void Utils::quicksort(double *x, int low, int high, double **p) {
  int newIndex;

  if (low<high) {
    int index = (low + high) / 2.;
    newIndex = partition(x, low, high, index, p);
    quicksort(x, low, newIndex - 1, p);
    quicksort(x, newIndex + 1, high, p);
  }
}

int Utils::partition(double *x, int low, int high, int index, double **p) {
  int   store;
  double ind_val = x[index], temp;

  SWAP(x[index], x[high], temp);
  swap_population(p[index], p[high]);

  store = low;

  for (int i = low; i < high; i++)
    if (x[i] <= ind_val) {
      SWAP(x[i], x[store], temp);
      swap_population(p[i], p[store]);
      store++;
    }
  SWAP(x[store], x[high], temp);
  swap_population(p[store], p[high]);

  return store;
}

void Utils::swap_population(double *a, double *b) {
  double temp;

  for (int i = 0; i < potential->num_free_params + 2; i++) {
    SWAP(a[i], b[i], temp);
  }

  return;
}

void Utils::start_timer(void) {
  time(&t_begin);

  return;
}

void Utils::end_timer(void) {
  time(&t_end);

  return;
}

int Utils::timediff(void) {
  return (int)difftime(t_end, t_begin);
}

void Utils::set_flagfile(std::string str) {
  flagfile = str;

  return;
}

int Utils::check_for_flagfile(void) {
  return 0;
}
