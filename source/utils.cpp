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

#include <cstring>
#include <fstream>

#include <stdint.h>
#include <unistd.h>

/* 32-bit */
#if UINTPTR_MAX == 0xffffffff
#ifndef ACML
#include <mkl_vml.h>
#endif /* ACML */
#define _32BIT

/* 64-bit */
#elif UINTPTR_MAX == 0xffffffffffffffff
#ifndef ACML
#include <mkl_vml.h>
#elif defined ACML4
#include <acml_mv.h>
#elif defined ACML5
#include <amdlibm.h>
#endif /* ACML */

/* wtf */
#else
#include <mkl_vml.h>
//#error Unknown integer size

#endif /* UINTPTR_MAX */

#include "potential.h"
#include "utils.h"

using namespace POTFIT_NS;

Utils::Utils(POTFIT *ptf) :
  Pointers(ptf),
  t_begin(0),
  t_end(0),
  flagfile("STOP")
{}

Utils::~Utils() {}

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

double Utils::vect_dist(const vector &a, const vector &b) {
  return sqrt(square(b.x-a.x)+square(b.y-a.y)+square(b.z-a.z));
}


void Utils::start_timer(void) {
  time(&t_begin);

  return;
}

void Utils::end_timer(void) {
  time(&t_end);

  return;
}

const int Utils::timediff(void) {
  return (int)difftime(t_end, t_begin);
}

void Utils::set_flagfile(const std::string &str) {
  flagfile = str;

  return;
}

const int Utils::check_for_flagfile(void) {
  std::ifstream infile(flagfile.c_str());
  if (infile) {
    infile.close();
    unlink(flagfile.c_str());
    return 1;
  }

  return 0;
}

void POTFIT_NS::power_1(double &result, const double &x, const double &y) {
#ifdef _32BIT
  *result = pow(*x, *y);
#else
#ifndef ACML
  vdPow(1, &x, &y, &result);
#elif defined ACML4
  *result = fastpow(*x, *y);
#elif defined ACML5
  *result = pow(*x, *y);
#endif /* ACML */
#endif /* _32BIT */

  return;
}

void POTFIT_NS::power_m(int const& dim, double *result, const double *x, const double *y) {
#ifdef _32BIT
  int   i = 0;
  for (i = 0; i < dim; i++)
    result[i] = pow(x[i], y[i]);
#else
#ifndef ACML
  vdPow(dim, x, y, result);
#elif defined ACML4
  int   i;
  for (i = 0; i < dim; i++)
    *(result + i) = fastpow(*(x + i), *(y + i));
#elif defined ACML5
  int   i;
  for (i = 0; i < dim; i++)
    *(result + i) = pow(*(x + i), *(y + i));
#endif /* ACML */
#endif /* _32BIT */


  return;
}
