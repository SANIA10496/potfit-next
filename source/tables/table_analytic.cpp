/****************************************************************
 *
 * table_analytic.cpp:
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
#include <cstring>

#include "table_analytic.h"

#include "../io.h"
#include "../memory.h"
#include "../settings.h"

#include "../functions/list_functions.h"

using namespace POTFIT_NS;

TableAnalytic::TableAnalytic(POTFIT *ptf) : Table(ptf) {
  bare = 0;

  name = NULL;
  begin = 0.0;
  end = 0.0;
  n_par = 0;
  invar = 0;

  values = NULL;
  param_name = NULL;
  val_min = NULL;
  val_max = NULL;

  function = NULL;
}

TableAnalytic::~TableAnalytic() {
}

void TableAnalytic::init(const char *fname) {
  if (strcmp(fname,"none") == 0)
    return;
#define FUNCTION_TYPE
#define FunctionType(key,Class) \
  else if (strcmp(fname,#key) == 0) \
    function = new Class();
#include "../functions/list_functions.h"
#undef FUNCTION_TYPE
  return;
}

void TableAnalytic::init_bare(const char *name_str, int num_par) {
  bare = 1;

  n_par = num_par;

  name = (char *)malloc(40 * sizeof(char));
  strcpy(name,name_str);
  param_name = (char **)malloc(n_par * sizeof(char *));
  if (NULL == param_name)
    io->error("Could not allocate memory for potential name\n");
  for (int i=0; i<n_par; i++) {
    param_name[i] = (char *)malloc(20 * sizeof(char));
    if (NULL == param_name[i])
      io->error("Could not allocate memory for parameter names\n");
    strcpy(param_name[i],"\0");
  }
  memory->create(values,n_par,"potential values");
  memory->create(val_min,n_par,"minimum potential values");
  memory->create(val_max,n_par,"maximum potential values");
  memory->create(invar_par,n_par,"invariant switch");

  return;
}

const char *TableAnalytic::get_param_name(int index) {
  if (index>=n_par)
    io->error("There are only %d parameters in this potential table (%s).\n",n_par,name);
  return param_name[index];
}

void TableAnalytic::set_value(int a, double b, double c) {
}

void TableAnalytic::set_value(int number, const char *name_str, double val, double min, double max) {
  char msg[255];

  /* check for invariance and proper value (respect boundaries) */
  /* parameter will not be optimized if min==max */
  if (min == max) {
    invar_par[number] = 1;
// what is this good for?
//        apt->invar_par[i][apt->globals]++;
  } else if (min > max) {
    double temp = min;
    min = max;
    max = temp;
  } else if ((val < min) || (val > max)) {
    /* Only print warning if we are optimizing */
    if (settings->opt) {
      if (val < min)
        val = min;
      if (val > max)
        val = max;
      sprintf(msg, "Starting value for %s #%d is ", name, number + 1);
      sprintf(msg, "%soutside of specified adjustment range.\n",msg);
      io->warning("%sResetting it to %f.\n", msg, number + 1, val);
      if (val == 0)
        io->warning("New value is 0 ! Please be careful about this.\n");
    }
  }
  strcpy(param_name[number],name);
  values[number] = val;
  val_min[number] = min;
  val_max[number] = max;

  return;
}
