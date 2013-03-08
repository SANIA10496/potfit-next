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
#include "../potential.h"
#include "../settings.h"

#include "../functions/list_functions.h"

using namespace POTFIT_NS;

TableAnalytic::TableAnalytic(POTFIT *ptf) : Table(ptf) {
  bare = 0;

  name = NULL;
  begin = 0.0;
  end = 0.0;
  n_par = 0;
  n_invar = 0;
  init_done = 0;
  smooth_pot = 0;
  pot_number = 0;

  values = NULL;
  param_name = NULL;
  val_min = NULL;
  val_max = NULL;
  invar_par = NULL;

  function = NULL;
}

TableAnalytic::~TableAnalytic() {
}

void TableAnalytic::init(const char *fname, int index) {
  char *token;
  char buffer[255];

  if (!init_done) {
    init_done = 1;

    pot_number = index + 1;

    // split name and _sc
    token = strrchr((char *)fname, '_');
    if (token != NULL && strcmp(token + 1, "sc") == 0) {
      strncpy(buffer, fname, strlen(fname) - 3);
      buffer[strlen(fname) - 3] = '\0';
      strcpy((char *)fname, buffer);
      smooth_pot = 1;
    }

    name = (char *)malloc(40 * sizeof(char));
    strcpy(name, fname);

    if (strcmp(fname,"none") == 0)
      return;
#define FUNCTION_TYPE
#define FunctionType(key,Class) \
  else if (strcmp(fname,#key) == 0) \
    function = new Class();
#include "../functions/list_functions.h"
#undef FUNCTION_TYPE
    if (!function)
      io->error("Could not create an analytic potential of type \"%s\".\n",fname);

    n_par = function->num_params();

    // add one parameter for cutoff function if _sc is found
    if (smooth_pot == 1)
      n_par++;
    potential->total_par += n_par;

    memory->create(values,n_par,"potential values");
    memory->create(val_min,n_par,"minimum potential values");
    memory->create(val_max,n_par,"maximum potential values");
    memory->create(invar_par,n_par,"maximum potential values");

    for (int i=0; i<n_par; i++) {
      values[i] = 0.0;
      val_min[i] = 0.0;
      val_max[i] = 0.0;
      invar_par[i] = 0;
    }

    param_name = (char **)malloc(n_par * sizeof(char *));
    for (int i=0; i<n_par; i++) {
      param_name[i] = (char *)malloc(30 * sizeof(char));
      strcpy(param_name[i],"\0");
    }

    return;
  } else
    io->error("This potential is already initialized.\n");
}

void TableAnalytic::init_bare(const char *name_str, int num_par) {
  if (!init_done) {
    init_done = 1;
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

    param_name = (char **)malloc(n_par * sizeof(char *));
    for (int i=0; i<n_par; i++) {
      param_name[i] = (char *)malloc(30 * sizeof(char));
      strcpy(param_name[i],"\0");
    }

    return;
  } else
    io->error("This potential is already initialized.\n");
}

void TableAnalytic::read_potential(FILE *infile) {
  char buffer[255];
  int j, ret_val;
  fpos_t filepos;

  if (bare)
    io->error("A bare potential cannot be read!\n");

  if (!init_done)
    io->error("Please initialize the potential before reading any potentials.\n");

  // read cutoff */
  if (2 > fscanf(infile, "%s %lf", buffer, &end))
    io->error("Could not read cutoff for potential #%d in potential file.\n", pot_number);
  if (strcmp(buffer, "cutoff") != 0)
    io->error("No cutoff found for the %d. potential.\n", pot_number);
  // set very small begin, needed for EAM embedding function
  begin = .0001;

  // check for comments
  do {
    j = fgetc(infile);
  } while (j != 10);

  fgetpos(infile, &filepos);
  fgets(buffer, 255, infile);
  while (buffer[0] == '#') {
    fgetpos(infile, &filepos);
    fgets(buffer, 255, infile);
  }
  fsetpos(infile, &filepos);

  // read parameters */
  n_invar = 0;
  for (j = 0; j < n_par; j++) {
    param_name[j] = (char *)malloc(30 * sizeof(char));
    if (NULL == param_name[j])
      io->error("Error in allocating memory for parameter name");
    strcpy(param_name[j], "\0");
    fgetpos(infile, &filepos);
    ret_val = fscanf(infile, "%s %lf %lf %lf", buffer, &values[j], &val_min[j], &val_max[j]);
    strncpy(param_name[j], buffer, 30);

    /* if last char of name is "!" we have a global parameter
    if (strrchr(param_name[j], '!') != NULL) {
      param_name[j][strlen(param_name[j]) - 1] = '\0';
      l = -1;
      for (k = 0; k < apt->globals; k++) {
        if (strcmp(apt->param_name[i][j], apt->param_name[global_pot][k]) == 0)
          l = k;
      }
      if (l == -1)
        error(1, "Could not find global parameter %s!\n", apt->param_name[i][j]);
      sprintf(apt->param_name[i][j], "%s!", apt->param_name[i][j]);

      /* write index array for global parameters
      if (++apt->n_glob[l] > 1) {
        apt->global_idx[l] = (int **)realloc(apt->global_idx[l], apt->n_glob[l] * sizeof(int *));
      } else {
        apt->global_idx[l] = (int **)malloc(1 * sizeof(int *));
      }
      apt->global_idx[l][apt->n_glob[l] - 1] = (int *)malloc(2 * sizeof(int));
      apt->global_idx[l][apt->n_glob[l] - 1][0] = i;
      apt->global_idx[l][apt->n_glob[l] - 1][1] = j;

      apt->values[i][j] = apt->values[global_pot][l];
      apt->pmin[i][j] = apt->pmin[global_pot][l];
      apt->pmax[i][j] = apt->pmax[global_pot][l];
      apt->invar_par[i][j] = 1;
      apt->invar_par[i][apt->n_par[i]]++;
    } else {
      /* this is no global parameter
      if (4 > ret_val) {
        if (smooth_pot[i] && j == apot_parameters(apt->names[i])) {
          if (strcmp(apt->param_name[i][j], "type") == 0 || feof(infile)) {
            warning(1, "No cutoff parameter given for potential #%d: adding one parameter.", i);
            strcpy(apt->param_name[i][j], "h");
            apt->values[i][j] = 1;
            apt->pmin[i][j] = 0.5;
            apt->pmax[i][j] = 2;
            fsetpos(infile, &filepos);
          }
        } else {
          if (strcmp(apt->param_name[i][j], "type") == 0) {
            error(0, "Not enough parameters for potential #%d (%s) in file %s!\n", i + 1, apt->names[i],
                  filename);
            error(1, "You specified %d parameters, but needed are %d.\n", j, apt->n_par[i]);
          }
          error(1, "Could not read parameter #%d of potential #%d in file %s", j + 1, i + 1, filename);
        }
      }
    }
    */
  }
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
// what is this good for? global counter for invar_pars per potential^^
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
