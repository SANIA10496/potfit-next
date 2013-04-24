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
  // from table.h
  init_done = 0;

  name = NULL;

  begin = 0.0;
  end = 0.0;

  pot_number = 0;
  num_params = 0;
  num_free_params = 0;

  // from table_analytic.h
  smooth_pot = 0;

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

    pot_number = index;

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

    num_params = function->num_params();

    // add one parameter for cutoff function if _sc is found
    if (smooth_pot == 1)
      num_params++;

    memory->create(values, num_params, "potential values");
    memory->create(val_min, num_params, "minimum potential values");
    memory->create(val_max, num_params, "maximum potential values");
    memory->create(invar_par, num_params, "maximum potential values");

    for (int i=0; i<num_params; i++) {
      values[i] = 0.0;
      val_min[i] = 0.0;
      val_max[i] = 0.0;
      invar_par[i] = 0;
    }

    param_name = (char **)malloc(num_params * sizeof(char *));
    for (int i=0; i<num_params; i++) {
      param_name[i] = (char *)malloc(30 * sizeof(char));
      strcpy(param_name[i],"\0");
    }

    return;
  } else
    io->error("This potential is already initialized.\n");

  return;
}

void TableAnalytic::read_potential(FILE *infile) {
  char buffer[255], msg[255];
  int i, j, ret_val;
  fpos_t filepos;

  if (!init_done)
    io->error("Please initialize the potential before reading any potentials.\n");

  // read cutoff
  if (2 > fscanf(infile, "%s %lf", buffer, &end))
    io->error("Could not read cutoff for potential #%d in potential file.\n", pot_number);
  if (strcmp(buffer, "cutoff") != 0)
    io->error("No cutoff found for the %d. potential.\n", pot_number);
  // set very small begin, needed for EAM embedding function
  begin = .0001;

  // check for comments
  do {
    i = fgetc(infile);
  } while (i != 10);

  fgetpos(infile, &filepos);
  fgets(buffer, 255, infile);
  while (buffer[0] == '#') {
    fgetpos(infile, &filepos);
    fgets(buffer, 255, infile);
  }
  fsetpos(infile, &filepos);

  // read parameters
  num_free_params = num_params;
  for (i = 0; i < num_params; i++) {
    param_name[i] = (char *)malloc(30 * sizeof(char));
    if (NULL == param_name[i])
      io->error("Error in allocating memory for parameter name");
    strcpy(param_name[i], "\0");
    fgetpos(infile, &filepos);
    ret_val = fscanf(infile, "%s %lf %lf %lf", buffer, &values[i], &val_min[i], &val_max[i]);
    strncpy(param_name[i], buffer, 30);

    // if last char of name is "!" we have a global parameter
    if (strrchr(param_name[i], '!') != NULL) {
      if (potential->num_globals == 0)
        io->error("You need to define a global parameter before using it!");
      param_name[i][strlen(param_name[i]) - 1] = '\0';
      j = potential->global_params->get_index(param_name[i]);
      if (j<0)
        io->error("Could not find global parameter %s!", param_name[i]);
      sprintf(param_name[i], "%s!", param_name[i]);

      // register global parameter
      if (potential->invar_pot[pot_number] == 0)
        potential->global_params->add_usage(j, pot_number, i);

      // get value from global parameter table
      double temp[3];
      potential->global_params->get_value(j, temp);
      values[i] = temp[0];
      val_min[i] = temp[1];
      val_max[i] = temp[2];
      invar_par[i] = 1;
      num_free_params--;
    } else {
      // this is not a global parameter
      if (4 > ret_val) {
        if (smooth_pot && i == function->num_params()) {
          if (strcmp(param_name[i], "type") == 0 || feof(infile)) {
            io->warning("No cutoff parameter given for potential #%d: adding one parameter.", pot_number + 1);
            strcpy(param_name[i], "h");
            values[i] = 1;
            val_min[i] = 0.5;
            val_max[i] = 2;
            fsetpos(infile, &filepos);
          }
        } else {
          if (strcmp(param_name[i], "type") == 0) {
            io->error("Not enough parameters for potential #%d (%s) in potential file.", pot_number + 1, name);
          }
          io->error("Could not read parameter #%d of potential #%d in potential file.", i + 1, pot_number + 1);
        }
      }

      // check for invariance and proper value (respect boundaries)
      // parameter will not be optimized if min==max
      if (val_min[i] == val_max[i]) {
        invar_par[i] = 1;
	num_free_params--;
      } else if (val_min[i] > val_max[i]) {
        double temp = val_min[i];
        val_min[i] = val_max[i];
        val_max[i] = temp;
      } else if ((values[i] < val_min[i]) || (values[i] > val_max[i])) {
        // Only print warning if we are optimizing */
        if (settings->opt) {
          if (values[i] < val_min[i])
            values[i] = val_min[i];
          if (values[i] > val_max[i])
            values[i] = val_max[i];
          sprintf(msg, "Starting value for parameter %s #%d is ", param_name[i], pot_number + 1);
          sprintf(msg, "%soutside of specified adjustment range.\n",msg);
          io->warning("%sResetting it to %f.", msg, values[i]);
          if (values[i] == 0)
            io->warning("New value is 0 ! Please be careful about this.");
        }
      }
    }
  }

  // create indirect indexing array

}

int TableAnalytic::get_number_params(void) {
  return num_params;
}

int TableAnalytic::get_number_free_params(void) {
  return num_free_params;
}

double TableAnalytic::get_cutoff(void) {
  return end;
}

double TableAnalytic::get_rmin(void) {
  return begin;
}

void TableAnalytic::set_params(double *) {

}

void TableAnalytic::write_potential(FILE *outfile) {
  io->writef(outfile, "\n");
  io->writef(outfile, "type %s",name);
  if (smooth_pot) {
    io->writef(outfile, "_sc\n");
  } else {
    io->writef(outfile, "\n");
  }
  io->writef(outfile, "cutoff %f\n",end);
  io->writef(outfile, "# r_min 1.2.3\n");
  for (int i=0; i<num_params; i++) {
    if (param_name[i][strlen(param_name[i]) - 1] != '!') {
      io->writef(outfile, "%s %f %f %f\n",param_name[i],values[i],val_min[i],val_max[i]);
    } else {
      io->writef(outfile, "%s\n",param_name[i]);
    }
  }
}

void TableAnalytic::write_plot(FILE *outfile) {
}

void TableAnalytic::write_plotpoint(FILE *outfile) {
}
