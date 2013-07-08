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
#include <iomanip>

#include "splines.h"
#include "table_analytic.h"

#include "../io.h"
#include "../output.h"
#include "../potential.h"
#include "../settings.h"
#include "../structures.h"

#include "../functions/list_functions.h"

using namespace POTFIT_NS;

TableAnalytic::TableAnalytic(POTFIT *ptf) :
  Table(ptf),
  smooth_pot(0),
  function(NULL)
{
  grad[0] = 0.0;
  grad[1] = 0.0;

  return;
}

TableAnalytic::~TableAnalytic() {
  param_name.clear();
  values.clear();
  val_min.clear();
  val_max.clear();
  stored_values.clear();
  invar_par.clear();
  idx.clear();

  return;
}

void TableAnalytic::init(const std::string &fname, const int &index) {
  std::size_t pos;

  if (!init_done) {
    init_done = 1;

    pot_number = index;

    // split name and _sc
    pos = fname.find_last_of('_');
    if (std::string::npos != pos) {
      if (fname.substr(pos+1).compare("sc") == 0) {
        name = fname.substr(0,pos);
	smooth_pot = 1;
      } else {
	name = fname;
      }
    } else {
      name = fname;
    }

    if (name.compare("none") == 0)
      return;
#define FUNCTION_TYPE
#define FunctionType(key,Class) \
  else if (name.compare(#key) == 0) \
    function = new Class();
#include "../functions/list_functions.h"
#undef FUNCTION_TYPE
    if (!function) {
      io->error << "Could not create an analytic potential of type \"" << fname << "\"." << std::endl;
      io->pexit(EXIT_FAILURE);
    }

    num_params = function->num_params();

    // add one parameter for cutoff function if _sc is found
    if (smooth_pot == 1)
      num_params++;

    values.resize(num_params,0.0);
    val_min.resize(num_params,0.0);
    val_max.resize(num_params,0.0);
    stored_values.resize(num_params,0.0);
    invar_par.resize(num_params,0.0);
    param_name.resize(num_params);

    return;
  } else {
    io->error << "This potential is already initialized." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  return;
}

void TableAnalytic::read_potential(FILE *infile) {
  char buffer[255], msg[255];
  int i, j, ret_val;
  fpos_t filepos;

  if (!init_done) {
    io->error << "Please initialize the potential before reading any potentials." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  // read cutoff
  if (2 > fscanf(infile, "%s %lf", buffer, &end)) {
    io->error << "Could not read cutoff for potential #" << pot_number << " in potential file." << std::endl;
    io->pexit(EXIT_FAILURE);
  }
  if (strcmp(buffer, "cutoff") != 0) {
    io->error << "No cutoff found for the " << pot_number << ". potential." << std::endl;
    io->pexit(EXIT_FAILURE);
  }
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
    fgetpos(infile, &filepos);
    ret_val = fscanf(infile, "%s %lf %lf %lf", buffer, &values[i], &val_min[i], &val_max[i]);
    param_name[i] = buffer;

    // if last char of name is "!" we have a global parameter
    if (param_name[i][param_name[i].length()] == '!') {
      if (potential->num_globals == 0) {
        io->error  << "You need to define a global parameter before using it!" << std::endl;
        io->pexit(EXIT_FAILURE);
      }
      param_name[i].erase(param_name[i].end() - 1);
      j = potential->global_params->get_index(param_name[i]);
      if (j<0) {
        io->error << "Could not find global parameter " << param_name[i] << "!" << std::endl;
        io->pexit(EXIT_FAILURE);
      }
      param_name[i] += "!";

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
          if (param_name[i].compare("type") == 0 || feof(infile)) {
            io->warning << "No cutoff parameter given for potential #" << pot_number + 1 << ": adding one parameter." << std::endl;
            param_name[i] = "h";
            values[i] = 1;
            val_min[i] = 0.5;
            val_max[i] = 2;
            fsetpos(infile, &filepos);
          }
        } else {
          if (param_name[i].compare("type") == 0) {
            io->error << "Not enough parameters for potential #" << pot_number + 1 << " (" << name << ") in potential file." << std::endl;
            io->pexit(EXIT_FAILURE);
          }
          io->error << "Could not read parameter #" << i + 1 << " of potential #" << pot_number + 1 << " in potential file." << std::endl;
          io->pexit(EXIT_FAILURE);
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
        if (settings->get_opt()) {
          if (values[i] < val_min[i])
            values[i] = val_min[i];
          if (values[i] > val_max[i])
            values[i] = val_max[i];
          io->warning << "Starting value for parameter " << param_name[i] << " of potential #" << pot_number + 1 << " is " << std::endl;
          io->warning << "outside of specified adjustment range." << std::endl;
          io->warning << "Resetting it to " << values[i] << "." << std::endl;
          if (values[i] == 0)
            io->warning << "New value is 0 ! Please be careful about this." << std::endl;
        }
      }
    }
  }

  idx.resize(num_free_params,-1);

  init_calc_table();

  return;
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

double TableAnalytic::get_val_min(const int &n) {
  return val_min[n];
}

double TableAnalytic::get_val_max(const int &n) {
  return val_max[n];
}

double TableAnalytic::get_plotmin(void) {
  return output->get_plotmin();
}

double TableAnalytic::get_value(const double &r) {
  double temp = 0.0;

  function->calc(r, values, &temp);
  if (1 == smooth_pot)
    temp *= cutoff(r, end, values[num_params - 1]);

  return temp;
}

void TableAnalytic::init_calc_table(void) {
  double fval, h;

  len = POT_STEPS;
  step = (end - begin) / (len - 1);
  invstep = 1. / step;
  grad[0] = 10e30;
  grad[1] = 0.0;
  xcoord = potential->xcoord + pot_number * POT_STEPS;
  table = potential->table + pot_number * POT_STEPS;
  d2tab = potential->d2tab + pot_number * POT_STEPS;

  h = values[num_params - 1];
  for (int i=0; i<len; i++) {
    xcoord[i] = begin + i * step;
    function->calc(xcoord[i], values , &fval);
    table[i] = smooth_pot ? fval * cutoff(xcoord[i], end, h) : fval;
  }

  for (int i=0; i<num_params; i++)
    stored_values[i] = values[i];

  splines.spline_ed(step, table, len, 10e30, 0, d2tab);

  return;
}

void TableAnalytic::set_param(const int &i, const double& val) {
  values[i] = val;

  return;
}

void TableAnalytic::update_potential(const int &update) {
  update_values();
  update_calc_table(update);

  return;
}

void TableAnalytic::update_values(void) {
  for (int i=0;i<num_free_params;i++) {
    values[idx[i]] = potential->opt->val_p[opt_pot_start + i];
  }

  return;
}

void TableAnalytic::update_calc_table(const int &update) {
  double f, h;
  int change = 0;

  for (int i=0; i<num_params; i++)
    if (stored_values[i] != values[i])
      change = 1;

  if (0 == change && 0 == update)
    return;

  if (1 == update)
    for (int i=0; i<len; i++)
      xcoord[i] = begin + i * step;

  h = values[num_params - 1];
  for (int i=0; i<len; i++) {
    function->calc(xcoord[i], values, &f);
    table[i] = smooth_pot ? f * cutoff(xcoord[i], end, h) : f;
  }

  splines.spline_ed(step, table, len, 10e30, 0, d2tab);

  return;
}

void TableAnalytic::update_slots(void) {
  double r = 0.0, rr = 0.0;
  Config *conf = NULL;
  Atom *atom = NULL;
  Neighbor *neigh = NULL;
  Table *pot = NULL;

  for (int i = 0; i < structures->get_num_total_configs(); i++) {
    conf = structures->config[i];
    for (int j = 0; j < conf->num_atoms; j++) {
      atom = &conf->atoms[j];
      for (int k = 0; k < atom->num_neighbors; k++) {
        neigh = &atom->neighs[k];
        r = neigh->r;

        // update slots for pair potential part, slot 0
        pot = potential->pots[neigh->col[0]];
        if (r < pot->end) {
          rr = r - pot->begin;
          neigh->slot[0] = (int)(rr * pot->invstep);
          neigh->step[0] = pot->step;
          neigh->shift[0] = (rr - neigh->slot[0] * pot->step) * pot->invstep;
        }
      }
    }
  }

  return;
}

void TableAnalytic::write_potential(std::ofstream &outfile) {
  outfile << std::endl;
  outfile << "type " << name;
  if (smooth_pot) {
    outfile << "_sc" << std::endl;
  } else {
    outfile << std::endl;
  }
  outfile << "cutoff " << end << std::endl;
  outfile << "# r_min " << get_rmin() << std::endl;
  for (int i=0; i<num_params; i++) {
    if (param_name[i][param_name[i].length()-1] != '!') {
      outfile << param_name[i] << "\t\t";
      outfile << std::fixed << std::setw(12) << std::setprecision(5) << values[i];
      outfile << "\t" << std::setw(12) << std::setprecision(5) << val_min[i];
      outfile << "\t" << std::setw(12) << std::setprecision(5) << val_max[i] << std::endl;
    } else {
      outfile << param_name[i] << std::endl;
    }
  }

  return;
}

void TableAnalytic::write_plot(FILE *outfile) {
  return;
}

void TableAnalytic::write_plotpoint(FILE *outfile) {
  return;
}

double TableAnalytic::cutoff(const double& x, const double& rcut, const double& h) {
  if ((x-rcut)>0)
    return 0.0;

  fval = (x-rcut) / h;
  fval *= fval;
  fval *= fval;

  return fval / (1. + fval);
}
