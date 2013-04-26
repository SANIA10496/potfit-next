/****************************************************************
 *
 * io.cpp: misc io routines
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

#include <mpi.h>
#include <cstdlib>

#include "force.h"
#include "input.h"
#include "interaction.h"
#include "io.h"
#include "potential.h"
#include "output.h"
#include "random.h"
#include "settings.h"
#include "structures.h"
#include "types.h"
#include "utils.h"

using namespace POTFIT_NS;

Input::Input(POTFIT *ptf, int argc, char **argv) : Pointers(ptf) {
  strcpy(param_file,"\0");
  strcpy(config_file,"\0");
  strcpy(maxchfile,"\0");
  strcpy(startpot,"\0");

  curline = 0;
  enable_maxch_file = 0;

  settings->myid = MPI::COMM_WORLD.Get_rank();
  settings->num_cpus = MPI::COMM_WORLD.Get_size();

  if (settings->myid == 0) {
    io->screen = 1;
  }

  if (argc < 2)
    io->print_help();
  int version = 0;
  for (int i=1; i<argc; i++) {
    if (strcmp(argv[i],"-h") == 0) {
      io->print_help();
      continue;
    } else if (strcmp(argv[i],"-v") == 0) {
      version = 1;
      continue;
    } else if (strncmp(argv[i],"-",1) == 0) {
      io->error("Argument %s unknown.",argv[i]);
    }
  }
  if (version == 1)
    io->print_version();

  strcpy(param_file,argv[1]);

  return;
}

Input::~Input() {
  return;
}

void Input::read_parameter_file() {
  char  buffer[1024];
  char *token, *res;
  FILE *params;

  int seed = 0;
  curline = 0;

  io->write("Reading parameter file \"%s\" ... ",param_file);
  params = fopen(param_file, "r");
  if (NULL == params)
    io->error("Could not open parameter file \"%s\".\n",param_file);

  do {
    res = fgets(buffer, 1024, params);
    if (NULL == res)
      break;			/* probably EOF reached */
    curline++;
    token = strtok(buffer, " \t\r\n");
    if (NULL == token)
      continue;			/* skip blank lines */
    if (token[0] == '#')
      continue;			/* skip comments */

    /* number of atom types */
    if (strcasecmp(token, "ntypes") == 0) {
      get_param("ntypes", &structures->ntypes, PARAM_INT, 1, 1);
    }
    /* file with start potential */
    else if (strcasecmp(token, "startpot") == 0) {
      get_param("startpot", input->startpot, PARAM_STR, 1, 255);
    }
    /* file for end potential */
    else if (strcasecmp(token, "endpot") == 0) {
      get_param("endpot", output->endpot, PARAM_STR, 1, 255);
    }
    /* prefix for all output files */
    else if (strcasecmp(token, "output_prefix") == 0) {
      get_param("output_prefix", output->output_prefix, PARAM_STR, 1, 255);
      if (strcmp(output->output_prefix, "") != 0)
        output->enable_output_files = 1;
    }
    /* file for IMD potential */
    else if (strcasecmp(token, "imdpot") == 0) {
      get_param("imdpot", output->imdpot, PARAM_STR, 1, 255);
      output->enable_imd_pot = 1;
    }
    /* file for plotting */
    else if (strcasecmp(token, "plotfile") == 0) {
      get_param("plotfile", output->plotfile, PARAM_STR, 1, 255);
      output->enable_plot_file = 1;
    }
    /* file for maximal change */
    else if (strcasecmp(token, "maxchfile") == 0) {
      get_param("maxchfile", input->maxchfile, PARAM_STR, 1, 255);
      input->enable_maxch_file = 1;
    }
    /* file for pair distribution */
    else if (strcasecmp(token, "distfile") == 0) {
      get_param("distfile", output->distfile, PARAM_STR, 1, 255);
      output->enable_distfile = 1;
    }
    /* number of steps in IMD potential */
    else if (strcasecmp(token, "imdpotsteps") == 0) {
      get_param("imdpotsteps", &output->imdpotsteps, PARAM_INT, 1, 1);
    }
    /* minimum for plotfile */
    else if (strcasecmp(token, "plotmin") == 0) {
      get_param("plotmin", &output->plotmin, PARAM_DOUBLE, 1, 1);
    }
    /* exclude chemical potential from energy calculations */
    else if (strcasecmp(token, "enable_cp") == 0) {
      get_param("enable_cp", &potential->enable_cp, PARAM_INT, 1, 1);
    }
    /* how far should the imd pot be extended */
    else if (strcasecmp(token, "extend") == 0) {
      get_param("extend", &settings->extend, PARAM_DOUBLE, 1, 1);
    }
    /* file with atom configuration */
    else if (strcasecmp(token, "config") == 0) {
      get_param("config", input->config_file, PARAM_STR, 1, 255);
    }
    /* Optimization flag */
    else if (strcasecmp(token, "opt") == 0) {
      get_param("opt", &settings->opt, PARAM_INT, 1, 1);
    }
    /* break flagfile */
    else if (strcasecmp(token, "flagfile") == 0) {
      get_param("flagfile", utils->flagfile, PARAM_STR, 1, 255);
    }
    /* write radial pair distribution ? */
    else if (strcasecmp(token, "write_pair") == 0) {
      get_param("write_pair", &output->enable_pair_dist, PARAM_INT, 1, 1);
    }
    /* plotpoint file */
    else if (strcasecmp(token, "plotpointfile") == 0) {
      get_param("plotpointfile", output->plotpointfile, PARAM_STR, 1, 255);
    }
    /* temporary potential file */
    else if (strcasecmp(token, "tempfile") == 0) {
      get_param("tempfile", output->tempfile, PARAM_STR, 1, 255);
    }
    /* seed for RNG */
    else if (strcasecmp(token, "seed") == 0) {
      get_param("seed", &seed, PARAM_INT, 1, 1);
      random->set_seed(seed);
    }
    /* stopping criterion for differential evolution */
    else if (strcasecmp(token, "evo_threshold") == 0) {
      get_param("evo_threshold", &settings->evo_threshold, PARAM_DOUBLE, 1, 1);
    }
    /* starting temperature for annealing */
    else if (strcasecmp(token, "anneal_temp") == 0) {
      get_param("anneal_temp", &settings->anneal_temp, PARAM_STR, 1, 20);
    }
    /* Scaling Constant for APOT Punishment */
    else if (strcasecmp(token, "apot_punish") == 0) {
      get_param("apot_punish", &settings->apot_punish_value, PARAM_DOUBLE, 1, 1);
    }
    /* Energy Weight */
    else if (strcasecmp(token, "eng_weight") == 0) {
      get_param("eng_weight", &settings->eweight, PARAM_DOUBLE, 1, 1);
    }
    /* error margin delta epsilon */
    else if (strcasecmp(token, "d_eps") == 0) {
      get_param("d_eps", &settings->d_eps, PARAM_DOUBLE, 1, 1);
    }
    /* Energy Weight */
    else if (strcasecmp(token, "stress_weight") == 0) {
      get_param("stress_weight", &settings->sweight, PARAM_DOUBLE, 1, 1);
    }
    /* write final potential in lammps format */
    else if (strcasecmp(token, "write_lammps") == 0) {
      get_param("write_lammps", &output->enable_lammps_pot, PARAM_INT, 1, 1);
    }
    /* cutoff-radius for long-range interactions */
    else if (strcasecmp(token, "dp_cut") == 0) {
      get_param("dp_cut", &settings->dp_cut, PARAM_DOUBLE, 1, 1);
    }
    /* dipole iteration precision */
    else if (strcasecmp(token, "dp_tol") == 0) {
      get_param("dp_tol", &settings->dp_tol, PARAM_DOUBLE, 1, 1);
    }
    /* mixing parameter for damping dipole iteration loop */
    else if (strcasecmp(token, "dp_mix") == 0) {
      get_param("dp_mix", &settings->dp_mix, PARAM_DOUBLE, 1, 1);
    }
    /* log file */
    else if (strcasecmp(token, "write_log") == 0) {
      get_param("write_log", &io->write_logfile, PARAM_INT, 1, 1);
    }
    // interaction type
    else if (strcasecmp(token, "interaction") == 0) {
      get_param("interaction", interaction->type, PARAM_STR, 1, 255);
    }
    /* unknown tag */
    else {
      fprintf(stderr, "Unknown tag <%s> ignored!\n", token);
      fflush(stderr);
    }
  } while (!feof(params));

  io->set_logfile("potfit.log");
  io->write("done\n");
  fclose(params);
  check_params();
  interaction->init();

  return;
}

int Input::get_param(const char *param_name, void *param, param_t ptype, int pnum_min, int pnum_max)
{
  char *str;
  int   i;
  int   numread;

  numread = 0;
  switch (ptype) {
  case PARAM_STR:
    str = strtok(NULL, " \t\r\n");
    if (str == NULL)
      io->error("Parameter for %s missing in line %d\nstring expected!", param_name, curline);
    else
      strncpy((char *)param, str, pnum_max);
    numread++;
    break;
  case PARAM_INT:
    for (i = 0; i < pnum_min; i++) {
      str = strtok(NULL, " \t\r\n");
      if (str == NULL)
        io->error("Parameter for %s missing in line %d!\nInteger vector of length %u expected!",
                  param_name, curline, (unsigned)pnum_min);
      else
        ((int *)param)[i] = atoi(str);
      numread++;
    }
    for (i = pnum_min; i < pnum_max; i++) {
      if ((str = strtok(NULL, " \t\r\n")) != NULL) {
        ((int *)param)[i] = atoi(str);
        numread++;
      } else
        break;
    }
    break;
  case PARAM_DOUBLE:
    for (i = 0; i < pnum_min; i++) {
      str = strtok(NULL, " \t\r\n");
      if (str == NULL)
        io->error("Parameter for %s missing in line %d!\nDouble vector of length %u expected!",
                  param_name, curline, (unsigned)pnum_min);
      else
        ((double *)param)[i] = atof(str);
      numread++;
    }
    for (i = pnum_min; i < pnum_max; i++) {
      if ((str = strtok(NULL, " \t\r\n")) != NULL) {
        ((double *)param)[i] = atof(str);
        numread++;
      } else
        break;
    }
    break;
  }
  return numread;
}

/****************************************************************
 *
 *  check all parameters for reasonable values and completeness
 *
 ****************************************************************/

void Input::check_params()
{
  if (structures->ntypes <= 0)
    io->error("Missing parameter or invalid value in %s : ntypes is \"%d\"",
              param_file, structures->ntypes);

  if (strcmp(startpot, "\0") == 0)
    io->error("Missing parameter or invalid value in %s : startpot is \"%s\"",
              param_file, startpot);

  if (strcmp(output->endpot, "\0") == 0) {
    io->warning("endpot is missing in %s, setting it to %s_end", param_file, startpot);
    sprintf(output->endpot, "%s_end", startpot);
  }

  if (strcmp(interaction->type, "\0") == 0)
    io->error("Error in %s - 'interaction' keyword is missing!", param_file);

  if (strcmp(config_file, "\0") == 0)
    io->error("Missing parameter or invalid value in %s: config is \"%s\"",
              param_file, config_file);

  if (strcmp(output->tempfile, "\0") == 0)
    io->error("Missing parameter or invalid value in %s: tempfile is \"%s\"",
              param_file, output->tempfile);

  if (settings->eweight < 0)
    io->error("Missing parameter or invalid value in %s: eng_weight is \"%f\"",
              param_file, settings->eweight);

  if (settings->sweight < 0)
    io->error("Missing parameter or invalid value in %s: stress_weight is \"%f\"",
              param_file, settings->sweight);

  if (output->enable_imd_pot && output->imdpotsteps <= 0)
    io->error("Missing parameter or invalid value in %s: imdpotsteps is \"%d\"",
              param_file, output->imdpotsteps);

  if (output->plotmin < 0)
    io->error("Missing parameter or invalid value in %s: plotmin is \"%f\"",
              param_file, output->plotmin);

  if (potential->enable_cp != 0 && potential->enable_cp != 1)
    io->error("Missing parameter or invalid value in %s: enable_cp is \"%d\"",
              param_file, potential->enable_cp);

  return;
}

void Input::read_potential_file() {
  FILE *infile;
  char  buffer[1024], *res, *str;
  char  temp_char[255];
  int   have_format = 0, have_type = 0, have_elements = 0, end_header = 0;
  int   format, i, size;
  fpos_t after_header;

  // open file
  io->write("Starting to read potential file \"%s\":\n",startpot);
  infile = fopen(startpot, "r");
  if (NULL == infile)
    io->error("Could not open file \"%s\"\n", startpot);

  // read the header
  do {
    // read one line
    res = fgets(buffer, 1024, infile);
    if (NULL == res)
      io->error("Unexpected end of file in >%s<", startpot);
    // check if it is a header line
    if (buffer[0] != '#')
      io->error("Header corrupt in file >%s<", startpot);

    // see if it is the format line
    if (buffer[1] == 'F') {
      // format complete?
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size))
        io->error("Corrupt format header line in file %s", startpot);

      if (size != interaction->force->cols()) {
        sprintf(buffer, "Wrong number of data columns in file \"%s\",\n", startpot);
        io->error("%sshould be %d for %s, but are %d.", buffer,
                  interaction->force->cols(), interaction->type, size);
      }
      have_format = 1;
      potential->init(size);
    } else if (!have_format && buffer[1] != '#') {
      io->error("Format line missing or not first line in potential file.", startpot);
    }

    // read the C line
    if (buffer[1] == 'C') {
      // remove newline at the end of buffer
      if ((str = strchr(buffer + 3, '\n')) != NULL)
        *str = '\0';
      res = strtok(buffer, " ");
      for (int i=0;i<structures->ntypes;i++) {
        res = strtok(NULL, " ");
	if (res == NULL)
	  io->error("Not enough items in #C header line.");
	if (strlen(res) > 10) {
	  strncpy(potential->elements[i], res, 10);
	  potential->elements[i][10] = '\0';
	} else {
	  strcpy(potential->elements[i], res);
	}
      }
      have_elements = 1;
    }

    // read the TYPE line
    else if (buffer[1] == 'T') {
      // remove newline at the end of buffer
      if ((str = strchr(buffer + 3, '\n')) != NULL)
        *str = '\0';
      if (strcmp(utils->tolowercase(buffer + 3), utils->tolowercase(interaction->type)) != 0) {
        io->error("The potential types in your parameter and potential file do not match!\n");
      }
      io->write("- Using %d %s potentials to calculate forces\n", size, interaction->type);
      have_type = 1;
    }

    // invariant potentials
    else if (buffer[1] == 'I') {
      if (have_format && have_type) {
        // check for enough items
        for (i = 0; i < size; i++) {
          str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
          if (str == NULL) {
            io->error("Not enough items in #I header line.");
          } else {
            potential->invar_pot[i] = atoi(str);
            if (potential->invar_pot[i] == 1)
              potential->num_free_pots--;
          }
        }
      } else
        io->error("#I needs to be specified after the #F and #T lines in the potential file \"%s\"", startpot);
    }

    // stop after last header line
    if (buffer[1] == 'E') {
      end_header = 1;
    }

  } while (!end_header);

  /* do we have a format in the header? */
  if (!have_format)
    io->error("Format not specified in header of file \"%s\"", startpot);

  if (!have_elements)
    io->error("The elements are not specified in the potential file (#C tag is missing).");

  fgetpos(infile, &after_header);

  fsetpos(infile, &after_header);
  potential->read_globals(infile);

  fsetpos(infile, &after_header);
  interaction->force->read_additional_data(infile);

  fsetpos(infile, &after_header);
  potential->read_potentials(infile);

  potential->global_params->check_usage();

  io->write("Reading potential file \"%s\" ... done\n", startpot);

  fclose(infile);

  return;
}

void Input::read_config_file() {
  FILE *infile;

  // open file
  io->write("Reading config file \"%s\" ... ", config_file);
  infile = fopen(config_file, "r");
  if (NULL == infile)
    io->error("Could not open file \"s\"\n", config_file);

  structures->init();
  structures->read_config(infile);

  fclose(infile);

  io->write("done\n\n");

  io->write("Read %d configurations (%d with forces, %d with stresses)\n",
            structures->num_conf, structures->using_forces, structures->using_stresses);
  io->write("with a total of %d atoms (", structures->num_atoms);
  for (int i=0; i<structures->ntypes; i++) {
    io->write("%d %s [%.2f%]",
              structures->num_per_type[i], potential->elements[i],
              100. * structures->num_per_type[i]/structures->num_atoms);
    if (i != (structures->ntypes - 1))
      io->write(", ");
  }
  io->write(").\n\n");

  structures->print_mindist();

  return;
}
