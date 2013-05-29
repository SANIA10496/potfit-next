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

#include <cstdlib>

#include "force.h"
#include "input.h"
#include "interaction.h"
#include "io.h"
#include "potential.h"
#include "optimization.h"
#include "output.h"
#include "random.h"
#include "settings.h"
#include "structures.h"
#include "types.h"
#include "utils.h"

using namespace POTFIT_NS;

Input::Input(POTFIT *ptf, int argc, char **argv) : Pointers(ptf) {
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

  param_file.assign(argv[1]);

  return;
}

Input::~Input() {
  return;
}

void Input::read_parameter_file() {
  char  buffer[1024];
  int seed = 0;
  std::string temp;

  io->write("Reading parameter file \"%s\" ... ",param_file.c_str());

  if (infile.is_open())
    io->error("File is already open?\n");

  try {
    infile.open(param_file.c_str(), std::ios::in);
  } catch (std::ifstream::failure e) {
    io->error("Could not open parameter file \"%s\".\n",param_file.c_str());
  }

  while (infile.good()) {
    infile >> temp;

    if (temp[0] == '#') {
      infile.getline(buffer, 1024);
      continue;
    }

    // number of atom types
    if (temp.compare("ntypes") == 0) {
      get_param_int("ntypes", structures->ntypes);
    }
    // file with start potential
    else if (temp.compare("startpot") == 0) {
      get_param_str("startpot", input->startpot);
    }
    // file for end potential
    else if (temp.compare("endpot") == 0) {
      get_param_str("endpot", output->endpot);
    }
    // prefix for all output files
    else if (temp.compare("output_prefix") == 0) {
      get_param_str("output_prefix", output->output_prefix);
      if (output->output_prefix.compare("") != 0)
        output->enable_output_files = 1;
    }
    // file for IMD potential
    else if (temp.compare("imdpot") == 0) {
      get_param_str("imdpot", output->imdpot);
      output->enable_imd_pot = 1;
    }
    // file for plotting
    else if (temp.compare("plotfile") == 0) {
      get_param_str("plotfile", output->plotfile);
      output->enable_plot_file = 1;
    }
    // file for maximal change */
    else if (temp.compare("maxchfile") == 0) {
      get_param_str("maxchfile", input->maxchfile);
      input->enable_maxch_file = 1;
    }
    // file for pair distribution
    else if (temp.compare("distfile") == 0) {
      get_param_str("distfile", output->distfile);
      output->enable_distfile = 1;
    }
    // number of steps in IMD potential
    else if (temp.compare("imdpotsteps") == 0) {
      get_param_int("imdpotsteps", output->imdpotsteps);
    }
    // minimum for plotfile
    else if (temp.compare("plotmin") == 0) {
      get_param_dbl("plotmin", output->plotmin);
    }
    // exclude chemical potential from energy calculations
    else if (temp.compare("enable_cp") == 0) {
      get_param_int("enable_cp", potential->enable_cp);
    }
    // how far should the imd pot be extended
    else if (temp.compare("extend") == 0) {
      get_param_dbl("extend", settings->extend);
    }
    // file with atom configuration
    else if (temp.compare("config") == 0) {
      get_param_str("config", input->config_file);
    }
    // Optimization flag
    else if (temp.compare("opt") == 0) {
      get_param_int("opt", settings->opt);
    }
    // break flagfile
    else if (temp.compare("flagfile") == 0) {
      get_param_str("flagfile", utils->flagfile);
    }
    // write radial pair distribution
    else if (temp.compare("write_pair") == 0) {
      get_param_int("write_pair", output->enable_pair_dist);
    }
    // plotpoint file
    else if (temp.compare("plotpointfile") == 0) {
      get_param_str("plotpointfile", output->plotpointfile);
    }
    // temporary potential file
    else if (temp.compare("tempfile") == 0) {
      get_param_str("tempfile", output->tempfile);
    }
    // seed for RNG
    else if (temp.compare("seed") == 0) {
      get_param_int("seed", seed);
      random->set_seed(seed);
    }
    // Scaling Constant for APOT Punishment
    else if (temp.compare("apot_punish") == 0) {
      get_param_dbl("apot_punish", settings->apot_punish_value);
    }
    // Energy Weight
    else if (temp.compare("eng_weight") == 0) {
      get_param_dbl("eng_weight", settings->eweight);
    }
    // Energy Weight
    else if (temp.compare("stress_weight") == 0) {
      get_param_dbl("stress_weight", settings->sweight);
    }
    // write final potential in lammps format
    else if (temp.compare("write_lammps") == 0) {
      get_param_int("write_lammps", output->enable_lammps_pot);
    }
    // cutoff-radius for long-range interactions
    else if (temp.compare("dp_cut") == 0) {
      get_param_dbl("dp_cut", settings->dp_cut);
    }
    // dipole iteration precision
    else if (temp.compare("dp_tol") == 0) {
      get_param_dbl("dp_tol", settings->dp_tol);
    }
    // mixing parameter for damping dipole iteration loop
    else if (temp.compare("dp_mix") == 0) {
      get_param_dbl("dp_mix", settings->dp_mix);
    }
    // log file
    else if (temp.compare("write_log") == 0) {
      get_param_int("write_log", io->write_logfile);
    }
    // interaction type
    else if (temp.compare("interaction") == 0) {
      get_param_str("interaction", interaction->type);
    }
    // optimization algorithm
    else if (temp.compare("opt_alg") == 0) {
      optimization->add_algorithm(infile);
    }
    // if tag could not be identified
    else {
      std::cerr << "Unknown tag " << temp << " ignored!" << std::endl;
    }

  // end of while loop
  }

  infile.close();

  io->set_logfile("potfit.log");
  io->new_write << "done." << std::endl;
  check_params();
  interaction->init();

  return;
}

void Input::get_param_dbl(const std::string &name, double &dest) {
  std::string temp;

  if (!infile.is_open())
    io->error("File is not open?\n");

  infile >> temp;

  dest = atof(temp.c_str());

  return;
}

void Input::get_param_int(const std::string &name, int &dest) {
  std::string temp;

  if (!infile.is_open())
    io->error("File is not open?\n");

  infile >> temp;

  dest = atoi(temp.c_str());

  return;
}

void Input::get_param_str(const std::string &name, std::string &dest) {
  if (!infile.is_open())
    io->error("File is not open?\n");

  infile >> dest;

  return;
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
              param_file.c_str(), structures->ntypes);

  if (strcmp(startpot.c_str(), "\0") == 0)
    io->error("Missing parameter or invalid value in %s : startpot is \"%s\"",
              param_file.c_str(), startpot.c_str());

  if (strcmp(output->endpot.c_str(), "\0") == 0) {
    io->warning("endpot is missing in %s, setting it to %s_end", param_file.c_str(), startpot.c_str());
    output->endpot = startpot + "_end";
  }

  if (strcmp(interaction->type.c_str(), "\0") == 0)
    io->error("Error in %s - 'interaction' keyword is missing!", param_file.c_str());

  if (strcmp(config_file.c_str(), "\0") == 0)
    io->error("Missing parameter or invalid value in %s: config is \"%s\"",
              param_file.c_str(), config_file.c_str());

  if (strcmp(output->tempfile.c_str(), "\0") == 0)
    io->error("Missing parameter or invalid value in %s: tempfile is \"%s\"",
              param_file.c_str(), output->tempfile.c_str());

  if (settings->eweight < 0)
    io->error("Missing parameter or invalid value in %s: eng_weight is \"%f\"",
              param_file.c_str(), settings->eweight);

  if (settings->sweight < 0)
    io->error("Missing parameter or invalid value in %s: stress_weight is \"%f\"",
              param_file.c_str(), settings->sweight);

  if (output->enable_imd_pot && output->imdpotsteps <= 0)
    io->error("Missing parameter or invalid value in %s: imdpotsteps is \"%d\"",
              param_file.c_str(), output->imdpotsteps);

  if (output->plotmin < 0)
    io->error("Missing parameter or invalid value in %s: plotmin is \"%f\"",
              param_file.c_str(), output->plotmin);

  if (potential->enable_cp != 0 && potential->enable_cp != 1)
    io->error("Missing parameter or invalid value in %s: enable_cp is \"%d\"",
              param_file.c_str(), potential->enable_cp);

  return;
}

void Input::read_potential_file() {
  FILE *infile;
  char  buffer[1024], *res, *str;
  int   have_format = 0, have_type = 0, have_elements = 0, end_header = 0;
  int   format, i, size;
  fpos_t after_header;

  // open file
  io->write("Starting to read potential file \"%s\":\n",startpot.c_str());
  infile = fopen(startpot.c_str(), "r");
  if (NULL == infile)
    io->error("Could not open file \"%s\"\n", startpot.c_str());

  // read the header
  do {
    // read one line
    res = fgets(buffer, 1024, infile);
    if (NULL == res)
      io->error("Unexpected end of file in >%s<", startpot.c_str());
    // check if it is a header line
    if (buffer[0] != '#')
      io->error("Header corrupt in file >%s<", startpot.c_str());

    // see if it is the format line
    if (buffer[1] == 'F') {
      // format complete?
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size))
        io->error("Corrupt format header line in file %s", startpot.c_str());

      if (size != interaction->force->cols()) {
        sprintf(buffer, "Wrong number of data columns in file \"%s\",\n", startpot.c_str());
        io->error("%sshould be %d for %s, but are %d.", buffer,
                  interaction->force->cols(), interaction->type.c_str(), size);
      }
      have_format = 1;
      potential->init(size);
    } else if (!have_format && buffer[1] != '#') {
      io->error("Format line missing or not first line in potential file.", startpot.c_str());
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
      //      TODO port to C++
//      if (strcmp(utils->tolowercase(buffer + 3), utils->tolowercase(interaction->type.c_str())) != 0) {
//        io->error("The potential types in your parameter and potential file do not match!\n");
//      }
      io->write("- Using %d %s potentials to calculate forces\n", size, interaction->type.c_str());
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
        io->error("#I needs to be specified after the #F and #T lines in the potential file \"%s\"", startpot.c_str());
    }

    // stop after last header line
    if (buffer[1] == 'E') {
      end_header = 1;
    }

  } while (!end_header);

  /* do we have a format in the header? */
  if (!have_format)
    io->error("Format not specified in header of file \"%s\"", startpot.c_str());

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

  io->write("Reading potential file \"%s\" ... done\n", startpot.c_str());

  fclose(infile);

  return;
}

void Input::read_config_file() {
  FILE *infile;

  // open file
  io->write("Reading config file \"%s\" ... ", config_file.c_str());
  infile = fopen(config_file.c_str(), "r");
  if (NULL == infile)
    io->error("Could not open file \"s\"\n", config_file.c_str());

  structures->init();
  structures->read_config(infile);

  fclose(infile);

  io->write("done\n\n");

  io->write("Read %d configurations (%d with forces, %d with stresses)\n",
            structures->get_num_total_configs(), structures->using_forces, structures->using_stresses);
  io->write("with a total of %d atoms (", structures->get_num_total_atoms());
  for (int i=0; i<structures->ntypes; i++) {
    io->write("%d %s [%.2f%]",
              structures->num_per_type[i], potential->elements[i],
              100. * structures->num_per_type[i]/structures->get_num_total_atoms());
    if (i != (structures->ntypes - 1))
      io->write(", ");
  }
  io->write(").\n\n");

  structures->print_mindist();

  return;
}
