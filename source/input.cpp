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

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <sstream>

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

Input::Input(POTFIT *ptf, int argc, char **argv) :
  Pointers(ptf),
  curline(0),
  enable_maxch_file(0)
{
  settings->set_myid(MPI::COMM_WORLD.Get_rank());
  settings->set_num_cpus(MPI::COMM_WORLD.Get_size());

  if (settings->get_myid() == 0) {
    io->set_screen(1);
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
      io->error << "Argument " << argv[i] << " unknown." << std::endl;
      io->pexit(EXIT_FAILURE);
    }
  }
  if (argc > 2) {
    io->error << "Please provide single filename." << std::endl;
    io->print_help();
    io->pexit(EXIT_FAILURE);
  }
  if (1 == version)
    io->print_version();

  param_file.assign(argv[1]);

  return;
}

Input::~Input() {}

void Input::read_parameter_file() {
  int counter = 0;
  std::ifstream infile;
  std::istringstream iss;
  std::string line;
  std::vector<std::string> tokens;

  io->write << "Reading parameter file \"" << param_file << "\" ... ";

  infile.exceptions (std::ifstream::failbit);

  try {
    infile.open(param_file.c_str(), std::ios::in);
  } catch (std::ifstream::failure &e) {
    io->error << "Could not open parameter file \"" << param_file << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  infile.exceptions ( std::ifstream::goodbit );

  while (!infile.eof()) {
    getline(infile,line);
    counter++;

    if (line[0] == '#' || line[0] == ' ' || line.empty()) {
      continue;
    }

    tokens.clear();
    iss.clear();
    iss.str(line);
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter<std::vector<std::string> >(tokens));

    // number of atom types
    if (tokens[0].compare("ntypes") == 0) {
      structures->set_ntypes(get_param_int("ntypes", tokens[1]));
    }
    // file with start potential
    else if (tokens[0].compare("startpot") == 0) {
      startpot = get_param_str("startpot", tokens[1]);
    }
    // file for end potential
    else if (tokens[0].compare("endpot") == 0) {
      output->set_endpot(get_param_str("endpot", tokens[1]));
    }
    // prefix for all output files
    else if (tokens[0].compare("output_prefix") == 0) {
      output->set_output_prefix(get_param_str("output_prefix", tokens[1]));
    }
    // file for IMD potential
    else if (tokens[0].compare("imdpot") == 0) {
      output->set_imdpot(get_param_str("imdpot", tokens[1]));
    }
    // file for plotting
    else if (tokens[0].compare("plotfile") == 0) {
      output->set_plotfile(get_param_str("plotfile", tokens[1]));
    }
    // file for maximal change */
    else if (tokens[0].compare("maxchfile") == 0) {
      maxchfile = get_param_str("maxchfile", tokens[1]);
      enable_maxch_file = 1;
    }
    // file for pair distribution
    else if (tokens[0].compare("distfile") == 0) {
      output->set_distfile(get_param_str("distfile", tokens[1]));
    }
    // number of steps in IMD potential
    else if (tokens[0].compare("imdpotsteps") == 0) {
      output->set_imdpotsteps(get_param_int("imdpotsteps", tokens[1]));
    }
    // minimum for plotfile
    else if (tokens[0].compare("plotmin") == 0) {
      output->set_plotmin(get_param_dbl("plotmin", tokens[1]));
    }
    // exclude chemical potential from energy calculations
    else if (tokens[0].compare("enable_cp") == 0) {
      potential->set_enable_cp(get_param_int("enable_cp", tokens[1]));
    }
    // how far should the imd pot be extended
    else if (tokens[0].compare("extend") == 0) {
      settings->set_extend(get_param_dbl("extend", tokens[1]));
    }
    // file with atom configuration
    else if (tokens[0].compare("config") == 0) {
      config_file = get_param_str("config", tokens[1]);
    }
    // Optimization flag
    else if (tokens[0].compare("opt") == 0) {
      settings->set_opt(get_param_int("opt", tokens[1]));
    }
    // break flagfile
    else if (tokens[0].compare("flagfile") == 0) {
      utils->set_flagfile(get_param_str("flagfile", tokens[1]));
    }
    // write radial pair distribution
    else if (tokens[0].compare("write_pair") == 0) {
      output->set_enable_pairdist(get_param_int("write_pair", tokens[1]));
    }
    // plotpoint file
    else if (tokens[0].compare("plotpointfile") == 0) {
      output->set_plotpointfile(get_param_str("plotpointfile", tokens[1]));
    }
    // temporary potential file
    else if (tokens[0].compare("tempfile") == 0) {
      output->set_tempfile(get_param_str("tempfile", tokens[1]));
    }
    // seed for RNG
    else if (tokens[0].compare("seed") == 0) {
      random->set_seed(get_param_int("seed", tokens[1]));
    }
    // Scaling Constant for APOT Punishment
    else if (tokens[0].compare("apot_punish") == 0) {
      settings->set_apot_punish(get_param_dbl("apot_punish", tokens[1]));
    }
    // Energy Weight
    else if (tokens[0].compare("eng_weight") == 0) {
      settings->set_eweight(get_param_dbl("eng_weight", tokens[1]));
    }
    // Energy Weight
    else if (tokens[0].compare("stress_weight") == 0) {
      settings->set_sweight(get_param_dbl("stress_weight", tokens[1]));
    }
    // write final potential in lammps format
    else if (tokens[0].compare("write_lammps") == 0) {
      output->set_lammps_pot(get_param_int("write_lammps", tokens[1]));
    }
    // cutoff-radius for long-range interactions
    else if (tokens[0].compare("dp_cut") == 0) {
      settings->set_dp_cut(get_param_dbl("dp_cut", tokens[1]));
    }
    // dipole iteration precision
    else if (tokens[0].compare("dp_tol") == 0) {
      settings->set_dp_tol(get_param_dbl("dp_tol", tokens[1]));
    }
    // mixing parameter for damping dipole iteration loop
    else if (tokens[0].compare("dp_mix") == 0) {
      settings->set_dp_mix(get_param_dbl("dp_mix", tokens[1]));
    }
    // log file
    else if (tokens[0].compare("write_log") == 0) {
      io->set_write_logfile(get_param_int("write_log", tokens[1]));
    }
    // interaction type
    else if (tokens[0].compare("interaction") == 0) {
      interaction->set_type(get_param_str("interaction", tokens[1]));
    }
    // optimization algorithm
    else if (tokens[0].compare("opt_alg") == 0) {
      optimization->add_algorithm(tokens);
    }
    // if tag could not be identified
    else {
      std::cerr << "Unknown tag " << tokens[0] << " ignored!" << std::endl;
    }

    // end of while loop
  }

  infile.close();

  io->set_logfile("potfit.log");
  io->write << "done" << std::endl;
  check_params();
  interaction->init();

  return;
}

double Input::get_param_dbl(const std::string &name, const std::string &token) {
  return atof(token.c_str());
}

int Input::get_param_int(const std::string &name, const std::string &token) {
  return atoi(token.c_str());
}

std::string Input::get_param_str(const std::string &name, const std::string &token) {
  return token;
}

/****************************************************************
 *
 *  check all parameters for reasonable values and completeness
 *
 ****************************************************************/

void Input::check_params(void)
{
  if (structures->get_ntypes() <= 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": ntypes is \"" << structures->get_ntypes() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (startpot.empty()) {
    io->error << "Missing parameter or invalid value in " << param_file << ": startpot is \"" << startpot << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (output->get_endpot().empty()) {
    io->warning << "endpot is missing in " << param_file << ", setting it to " << startpot << "_end" << std::endl;
    output->get_endpot() = startpot + "_end";
  }

  if (interaction->get_type().empty()) {
    io->error << "Error in " << param_file << " - 'interaction' keyword is missing!" << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (strcmp(config_file.c_str(), "\0") == 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": config is \"" << config_file << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (output->get_tempfile().empty()) {
    io->error << "Missing parameter or invalid value in " << param_file << ": tempfile is \"" << output->get_tempfile() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (settings->get_eweight() < 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": eng_weight is \"" << settings->get_eweight() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (settings->get_sweight() < 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": stress_weight is \"" << settings->get_sweight() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (!output->get_imdpot().empty() && output->get_imdpotsteps() <= 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": imdpotsteps is \"" << output->get_imdpotsteps() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (output->get_plotmin() < 0) {
    io->error << "Missing parameter or invalid value in " << param_file << ": plotmin is \"" << output->get_plotmin() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (potential->get_enable_cp() != 0 && potential->get_enable_cp() != 1) {
    io->error << "Missing parameter or invalid value in " << param_file << ": enable_cp is \"" << potential->get_enable_cp() << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  return;
}

void Input::read_potential_file() {
  FILE *infile;
  char  buffer[1024], *res, *str;
  int   have_format = 0, have_type = 0, have_elements = 0, end_header = 0;
  int   format, i, size;
  fpos_t after_header;

  // open file
  io->write << "Starting to read potential file \"" << startpot << "\":" << std::endl;
  infile = fopen(startpot.c_str(), "r");
  if (NULL == infile) {
    io->error << "Could not open file \"" << startpot << "\"" << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  // read the header
  do {
    // read one line
    res = fgets(buffer, 1024, infile);
    if (NULL == res) {
      io->error << "Unexpected end of file in \""<< startpot << "\"" << std::endl;
      io->pexit(EXIT_FAILURE);
    }
    // check if it is a header line
    if (buffer[0] != '#') {
      io->error << "Header corrupt in file \"" << startpot << "\"" << std::endl;
      io->pexit(EXIT_FAILURE);
    }

    // see if it is the format line
    if (buffer[1] == 'F') {
      // format complete?
      if (2 != sscanf((const char *)(buffer + 2), "%d %d", &format, &size)) {
        io->error << "Corrupt format header line in file " << startpot << std::endl;
        io->pexit(EXIT_FAILURE);
      }

      if (size != interaction->force->cols()) {
        io->error << "Wrong number of data columns in file \"" << startpot << "\"," << std::endl;
        io->error << "should be " << interaction->force->cols() << " for " << interaction->get_type() << ", but are " << size << "." << std::endl;
        io->pexit(EXIT_FAILURE);
      }
      have_format = 1;
      potential->init(size);
    } else if (!have_format && buffer[1] != '#') {
      io->error << "Format line missing or not first line in potential file." << std::endl;
      io->pexit(EXIT_FAILURE);
    }

    // read the C line
    if (buffer[1] == 'C') {
      // remove newline at the end of buffer
      if ((str = strchr(buffer + 3, '\n')) != NULL)
        *str = '\0';
      res = strtok(buffer, " ");
      for (int i=0; i<structures->get_ntypes(); i++) {
        res = strtok(NULL, " ");
        if (res == NULL) {
          io->error << "Not enough items in #C header line." << std::endl;
          io->pexit(EXIT_FAILURE);
        }
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
      io->write << "- Using " << size << " " << interaction->get_type() << " potentials to calculate forces" << std::endl;
      have_type = 1;
    }

    // invariant potentials
    else if (buffer[1] == 'I') {
      if (have_format && have_type) {
        // check for enough items
        for (i = 0; i < size; i++) {
          str = strtok(((i == 0) ? buffer + 2 : NULL), " \t\r\n");
          if (str == NULL) {
            io->error << "Not enough items in #I header line." << std::endl;
            io->pexit(EXIT_FAILURE);
          } else {
            potential->invar_pot[i] = atoi(str);
            if (potential->invar_pot[i] == 1)
              potential->num_free_pots--;
          }
        }
      } else {
        io->error << "#I needs to be specified after the #F and #T lines in the potential file \"" << startpot << "\"." << std::endl;
        io->pexit(EXIT_FAILURE);
      }
    }

    // stop after last header line
    if (buffer[1] == 'E') {
      end_header = 1;
    }

  } while (!end_header);

  /* do we have a format in the header? */
  if (!have_format) {
    io->error << "Format not specified in header of file \"" << startpot << "\"." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  if (!have_elements) {
    io->error << "The elements are not specified in the potential file (#C tag is missing)." << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  fgetpos(infile, &after_header);

  fsetpos(infile, &after_header);
  potential->read_globals(infile);

  fsetpos(infile, &after_header);
  interaction->force->read_additional_data(infile);

  fsetpos(infile, &after_header);
  potential->read_potentials(infile);

  potential->global_params->check_usage();

  io->write << "Reading potential file \"" << startpot << "\" ... done" << std::endl;

  fclose(infile);

  return;
}

void Input::read_config_file() {
  FILE *infile;

  // open file
  io->write << "Reading config file \"" << config_file << "\" ... ";
  infile = fopen(config_file.c_str(), "r");
  if (NULL == infile) {
    io->write << std::endl;
    io->error << "Could not open file \"" << config_file << "\"" << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  structures->init();
  structures->read_config(infile);

  fclose(infile);

  io->write << "done" << std::endl << std::endl;

  io->write << "Read " << structures->get_num_total_configs() << " configurations ";
  io->write << "(" << structures->using_forces << " with forces, " << structures->using_stresses << " with stresses)" << std::endl;
  io->write << "with a total of " << structures->get_num_total_atoms() << " atoms (";
  for (int i=0; i<structures->get_ntypes(); i++) {
    io->write << structures->num_per_type[i] << " " << potential->elements[i] << " [";
    io->write << 100. * structures->num_per_type[i]/structures->get_num_total_atoms() << "%]";
    if (i != (structures->get_ntypes() - 1))
      io->write << ", ";
  }
  io->write << ")." << std::endl << std::endl;

  structures->print_mindist();

  return;
}
