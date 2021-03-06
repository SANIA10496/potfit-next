/****************************************************************
 *
 * potfit.cpp: Contains main potfit program
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

//#include <stdexcept>

#include "potfit.h"
#include "io.h"
#include "input.h"
#include "output.h"
#include "force.h"
#include "interaction.h"
#include "optimization.h"
#include "random.h"
#include "structures.h"
#include "settings.h"
#include "potential.h"
#include "communication.h"
#include "utils.h"
#include "error.h"

using namespace POTFIT_NS;

POTFIT::POTFIT(int argc, char **argv) {
  io = new IO(this);
  output = new Output(this);
  interaction = new Interaction(this);
  optimization = new Optimization(this);
  random = new Random(this);
  structures = new Structures(this);
  settings = new Settings(this);
  potential = new Potential(this);
  communication = new Communication(this);
  utils = new Utils(this);
  input = new Input(this, argc, argv);
  error = new Error(this);

  return;
}

POTFIT::~POTFIT() {
  delete io;
  delete input;
  delete output;
  delete interaction;
  delete optimization;
  delete random;
  delete structures;
  delete settings;
  delete potential;
  delete communication;
  delete utils;
  delete error;

  return;
}

void POTFIT::run() {
  // print header before IO is established
  io->print_header();

  // set up MPI
  communication->init();

  // read all input files
  input->read_parameter_file();
  input->read_potential_file();
  input->read_config_file();

  // perform optimization
  if (0 == settings->get_myid())
    utils->start_timer();
  optimization->run();
  if (0 == settings->get_myid()) {
    utils->end_timer();
    // write potentials to disk
    output->write_output();
    // write error report
    error->write_report();
  }

  return;
}
