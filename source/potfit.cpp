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

#include "potfit.h"
#include "io.h"
#include "input.h"
#include "output.h"
#include "force.h"
#include "optimization.h"
#include "random.h"
#include "config.h"
#include "settings.h"
#include "potential.h"
#include "memory.h"
#include "communication.h"

using namespace POTFIT_NS;

POTFIT::POTFIT(int argc, char **argv) {
  io = new IO(this);
  input = new Input(this, argc, argv);
  output = new Output(this);
  force = new Force(this);
  optimization = new Optimization(this);
  random = new Random(this);
  config = new Config(this);
  settings = new Settings(this);
  potential = new Potential(this);
  memory = new Memory(this);
  communication = new Communication(this);
}

POTFIT::~POTFIT() {
  delete io;
  delete input;
  delete output;
  delete force;
  delete optimization;
  delete random;
  delete config;
  delete settings;
  delete potential;
  delete memory;
  delete communication;
}

void POTFIT::run() {
  io->print_header();
  input->read_parameter_file();
  input->read_potential_file();
  input->read_config_file();
  io->write("Setting up neighbor list ...\n");
  io->write("Running optimization algs ...\n");
  io->write("Calculating error reports ...\n");
  io->write("Cleaning up ...\n");
}
