/****************************************************************
 *
 * optimization.cpp:
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

#include "interaction.h"
#include "io.h"
#include "optimization.h"
#include "potential.h"
#include "settings.h"

#include "opt/list_opt.h"

using namespace POTFIT_NS;

Optimization::Optimization(POTFIT *ptf) :
  Pointers(ptf),
  opt(NULL)
{}

Optimization::~Optimization() {}

void Optimization::run(void) {

  if (1 == settings->get_opt()) {
    if (0 == algorithms.size()) {
      io->warning << "You have enabled the optimization in the parameter file" << std::endl;
      io->warning << "but did not specify any optimization algorithms with the opt_alg keyword!" << std::endl;
      io->write << std::endl;
    }
    io->write << "Starting optimization with the following parameters:" << std::endl;
    if (settings->get_eweight() != 0.0)
      io->write << " - Global energy weight " << settings->get_eweight() << std::endl;
    if (settings->get_sweight() != 0.0)
      io->write << " - Global stress weight "<< settings->get_sweight() << std::endl;
    io->write << " - " << potential->get_num_free_params() << " free parameters" << std::endl << std::endl;

    if (0 == settings->get_myid()) {
      for (unsigned int i=0; i<algorithms.size(); i++) {
        io->write << "Starting " << i + 1 << ". optimization algorithm:" << std::endl;
        opt = init_algorithm(algorithms[i]);
        if (NULL == opt) {
          io->error << "Could not find optimization algorithm \"" << algorithms[i] << "\"." << std::endl;
          io->pexit(EXIT_FAILURE);
        }
        opt->init(params[i]);
        io->write << " Type " << algorithms[i] << " - parameters: ";
        for (unsigned int j=0; j<params[i].size(); j++) {
          io->write << params[i][j];
          if (j < params[i].size()-1)
            io->write << ", ";
        }
        io->write << std::endl << std::endl;
        opt->run();
        io->write << std::endl << "Finished " << i + 1 << ". optimization algorithm." << std::endl << std::endl;
        delete opt;
      }
    } else {
      interaction->calc_forces();
    }
  } else {
    io->write << "Optimization has been disabled in the parameter file." << std::endl << std::endl;
  }

  // Calculate forces at least once for errors
  interaction->calc_forces();

  return;
}

void Optimization::add_algorithm(std::vector<std::string> &tokens) {
  std::vector<std::string> temp;

  algorithms.push_back(tokens[1]);

  for (unsigned int i=2; i<tokens.size(); i++) {
    temp.push_back(tokens[i]);
  }

  params.push_back(temp);

  return;
}

BaseOpt *Optimization::init_algorithm(const std::string &algo) {
  if (algo.compare("none") == 0 or algo.empty())
    return NULL;
#define OPTALG_TYPE
#define OptAlgType(key,Class) \
  else if (algo.compare(#key) == 0) \
    return new Class(ptf);
#include "opt/list_opt.h"
#undef OPTALG_TYPE
  return NULL;
}
