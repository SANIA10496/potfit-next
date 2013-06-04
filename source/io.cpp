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

#include "input.h"
#include "io.h"
#include "settings.h"
#include "version.h"

using namespace POTFIT_NS;

IO::IO(POTFIT *ptf) :
  Pointers(ptf),
  write("", std::cout, &logfile),
  warning("Warning", std::cout, &logfile),
  error("Error", std::cerr, &logfile),
  debug("Debug", std::cerr, &logfile),
  screen(0),
  init_done(0),
  write_logfile(0)
{}

IO::~IO() {}

void IO::print_header() {
  write << "This is potfit-next " << POTFIT_VERSION << " compiled on ";
  write << POTFIT_DATE << "." << std::endl << std::endl;

  return;
}

void IO::print_version() {
  MPI::COMM_WORLD.Barrier();
  if (screen) {
    printf("potfit-next version %s\n",POTFIT_VERSION);
    printf("Copyright © 2013 Institute for Theoretical and Applied Physics (ITAP).\n");
    printf("License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>.\n");
    printf("This is free software: you are free to change and redistribute it.\n");
    printf("There is absolutely NO WARRANTY.\n");
    printf("\nContributing authors are listed on <http://potfit.sourceforge.net>.\n");
  }
  pexit(EXIT_SUCCESS);
}

void IO::print_help() {
  MPI::COMM_WORLD.Barrier();
  if (screen) {
    printf("Usage: potfit PARAMETER_FILE\n\n");
    printf("Options:\n");
    printf("\t-h\tprint this screen\n");
    printf("\t-v\tprint version information\n");
  }
  pexit(EXIT_SUCCESS);
}

void IO::set_logfile(const char *filename) {
  if (write_logfile && screen) {
    logfile.open(filename);
    write.init_done(screen);
    warning.init_done(screen);
    error.init_done(screen);
    debug.init_done(screen);
  }
}

void IO::close_logfile() {
  if (write_logfile && screen) {
//    logfile.close();
  }
}

void IO::pexit(int status) {
  // broadcast error
  if (init_done == 1) {
    // MPI::BCast.force finish
  }
  // wait for others to arrive
  MPI::COMM_WORLD.Barrier();

  close_logfile();
  MPI::Finalize();
  exit(status);
}

void IO::set_screen(const int scr) {
  screen = scr;

  return;
}

void IO::set_write_logfile(int i) {
  if (1 == screen)
    write_logfile = i;

  return;
}

int IO::get_write_logfile(void) {
  return write_logfile;
}
