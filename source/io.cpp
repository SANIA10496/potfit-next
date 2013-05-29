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
#include <cstdio>
#include <cstdlib>

#include "input.h"
#include "io.h"
#include "settings.h"
#include "version.h"

using namespace POTFIT_NS;

IO::IO(POTFIT *ptf) : Pointers(ptf) {
  init_done = 0;
  screen = 0;
  write_logfile = 0;

  new_write.init("", &logfile);
  new_warning.init("Warning", &logfile);
  new_error.init("Error", &logfile);
  new_debug.init("Debug", &logfile);

  return;
}

IO::~IO() {

  return;
}

void IO::print_header() {
  write("This is potfit-next %s compiled on %s.\n\n",POTFIT_VERSION,POTFIT_DATE);
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
  MPI::Finalize();
  exit(EXIT_SUCCESS);
}

void IO::print_help() {
  MPI::COMM_WORLD.Barrier();
  if (screen) {
    printf("Usage: potfit PARAMETER_FILE\n\n");
    printf("Options:\n");
    printf("\t-h\tprint this screen\n");
    printf("\t-v\tprint version information\n");
  }
  MPI::Finalize();
  exit(EXIT_SUCCESS);
}

void IO::write(const char *msg, ...) {
  va_list ap;

  if (screen) {
    va_start(ap, msg);
    vfprintf(stdout, msg, ap);
    va_end(ap);
    // write to log file if enabled
//    if (write_logfile) {
//      va_start(ap, msg);
//      vfprintf(logfile, msg, ap);
//      va_end(ap);
//    }
  }
}

void IO::write_log(const char *msg, ...) {
  va_list ap;

  if (screen) {
    va_start(ap, msg);
//    if (write_logfile)
//      vfprintf(logfile, msg, ap);
//    va_end(ap);
  }

  return;
}

void IO::write_debug(const char *msg, ...) {
#ifdef DEBUG
  va_list ap;
  char temp[1024];

  if (screen) {
    va_start(ap, msg);
    sprintf(temp, msg, ap);
    va_end(ap);
    fprintf(stderr, "DEBUG: %s", temp);
  }
#endif // DEBUG
  return;
}

void IO::writef(FILE *outfile, const char *msg, ...) {
  va_list ap;

  if (screen) {
    va_start(ap, msg);
    vfprintf(outfile, msg, ap);
    va_end(ap);
  }
}

void IO::warning(const char *msg, ...) {
  va_list ap;

  if (screen) {
    fflush(stdout);
    std::cerr << "\n--> WARNING <--\n";
    va_start(ap, msg);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    std::cerr << "\n";
    fflush(stderr);
//    if (write_logfile) {
//      write_log("\n--> WARNING <--\n");
//      va_start(ap, msg);
//      vfprintf(logfile, msg, ap);
//      va_end(ap);
//    }
  }
}

void IO::error(const char *msg, ...) {
  va_list ap;

  fflush(stdout);
  // broadcast error
  if (init_done == 1) {
    // MPI::BCast.force finish
  }
  // wait for others to arrive
  MPI::COMM_WORLD.Barrier();
  if (screen) {
    fflush(stderr);
    std::cerr << "\n--> ERROR <--\n";
    va_start(ap, msg);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    fflush(stderr);
    std::cerr << std::endl << std::endl;
//    if (write_logfile) {
//      write_log("\n--> ERROR <--\n");
//      va_start(ap, msg);
//      vfprintf(logfile, msg, ap);
//      va_end(ap);
//      write_log("\n\n");
//    }
  }
  close_logfile();
  MPI::Finalize();
  exit(EXIT_FAILURE);
}

void IO::set_logfile(const char *filename) {
  if (write_logfile && screen) {
    logfile.open(filename);
    logfile << "This is potfit-next " << POTFIT_VERSION << " compiled on " << POTFIT_DATE << ".\n\n";
    logfile << "Reading parameter file \"" << input->param_file.c_str() << "\" ... ";
    new_write.init_done(screen);
    new_warning.init_done(screen);
    new_error.init_done(screen);
    new_debug.init_done(screen);
  }
}

void IO::close_logfile() {
  if (write_logfile && screen) {
    logfile.close();
  }
}
