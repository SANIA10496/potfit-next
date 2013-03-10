/****************************************************************
 *
 * io.h: misc io routines
 *
 ****************************************************************
 *
 * Copyright 2002-2012
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

#ifndef PTF_IO_H
#define PTF_IO_H

#include <cstdio>
#include <iostream>

#include "pointers.h"

namespace POTFIT_NS {

  class IO : protected Pointers {
  public:
    IO(class POTFIT *);
    ~IO();

    // output for standard runs
    void print_header();

    // output for flags
    void print_version();
    void print_help();

    // standard output
    void write(const char *, ...);
    void write_log(const char *, ...);

    // output to files
    void writef(FILE *, const char *, ...);

    // misc. stuff
    void warning(const char *, ...);
    void error(const char *, ...);

    void set_logfile(const char *);
    void close_logfile();

    int screen;
    int init_done;
    int write_logfile;

  private:
    FILE *logfile;
  };

}

#endif /* PTF_IO_H */
