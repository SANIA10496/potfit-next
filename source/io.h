/****************************************************************
 *
 * io.h: misc io routines
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

#ifndef PTF_IO_H
#define PTF_IO_H

#include <fstream>
#include <iostream>

#include "pointers.h"
#include "io/pstream.h"

namespace POTFIT_NS {

  class IO : protected Pointers {
  public:
    IO(class POTFIT *);
    ~IO();

    // output for standard runs
    void print_header(void);

    // output for command line flags
    void print_version(void);
    void print_help(void);

    // custom streams for output
    PStream write;
    PStream warning;
    PStream error;
    PStream debug;

    // initialize logfile and pass pointer to streams
    void set_logfile(const char *);
    void close_logfile();

    void set_write_logfile(const int &);
    int get_write_logfile(void);

    void pexit(const int &);

    void set_screen(const int &);

  private:
    int screen;
    int init_done;
    int write_logfile;

    std::ofstream logfile;
  };

}

#endif /* PTF_IO_H */
