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

#ifndef PTF_OUTPUT_H
#define PTF_OUTPUT_H

#include "pointers.h"

namespace POTFIT_NS {

  class Output : protected Pointers {
  public:
    Output(class POTFIT *);
    ~Output();

    void write_output(void);
    void write_tempfile(void);

    char distfile[255];
    char endpot[255];
    char imdpot[255];
    char lammpspot[255];
    char output_prefix[255];
    char plotfile[255];
    char plotpointfile[255];
    char tempfile[255];

    int imdpotsteps;
    int enable_distfile;
    int enable_output_files;
    int enable_imd_pot;
    int enable_lammps_pot;
    int enable_plot_file;
    int enable_pair_dist;
    int enable_log;

    double plotmin;
  private:
    void write_distfile(void);
    void write_endpot(void);
    void write_imdpot(void);
    void write_lammpspot(void);
    void write_output_files(void);
    void write_plotfile(void);
    void write_plotpointfile(void);
  };

}

#endif /* PTF_OUTPUT_H */
