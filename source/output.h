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

#include <string>

#include "pointers.h"

namespace POTFIT_NS {
  class Output : protected Pointers {
  public:
    Output(class POTFIT *);
    ~Output();

    void write_output(void);
    void write_tempfile(void);

    void set_endpot(const std::string &);
    std::string get_endpot(void);

    void set_output_prefix(const std::string &);
    std::string get_output_prefix(void);

    void set_imdpot(const std::string &);
    std::string get_imdpot(void);

    void set_plotfile(const std::string &);
    std::string get_plotfile(void);

    void set_distfile(const std::string &);
    std::string get_distfile(void);

    void set_imdpotsteps(const int &);
    int get_imdpotsteps(void);

    void set_plotmin(const double &);
    double get_plotmin(void);

    void set_tempfile(const std::string &);
    std::string get_tempfile(void);

    void set_enable_pairdist(const int &);
    int get_enable_pairdist(void);

    void set_plotpointfile(const std::string &);
    std::string get_plotpointfile(void);

    void set_lammps_pot(const int &);
    std::string get_lammps_pot(void);

  private:
    void write_distfile(void);
    void write_endpot(void);
    void write_imdpot(void);
    void write_lammpspot(void);
    void write_output_files(void);
    void write_plotfile(void);
    void write_plotpointfile(void);

    void write_output_file_forces(void);
    void write_output_file_energies(void);
    void write_output_file_stresses(void);
    void write_output_file_punishments(void);

    double plotmin;

    int imdpotsteps;
    int enable_distfile;
    int enable_output_files;
    int enable_imd_pot;
    int enable_lammps_pot;
    int enable_plotfile;
    int enable_pairdist;
    int enable_log;

    std::string distfile;
    std::string endpot;
    std::string imdpot;
    std::string lammpspot;
    std::string output_prefix;
    std::string plotfile;
    std::string plotpointfile;
    std::string tempfile;
  };
}

#endif /* PTF_OUTPUT_H */
