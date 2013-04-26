/****************************************************************
 *
 * table_tab3.cpp:
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

#include <cstdlib>
#include <cstring>

#include "table_tab3.h"

#include "../io.h"
#include "../settings.h"

using namespace POTFIT_NS;

TableTab3::TableTab3(POTFIT *ptf) : Table(ptf) {
}

TableTab3::~TableTab3() {
}

void TableTab3::init(const char *name, int index) {
  return;
}

void TableTab3::read_potential(FILE *a) {
}

int TableTab3::get_number_params(void) {
}

int TableTab3::get_number_free_params(void) {
}

double TableTab3::get_cutoff(void) {
  return end;
}

double TableTab3::get_rmin(void) {
  return begin;
}

void TableTab3::set_params(double *val) {
}

void TableTab3::update_calc_table(void) {

  return;
}

void TableTab3::write_potential(FILE *outfile) {
}

void TableTab3::write_plot(FILE *outfile) {
}

void TableTab3::write_plotpoint(FILE *outfile) {
}
