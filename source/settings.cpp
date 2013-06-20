/****************************************************************
 *
 * settings.cpp:
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

#include <cstring>

#include "settings.h"

using namespace POTFIT_NS;

Settings::Settings(POTFIT *ptf) :
  Pointers(ptf),
  myid(-1),
  num_cpus(0),
  opt(1),
  eweight(-1.0),
  sweight(-1.0),
  extend(0),
  apot_punish_value(0),
  d_eps(0),
  dp_cut(0),
  dp_tol(0),
  dp_mix(0)
{}

Settings::~Settings() {
  return;
}

void Settings::set_myid(const int &id) {
  myid = id;

  return;
}

const int Settings::get_myid(void) {
  return myid;
}

void Settings::set_num_cpus(const int &n) {
  num_cpus = n;

  return;
}

const int Settings::get_num_cpus(void) {
  return num_cpus;
}

void Settings::set_extend(const double &x) {
  extend = x;

  return;
}

const double Settings::get_extend(void) {
  return extend;
}

void Settings::set_opt(const int &i) {
  opt = i;

  return;
}

const int Settings::get_opt(void) {
  return opt;
}

void Settings::set_apot_punish(const double &x) {
  apot_punish_value = x;

  return;
}

const double Settings::get_apot_punish(void) {
  return apot_punish_value;
}

void Settings::set_eweight(const double &x) {
  eweight = x;

  return;
}

const double Settings::get_eweight(void) {
  return eweight;
}

void Settings::set_sweight(const double &x) {
  sweight = x;

  return;
}

const double Settings::get_sweight(void) {
  return sweight;
}

void Settings::set_deps(const double &x) {
  d_eps = x;

  return;
}

const double Settings::get_deps(void) {
  return d_eps;
}

void Settings::set_dp_cut(const double &x) {
  dp_cut = x;

  return;
}

const double Settings::get_dp_cut(void) {
  return dp_cut;
}

void Settings::set_dp_tol(const double &x) {
  dp_tol = x;

  return;
}

const double Settings::get_dp_tol(void) {
  return dp_tol;
}

void Settings::set_dp_mix(const double &x) {
  dp_mix = x;

  return;
}

const double Settings::get_dp_mix(void) {
  return dp_mix;
}
