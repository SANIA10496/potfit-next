/****************************************************************
 *
 * pstream.cpp:
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

#include "pstream.h"

using namespace POTFIT_NS;

PStream::PStream() {
  screen = 0;
  write_prefix = 0;
  write_logfile = 0;

  return;
}

PStream::~PStream() {
  return;
}

void PStream::init(const std::string &pref, std::ofstream *log) {
  prefix.assign(pref);
  write_prefix = 1;

  output = log;

  return;
}

void PStream::init_done(int &scr) {
  screen = scr;
  write_logfile = 1;

  return;
}

template <typename T>
std::ofstream& PStream::operator<<(const T& x) {
  (*output) << x;
  std::cout << x;

  return (*this);
}


//std::ofstream& PStream::operator<<(StandardEndLine manip) {
//  std::cout << manip;

//  return *this;
//}

//std::ofstream& PStream::operator<<(const std::string &text) {
//  (*output) << text;
//  std::cout << text;

//  return (*this);
//}

//std::ofstream& PStream::operator<<(const char *text) {
//  (*output) << text;
//  std::cout << text;

//  return (*this);
//}

//std::ostream& PStream::operator<<(std::ostream&(*f)(std::ostream&)) {
//  (*output) << f;
//  std::cout << f;

//  return (*this);
//}

