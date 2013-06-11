/****************************************************************
 *
 * utils.cpp:
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

#include "potential.h"
#include "utils.h"

using namespace POTFIT_NS;

Utils::Utils(POTFIT *ptf) :
  Pointers(ptf),
  t_begin(0),
  t_end(0)
{}

Utils::~Utils() {}

char *Utils::tolowercase(char *str) {
  int differ = 'A'-'a';
  char ch;
  int ii = strlen(str);

  for (int i=0; i <ii; i++)
  {
    strncpy(&ch,str+i,1);
    if (ch>='A' && ch<='Z')
    {
      ch = ch-differ;
      memcpy(str+i,&ch,1);
    }
  }
  return str;
}

double Utils::vect_dist(vector a, vector b) {
  return sqrt(square(b.x-a.x)+square(b.y-a.y)+square(b.z-a.z));
}


void Utils::start_timer(void) {
  time(&t_begin);

  return;
}

void Utils::end_timer(void) {
  time(&t_end);

  return;
}

int Utils::timediff(void) {
  return (int)difftime(t_end, t_begin);
}

void Utils::set_flagfile(std::string str) {
  flagfile = str;

  return;
}

int Utils::check_for_flagfile(void) {
  return 0;
}
