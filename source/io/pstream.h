/****************************************************************
 *
 * pstream.h
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

#ifndef PTF_PSTREAM_H
#define PTF_PSTREAM_H

#include <fstream>
#include <iostream>
#include <ostream>

#include "../pointers.h"

namespace POTFIT_NS {

  class PStream : public std::ofstream {
  public:
    PStream();
    ~PStream();

    void init(const std::string &, std::ofstream *);
    void init_done(int &);

    template <typename T>
    std::ofstream& operator<<(const T&);

    // function that takes a custom stream, and returns it
//    typedef std::ofstream& (*MyStreamManipulator)(std::ofstream&);

    // take in a function with the custom signature
//    std::ofstream& operator<<(MyStreamManipulator manip);

    // define the custom endl for this stream.
    // note how it matches the `MyStreamManipulator`
    // function signature
//    static std::ofstream& endl(std::ofstream& stream);

    // this is the type of std::cout
//    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

    // this is the function signature of std::endl
//    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
//    std::ofstream& operator<<(StandardEndLine manip);

//    std::ofstream& operator<<(const char *);
//    std::ostream& operator<<(std::ostream&(*f)(std::ostream&));

//    std::ofstream& operator<<(bool val);
//    std::ofstream& operator<<(short val);
//    std::ofstream& operator<<(unsigned short val);
//    std::ofstream& operator<<(int val);
//    std::ofstream& operator<<(unsigned int val);
//    std::ofstream& operator<<(long val);
//    std::ofstream& operator<<(unsigned long val);
//    std::ofstream& operator<<(float val);
//    std::ofstream& operator<<(double val);
//    std::ofstream& operator<<(long double val);
//    std::ofstream& operator<<(void* val);

//    std::ostream& operator<<(streambuf* sb);

//    std::ostream& operator<<(ostream& (*pf)(ostream&));
//    std::ostream& operator<<(ios& (*pf)(ios&));
//    std::ostream& operator<<(ios_base& (*pf)(ios_base&));
  private:
    int screen;
    int write_prefix;
    int write_logfile;

    std::ofstream *output;
    std::string prefix;
  };
}

#endif // PTF_PSTREAM_H
