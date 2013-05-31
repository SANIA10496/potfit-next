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
#include <sstream>
#include <string>

#include "../pointers.h"

namespace POTFIT_NS {

  class PStream: public std::ostream {

  public:
    PStream(const char *, std::ostream &, std::ofstream *, PStream *);
    ~PStream();

    void init_done(int &);

  private:
    // Write a stream buffer that prefixes each line
    class PStreamBuf: public std::stringbuf
    {
    public:
      PStreamBuf(std::ostream& ostr, std::ofstream& lstr, PStream *ps, const char *pref) :
        screen(0),
        write_logfile(0),
        prefix(pref),
        pstream(ps),
        output(ostr),
        logfile(lstr)
      {}

      // When we sync the stream with the output.
      // 1) Output prefix then the buffer
      // 2) Reset the buffer
      // 3) flush the actual output stream we are using.
      virtual int sync ( )
      {
        if (1 == screen) {
          if (!prefix.empty()) {
            output << "[" << prefix << "] " << str();
	    if (1 == write_logfile)
              logfile << "[" << prefix << "] " << str();
            str("");
          } else {
            output << str();
	    if (1 == write_logfile)
              logfile << str();
            str("");
          }
          output.flush();
          logfile.flush();
          return 0;
        } else
          return 0;
      }

      void init_done(int &);

    private:
      int screen;
      int write_logfile;
      std::string prefix;

      PStream *pstream;
      std::ostream& output;
      std::ofstream& logfile;
    };

    PStreamBuf outbuff;
  };
}

#endif // PTF_PSTREAM_H
