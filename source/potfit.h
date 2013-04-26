/****************************************************************
 *
 * potfit.h: Contains main potfit class
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

#ifndef PTF_POTFIT_H
#define PTF_POTFIT_H

#define POT_STEPS 500

namespace POTFIT_NS {

  class POTFIT {
  public:
    class Settings *settings; 		// Global settings
    class IO *io;  			// I/O for the console
    class Input *input; 		// Input from files
    class Output *output; 		// Output to files
    class Interaction *interaction;	// Force routines
    class Optimization *optimization; 	// Optimization algorithms
    class Random *random; 		// RNG
    class Structures *structures;	// Reference structures
    class Potential *potential; 	// Potential
    class Memory *memory; 		// Memory management
    class Communication *communication; // Process communication
    class Utils *utils; 		// Utilities

    POTFIT(int, char **);
    ~POTFIT();

    // main potfit procedure, calls all subroutines
    void run();
  };

}

#endif /* POTFIT_H */
