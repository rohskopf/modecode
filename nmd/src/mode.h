/*
 mode.h

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <string>
#include "mpi.h"

#include <iostream>
#include <new>
#include <cstdlib>
#include "pointers.h"

using namespace std;

namespace EM3_NS
{
  class Mode: protected Pointers
  {
  public:
    Mode(class EM3 *);
    ~Mode();

    FILE * fh_debug; // Debug file handle
    FILE * fh_em; // positions output in .xyz format, from realspace integration
    FILE * fh_fv; // positions output in .xyz format, from realspace integration

    void initialize();
    void mode2Atom(); // Convert mode quantities to atom quantities
    void atom2Mode(); // Convert atom quantities to mode quantities
    void calcModeEnergy(); // calclate mode energies

    int nmodes; // number of modes
    int natoms; // number of atoms
    double *xm; // Mode coordinates
    double **xmta; // Mode coordinates converted to atom coordinates
    double *vm; // Mode velocities
    double *fm; // Mode forces (accelerations)
    double *fm2; // Mode forces (accelerations) (2nd order)
    double *fm3; // Mode forces (accelerations) (3rd order)
    double *fm4; // Mode forces (accelerations) (4th order)
    double *em; // mode energies
    double fv; // mode power transfer for a single mode

    double **xa; // Atom positions
    double **va; // atom velocities
    double **fa; // atom forces
    double **fa2; // atom forces (2nd order)
    double **fa3; // atom forces (3rd order)
    double **fa4; // atom forces (4th order)

    double **fmta; // Mode forces converted to atom forces
    double **fmta2; // Mode forces converted to atom forces (2nd order)
    double **fmta3; // Mode forces converted to atom forces (3rd order)
    double **fmta4; // Mode forces converted to atom forces (4th order)
    

  };
}

