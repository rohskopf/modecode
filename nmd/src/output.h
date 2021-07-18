/*
 output.h

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
  class Output: protected Pointers
  {
  public:
    Output(class EM3 *);
    ~Output();

    FILE * fh_xyz_a; // positions output in .xyz format, from realspace integration
    FILE * fh_xyz_m; // positions output in .xyz format, from modespace integration
    FILE * fh_forces; // forces
    FILE * fh_forces2; // forces (2nd order)
    FILE * fh_forces3; // forces (3rd order)
    FILE * fh_forces4; // forces (3rd order)
    FILE * fh_forces_mta; // mode forces converted to atom forces
    FILE * fh_forces_mta2; // mode forces converted to atom forces (2nd order)
    FILE * fh_forces_mta3; // mode forces converted to atom forces (3rd order)
    FILE * fh_forces_mta4; // mode forces converted to atom forces (3rd order)
    FILE * fh_dfile; // alamode dfile
    FILE * fh_ffile; // alamode ffile

    void write_xyz(); // function to write .xyz format position file
    void write_forces(); // function to write forces
    

  };
}

