/*
 potential.h

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
  class Potential: protected Pointers
  {
  public:
    Potential(class EM3 *);
    ~Potential();

    FILE * fh_debug; // Debug file handle

    void calculate();

    double ***mcc3; // copy of MCC3 from input

    int natoms;
    double pea; // real space potential energy (atoms)
    double pem; // mode space potential energy
    double *em; // mode energies
    double *f; // Mode forces
    double p_xx; // xx stress tensor component
    double p_yy; // yy stress tensor component
    double p_zz; // zz stress tensor component
    

  };
}

