/*
 compute.h

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
  class Compute: protected Pointers
  {
  public:
    Compute(class EM3 *);
    ~Compute();

    FILE * fh_debug; // positions output in .xyz format

    void compute_ke(); //  computes kinetic energy, total energy, and temperature
    void compute_pressure(); // computes stress virial and total pressure
    void compute_msd(); // compute mean squared distance
    void compute_com(); // computes center of mass quantities (velocity, position, etc.)
    void compute_rmsd(); // computes RMSD of quantities
    
    double ke; // kinetic energy
    double etot; // total energy
    double temp; // temperature
    double pressure; // total pressure
    double msd; // mean squared distance
    double rmsd_etot; // RMSD of total energy

    double *em; // mode energies

  };
}

