/*
 update.h

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
  class Update: protected Pointers
  {
  public:
    Update(class EM3 *);
    ~Update();

    FILE * fh_debug; // Debug file handle

    void initialize();
    void integrate();
    void integrate2(); // Stormer-Verlet integration

    int natoms;
    //double dt;
    //double *xa; // atomic positions
    //double *va; // atomic velocities
    double **xa_old; // atomic positions at previous step
    double **aa; // atomic accelerations

    double *xm_old; // mode positions at previous step
    double *xm; // mode positions
    double *am; // mode accelerations
    double *vm; // mode velocities
    

  };
}

