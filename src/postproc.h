#pragma once

#include <vector>
#include <string>
#include "mpi.h"

#include <iostream>
#include <new>
#include <cstdlib>
#include "ptrs.h"

using namespace std;

namespace MC_NS
{
  class Postproc: protected Ptrs
  {
  public:
    Postproc(class MC *);
    ~Postproc();

    // Constructor.
    FILE * fh_debug;
    int rank;

    // Functions.
    void initialize();
    void readEmat();
    void calc1();
    
    // Settings
    string ensemble_dirname;

    /* Booleans */

    /* Arrays */
    double **emat; // eigenvectors



    int natoms;
    int nfc2;
    int nfc3;
    int nfc4;



  };
}

