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
    void task1();
    double integrate(double*,double*,int);
    
    // Settings
    string ensemble_dirname;

    /* Booleans */

    /* Arrays */
    double **emat; // eigenvectors


    // Task 1 variables
    int ntimesteps;
    double sampling_interval; // time (in ps typically) between measured data points, or MD output.

    int natoms;
    int nmodes;
    int nfc2;
    int nfc3;
    int nfc4;



  };
}

