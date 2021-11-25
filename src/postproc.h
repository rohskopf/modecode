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
    int task; // 1 - calculate power spectrum overlap.
    string ensemble_dirname;

    /* Booleans */

    /* Arrays */
    double **emat; // eigenvectors


    // Task 1 variables
    int ntimesteps;
    double sampling_interval; // time (in ps typically) between measured data points, or MD output.
    //int ncols; // number of columns, or modes, in the xm.dat and vm.dat files.
    int nind; // number of indices, or number of columns (not including the time column) in xm.dat and vm.dat files.
    int nens; // number of ensembles to average the PS overlap for.
    bool indices_bool; // true if INDICES file exists in current dir
    int *indices;
    int overlap_output_tag; // tag which uniquely names the overlaps.dat file.

    int natoms;
    int nmodes;
    int nfc2;
    int nfc3;
    int nfc4;



  };
}

