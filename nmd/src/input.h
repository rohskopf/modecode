/*
 input.h

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
  class Input: protected Pointers
  {
  public:
    Input(class EM3 *);
    ~Input();

    FILE * fh_debug; // Debug file handle

    // Declare member functions

    void readinput(); // function to read INPUT file
    void readconfig(); // function to read CONFIG and EQUIL files
    void initialize(); // function to initialize velocities and convert inputs
    void readParams(); // function to read IFCs and store in array

    // readinput() variables

    int nsteps; // number of timesteps 
    int space; // 0 for real space, 1 for mode space
    int nout; // output data this many timesteps
    int order; // Taylor expansion order
    double rc; // cutoff
    int neighcount;
    int newton;
    int offset; // whether or not to calculate energy offset
    double *mass; // mass
    double temp;
    double dt; // timestep
    double epsilon;
    double sigma;
   

    // readconfig() variables

    int natoms;
    int ntypes;
    double box[3];
    int *type;
    double **xa; // atomic positions
    double **xa0; // original atomic positions to store

    // initilize() variables

    double tau; // LJ time unit
    double velocity; // LJ velocity unit
    double force; // LJ force unit
    double pressure; // LJ pressure unit
    double temperature; // LJ temperature unit
    double **va; // velocities
    double volume; // box volume

    // readIfcs() variables

    //double ***fc2; // NxNx9 IFC2s
    //double ****fc3; // NxNxNx27 IFC3s

    double **emat; // eigenvectors
    double *freqs; // vector of frequencies

    int nmcc2;
    struct mcc2_struct{
        int i;
        double val;
    };
    mcc2_struct *mcc2; // array of MCC2 values

    int nmcc3;
    struct mcc3_struct{
        int i,j,k;
        double val;
    };
    mcc3_struct *mcc3; // array of MCC3 values

    int nmcc4;
    struct mcc4_struct{
        int i,j,k,l;
        double val;
    };
    mcc4_struct *mcc4; // array of MCC4 values

    double *xm; // mode amplitudes
    double *vm; // mode velocities

    int nfc2;
    struct fc2_struct{
        int i,j;
        int a,b;
        double val;
    };
    fc2_struct *fc2; // array of FC2 values
    
    int nfc3;
    struct fc3_struct{
        int i,j,k;
        int a,b,c;
        double val;
    };
    fc3_struct *fc3; // array of FC2 values

    int nfc4;
    struct fc4_struct{
        int i,j,k,l;
        int a,b,c,d;
        double val;
    };
    fc4_struct *fc4; // array of FC2 values

  };
}

