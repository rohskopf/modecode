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
  class Verify: protected Ptrs
  {
  public:
    Verify(class MC *);
    ~Verify();

    FILE * fh_debug;

    void readMCC(const char*);
    void compareRandom(int,double,int);
    double calcCohesiveEnergy();
    double calcEnergyLmp(double **);
    double calcEnergyMtep(double *);
    
    double ***mcc3; // array of MCC3 values
    int nmodes; // number of modes

  };
}

