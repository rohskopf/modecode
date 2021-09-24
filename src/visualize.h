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
  class Visualize: protected Ptrs
  {
  public:
    Visualize(class MC *);
    ~Visualize();

    FILE * fh_debug;

    void initialize();
    void calcInitialState();
    void calcTimeDependence();
    void readEmat();
    void readGV();

    double temperature; // Temperature in Kelvin
    int n_indx; // Index of mode "n" in Qnm

    int natoms;
    int ngv;

    double *freq; // frequencies
    double **emat; // eigenvectors
    struct gv_struct{
        int n1,n2;
        double val;
    };
    gv_struct *gv;


    double *xm; // mode amplitudes
    double *vm; // mode velocities

    double kb = 1.38064852e-23; // Boltzmann constant



  };
}

