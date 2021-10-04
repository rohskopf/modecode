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
    FILE * fh_disp;

    void initialize();
    void calcInitialState();
    void calcTimeDependence();
    void calcAmplitude(int, double);
    void calcVelocity(int, double);
    void convertMode2Cartesian();
    void readEmat();
    void readGV();

    // Input args
    double temperature; // Temperature in Kelvin
    int n_indx; // Index of mode "n" in Qnm
    double timestep; // timestep in ps
    double scale_factor; // factor to scale displacements by, for visualization purposes.
    
    double largest_gv; // largest gv within +/- 10 modes from n_indx

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
    double *xm0; // initial mode amplitudes
    double *vm0; // initial mode velocities
    
    double *mass;
    int *type;
    int *tag;
    double **x;

    double kb = 1.38064852e-23; // Boltzmann constant
    double pi = 3.1415926535;



  };
}

