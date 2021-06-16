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
  class Tic: protected Ptrs
  {
  public:
    Tic(class MC *);
    ~Tic();

    FILE * fh_debug;
    FILE * fh_fc2;
    FILE * fh_linewidths;

    void go();
    void readFcs();
    void readEmat();
    void readFrequencies();
    void readSMCC2();
    void calcLineWidths();
    double calcMCC3(int,int,int);
    double calcDistribution(int);
    //void calcTic(); 

    /* Booleans */
    bool pr_call;
    bool esp_call;

    /* Arrays */
    int *type; // LAMMPS
    int *tag; // LAMMPS
    double *mass; // LAMMPS
    double *linewidths_p; // linewidths for all 3N modes, on a single proc.
                          // modes not on this proc will have a linewidth of zero.
                          
    double *linewidths; // linewidths for all 3N modes.
    double *freq; // frequencies (THz)

    int rank;
    int natoms;
    int nmodes;
    int nfc2;
    int nfc3;
    int nfc4;
    int nsmcc2;

    double kb,hbar,pi;

    struct fc2_struct{
        int i,j;
        int a,b;
        double val;
    };
    fc2_struct *fc2;
    struct fc3_struct{
        int i,j,k;
        int a,b,c;
        double val;
    };
    fc3_struct *fc3;
    struct fc4_struct{
        int i,j,k,l;
        int a,b,c,d;
        double val;
    };
    fc4_struct *fc4;

    struct smcc2_struct{
        int n1,n2;
        double val;
    };
    smcc2_struct *smcc2;


    double **emat; // eigenvectors



  };
}

