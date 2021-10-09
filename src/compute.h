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
  class Compute: protected Ptrs
  {
  public:
    Compute(class MC *);
    ~Compute();

    FILE * fh_debug;
    FILE * fh_fc2;

    void participationRatio();
    void eigenvectorSpatialParameter();
    void readFcs();
    void readEmat();
    void sumDMCC(); // Calculate Kn = sum(Knm) for all modes n.

    /* Booleans */
    bool pr_call;
    bool esp_call;

    /* Arrays */
    double *prs; // participation ratios
    double *esp; // eigenvector spatial parameters
    double *kn; // Kn = sum(Knm), sum of directional mode coupling constants.

    int order;  
    int natoms;
    int nfc2;
    int nfc3;
    int nfc4;

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


    double **emat; // eigenvectors



  };
}

