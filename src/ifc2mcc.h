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
  class Ifc2mcc: protected Ptrs
  {
  public:
    Ifc2mcc(class MC *);
    ~Ifc2mcc();

    FILE * fh_debug;

    //double tol; // tolerance

    void go(double);
    void readFcs();
    void readEmat();
    unsigned long long int nChoosek(unsigned long long int,unsigned long long int);
    void extract();
    void average(int,int);
    void extractFew(int,int);
    void calcFew();
    void go_spatial();
    void go_gv(int); // generalized velocities
    void go_n1(int);
    void extract_smcc2(int);

    int task; // 0 - calculate 
              // 1 - extract
    int rank;
    int natoms;
    int order;
    int nfc2;
    int nfc3;
    int nfc4;
    int lowmode;
    int highmode;

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
    struct mcc3_struct{
        int i,j,k;
        double mcc;
    };
    mcc3_struct *mcc3;

    double **emat; // eigenvectors



  };
}

