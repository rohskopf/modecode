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
  class In: protected Ptrs
  {
  public:
    In(class MC *);
    ~In();

    void readInput();
    void readConfig();
    void evaluate();
    void calcFC2();
    double calcD1(int,int,int,int);
    double calcD2(int,int,int,int,int,int);
    double calcD3(int,int,int,int,int,int,int,int);
    double calcF(int*,int*,double*,int);
    double calcD3_MTEP(int,int,int);
    double calcE_MTEP(int*,double*,int);
    double calcRij(double,double,double,double,double,double);
    
    bool in_call = true; // tells whether in class was called.

    FILE * fh_data;
    FILE * fh_debug;
    FILE * fh_fc2;
    FILE * fh_fc3;
    FILE * fh_fc4;

    // User input variables
    int order;

    /* readinput() variables */
    int ndisps;
    double delta;
    int pop_size_in;
    double mut_rate_in;
    double eliteperc_in; 
    int unknowns_in;
    int unknowns_extra_in;
    int M_in;
    int N_rad;
    int M;
    double w_f_in;
    double N_ang;
    int N_types;
    double rc;
    double rc2;
    double rc3;
    double tol2;
    double tol3;
    double tol4;
    double w_p_in;
    double w_p_s_in;
    // Declare vector inputs
    vector<double> intervals_in;
    vector<int> ssl_in;
    vector<double> d_vec_in;
    vector<double> sym_vec_in;
    vector<int> N_vec;
    // Declare string inputs
    string sym_coeffs_in;
    string commands; // lammps commands
    string nontab_commands; // non-tabulated parameter commands

    int N, Ntypes, Nbonds, Nangles, Nbondtypes, Nangletypes;
    char* types;
    double** positions;
    double* box;

    /* LAMMPS variables */
    int inum;
    int *ilist,*jlist,*numneigh,**firstneigh, *jneighlist, *type, *tag;
    //double **x;

    double **emat; // eigenvectors
    double **x0; // equilibrium positions
    double *freqs; // vector of frequencies

    int rank; // proc rank

    struct fc2_struct{
        int i,j; // atom indices
        int a,b;
        double fc;
    };
    fc2_struct *fc2;
    fc2_struct *fc2_all;

    struct fc3_struct{
        int i,j,k;
        int a,b,c;
        double fc;
    };
    fc3_struct *fc3;
    struct fc4_struct{
        int i,j,k,l;
        int a,b,c,d;
        double fc;
    };
    fc4_struct *fc4;


  };
}

