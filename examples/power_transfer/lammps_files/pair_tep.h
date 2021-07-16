/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tep,PairTep)

#else

#ifndef LMP_PAIR_TEP_H
#define LMP_PAIR_TEP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTep : public Pair {
 public:
  PairTep(class LAMMPS *);
  virtual ~PairTep();
  virtual void compute(int, int);
	
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
  void read_fc2();
  void read_fc3();
  void read_fc4();
  void read_equil();

  int nfc2;
  int nfc3;
  int nfc4;

 protected:
  int order;
  double **cut;
  double **d0,**alpha,**r0;
  double **morse1;
  double **offset;
  //double ****fc3;
  double **x0;

  double ****phi;
  double ***psi; // index with [3*i+a][3*j+b][3*k+c], has dimensions (3*natoms)^3

  double **emat; // eigenvectors

  double **u_p; // atomic displacements on a single proc
  double **u; // atomic displacements
  double *xm_p; // mode amplitudes on a single proc
  double *xm; // mode amplitudes
  double *fm_p; // mode velocities on a single proc
  double *vm; // mode velocities
  double *vm_p; // mode forces on a single proc
  double *fm; // mode forces

  //int nfc2;
  struct fc2_struct{
    int i,j;
    int a,b;
    double fc;
  };
  fc2_struct *fc2; // array of FC2 values

  struct fc3_struct{
    int i,j,k;
    int a,b,c;
    double fc;
  };
  fc3_struct *fc3; // array of FC2 values

  struct fc4_struct{
    int i,j,k,l;
    int a,b,c,d;
    double fc;
  };
  fc4_struct *fc4; // array of FC4 values

  FILE * fh_pe;
  FILE * fh_pe3;
  FILE * fh_ht;
  FILE * fh_ht2;
  FILE * fh_ht3;
  FILE * fh_debug;
  

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
