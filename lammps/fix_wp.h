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

#ifdef FIX_CLASS

FixStyle(wp,FixWp)

#else

#ifndef LMP_FIX_WP_H
#define LMP_FIX_WP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWp : public Fix {
 public:
  FixWp(class LAMMPS *, int, char **);
  virtual ~FixWp();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  virtual void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_vector(int);

  double memory_usage();

  //double amp;
  int modeindx; // index of mode (starting from 0) to pump
  //int natoms; // number of atoms in box
  double time; // time of simulation in seconds
  double **emat; // eigenvectors
  double *freq; // frequencies

  double *tm; // mode temperatures
  double **x0; // equilibrium positions
  double **up; // packet displacements
  double **vp; // packet velocities

  double *xm; // mode amplitudes
  double *xm_p; // mode amplitudes on a single proc
  double *vm; // mode velocities
  double *vm_p; // mode forces on a single proc

  double *em; // mode energy
  double **u_p; // atomic displacements on a single proc


  FILE * fh_d; // debug

  double amp,k0,xi,vg,z0,a0; // WP parameters
  int direction; // WP direction (0,1, or 2 for Cartesian directions)

 protected:
  double xvalue,yvalue,zvalue;
  int varflag,iregion;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  double foriginal[3],foriginal_all[3];
  int force_flag;
  int nlevels_respa,ilevel_respa;

  int maxatom;
  double **sforce;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
