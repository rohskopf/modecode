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

#ifdef COMPUTE_CLASS

ComputeStyle(mode/heatflux,ComputeModeHeatflux)

#else

#ifndef LMP_COMPUTE_MODE_HEATFLUX_H
#define LMP_COMPUTE_MODE_HEATFLUX_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeModeHeatflux : public Compute {
 public:
  ComputeModeHeatflux(class LAMMPS *, int, char **);
  ~ComputeModeHeatflux();
  void init();
  void compute_vector();

  int rank;
  int nprocs;
  
  int *nepp; // number of elements per proc (number of nmcc2s per proc).

  double **emat; // eigenvectors
  double **x0; // equilibrium positions
  double **u_p; // atomic displacements on a single proc
  double **u; // atomic displacements
  double *xm_p; // mode amplitudes on a single proc
  double *xm; // mode amplitudes
  double *fm_p; // mode velocities on a single proc
  double *vm; // mode velocities
  double *vm_p; // mode forces on a single proc
  double *fm; // mode forces
  double *fmi_p; // mode force due to atom i on a single proc
  double *fmi; // mode force due to atom i
               // index with fmi[natoms*n+i] where n=mode index
               //                                  i=atom index  
  double *freq; // frequencies
  double *tm_p; // mode temperatures on a single proc
  double *tm; // mode temperatures
  double *em; // mode energies

  //double *mcc3; // mcc3s

  int *sidemap; // side map, each index is an atom type (starting from zero) that yields:
                // 1 - left side of interface
                // 2 - right side of interface

  double *zvals; // z values that define the regions
  int *regions; // regions[i] gives the region that atom i belongs to

  double iflux_p; // interface flux from single proc
  double iflux; // interface flux from single proc
  double *miflux_p; // modal interface flux from single proc
  double *miflux; // modal interface flux from single proc

  double *iflux_arr; // 1D array that holds the total interface flux, for MPI reduce.
                     // iflux = iflux_arr[0]

  double qtot_p; // total harmonic heat flux at each timestep, on a single proc.
  double qtot; // total harmonic heat flux at each timestep.
  double **hf_p; // 2D array holds the time-averaged contributions of mode heat fluxes on a single proc.
  double **hf; // 2D array that holds the time-averaged contributions of mode heat fluxes. 
  double *qnm_p; // 1D array that holds time-averaged mode contributions to the heat flux, on a single proc.
  double *qnm; // 1D array that holds time-averaged mode contributions to the heat flux.
  double normalize; // integer argument of the "thermo" command
  double testicle; 
  int nsteps; // number of timesteps
  // Setting variable for output. 
  // All settings output total heat flux as a function of time in qtot.dat.
  int setting1; // 0 to compute <Qnm> using ../GV
                // 1 to compute Qnm(t) using GV in current directory, used if you wanna study behavior of a few particular modes. ONLY USE WITH ONE PROC FOR NOW!!!!!!!!!
                // 2 to compute Qnm(t) using ../GV. NOT CODED YET


  int nmcc2;
  struct mcc2_struct{
      int i,j;
      double val;
  };
  mcc2_struct *mcc2; // array of MCC2 values

  int nmcc3;
  struct mcc3_struct{
      int i,j,k;
      double val;
  };
  mcc3_struct *mcc3; // array of MCC3 values

  int nfc2;
  struct fc2_struct{
      int i,j;
      int a,b;
      double val;
  };
  fc2_struct *fc2; // array of FC2s

  int nfc3;
  struct fc3_struct{
      int i,j,k;
      int a,b,c;
      double val;
  };
  fc3_struct *fc3; // array of FC2s


 private:
  char *id_ke,*id_pe,*id_stress;
  class Compute *c_ke,*c_pe,*c_stress;

  FILE * fh_disp; // Atomic displacement file handle.
  FILE * fh_xm; // Mode coordinates file handle. 
  FILE * fh_fm; // Mode force file handle. 
  FILE * fh_tm; // Mode temperature file handle. 
  FILE * fh_miflux; // Modal interface flux file handle.
  FILE * fh_iflux; // Interface flux file handle.
  FILE * fh_fv; // F*V file handle for mode-mode interactions
  FILE * fh_fvtot; // F*V file handle for total F*V
  FILE * fh_em; // Mode energy file handle
  FILE * fh_qtep; // TEP heat flow
  FILE * fh_etep; // TEP total energy
  FILE * fh_qnm; // Qnm(t), used if setting1==1

  FILE * fh_debug; // debug file

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute heat/flux compute ID

Self-explanatory.

E: Compute heat/flux compute ID does not compute ke/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute pe/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute stress/atom

Self-explanatory.

*/
