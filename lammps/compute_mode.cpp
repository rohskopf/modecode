/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Andrew Rohskopf (MIT)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Compute mode-level quantities such as contributions to PE, KE, etc. 
   for a group of atoms. 
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include "compute_mode.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "mpi.h"
#include "universe.h"

#include <sstream>
#include <vector>
#include <fstream>
#include <string>

using namespace LAMMPS_NS;
using namespace std;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

ComputeMode::ComputeMode(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute mode command");

  scalar_flag = 1;
  extscalar = 1;

  // Get output setting.

  out_setting = atoi(arg[3]);
  if (out_setting > 3){
    error->all(FLERR,"Illegal compute mode output setting > 3.");
  }
  //printf("%c\n", out_setting[0]);

  //rank = MPI_COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  fh_disp = fopen("disp.dat", "w");
  fh_xm = fopen("xm.dat", "w");
  fh_vm = fopen("vm.dat", "w");
  fh_fm = fopen("fm.dat", "w");
  fh_tm = fopen("tm.dat", "w");
  fh_em = fopen("em.dat", "w");
  fh_em2 = fopen("em2.dat", "w");
  fh_mc = fopen("mc.dat", "w");
  fh_jnm = fopen("jnm.dat", "w");

}

/* ---------------------------------------------------------------------- */

ComputeMode::~ComputeMode()
{


  memory->destroy(emat);
  memory->destroy(x0);
  memory->destroy(u_p);
  memory->destroy(u);
  memory->destroy(xm_p);
  memory->destroy(xm);
  memory->destroy(freq);
  memory->destroy(vm_p);
  memory->destroy(vm);
  memory->destroy(fm_p);
  memory->destroy(fm);
  memory->destroy(fmi_p);
  memory->destroy(fmi);
  memory->destroy(miflux_p);
  memory->destroy(miflux);
  memory->destroy(tm);
  memory->destroy(em);
  if (indices_bool){
    memory->destroy(indices);
  }
  if (pairs_bool){
    memory->destroy(indices);
  }
  //memory->destroy(mcc3);

  //memory->destroy(sidemap);
  //memory->destroy(regions);
  //memory->destroy(zvals);

  fclose(fh_disp);
  fclose(fh_xm);
  fclose(fh_vm);
  fclose(fh_fm);
  fclose(fh_tm);
  fclose(fh_em); 
  fclose(fh_em2); 
  fclose(fh_mc);
  fclose(fh_jnm);

}

/* ---------------------------------------------------------------------- */

void ComputeMode::init()
{


  int natoms = atom->natoms;
  int nlocal = atom->nlocal;

  /*
  if (universe->me == 0) {
    printf("This is proc 0.\n", universe->me);
  }
  */

  //printf("    -------------- Initializing Mode Computations on Proc %d. -------------------\n", universe->me);
  printf("Reading mode inputs for %d atoms on proc %d.\n", nlocal, universe->me);

  memory->create(emat,natoms*3,natoms*3,"mode:emat");
  memory->create(x0,natoms,3,"mode:x0");
  memory->create(u_p,natoms,3,"mode:u_p");
  memory->create(u,natoms,3,"mode:u");
  memory->create(xm_p,natoms*3,"mode:xm_p");
  memory->create(xm,natoms*3,"mode:xm");
  for (int n=0; n<3*natoms; n++){
    xm[n]=0.0;
    xm_p[n]=0.0;
  }
  memory->create(freq,natoms*3,"mode:freq");
  memory->create(vm_p,natoms*3,"mode:vm_p");
  memory->create(vm,natoms*3,"mode:vm");
  memory->create(fm_p,natoms*3,"mode:fm_p");
  memory->create(fm,natoms*3,"mode:fm");
  memory->create(fmi_p,natoms*natoms*3,"mode:fmi_p");
  memory->create(fmi,natoms*natoms*3,"mode:fmi");

  memory->create(miflux_p,natoms*3,"mode:miflux_p");
  memory->create(miflux,natoms*3,"mode:miflux");
  memory->create(tm,natoms*3,"mode:tm");
  memory->create(em,natoms*3,"mode:em");
  for (int n=0; n<3*natoms; n++){
    tm[n]=0.0;
    em[n]=0.0;
  }

  /* EQUIL */

  ifstream readfile;

  for (int i = 0; i < natoms; i++) {
      for (int a = 0; a < 3; a++) {
          x0[i][a] = 0.0;
      }
  }

  readfile.open("EQUIL");

  if (!readfile.is_open()) {
      printf("Unable to open EQUIL.\n");
      exit(1);
  }

  //printf("natoms: %d\n",  natoms);
  for (int i=0;i<natoms;i++){
      for (int a=0;a<3;a++){
          readfile>>x0[i][a];
      }
  }

  /* EMAT */

  ifstream readfile2;

  for (int i = 0; i < natoms*3; i++) {
      for (int j = 0; j < natoms*3; j++) {
          emat[i][j] = 0.0;
      }
  }

  readfile2.open("../EMAT");
  //readfile2.open("ev_real.txt");

  if (!readfile2.is_open()) {
      printf("Unable to open EMAT.\n");
      exit(1);
  }

  //printf("natoms: %d\n",  natoms);
  for (int i=0;i<3*natoms;i++){
      for (int j=0;j<3*natoms;j++){
          readfile2>>emat[i][j];
      }
  }

  /* FREQUENCIES */

  ifstream readfile3;

  for (int i = 0; i < natoms*3; i++) {
    freq[i]=0.0;
  }

  readfile3.open("FREQUENCIES");
  //readfile2.open("ev_real.txt");

  if (!readfile3.is_open()) {
      printf("Unable to open FREQUENCIES.\n");
      exit(1);
  }

  //printf("natoms: %d\n",  natoms);
  for (int i=0;i<3*natoms;i++){
    readfile3>>freq[i];
  }

  /* INDICES */
  ifstream fh_ind;
  fh_ind.open("INDICES");

  if (fh_ind.is_open()) {
      //printf("Unable to open INDICES file, won't output any mode quantities.\n");
    indices_bool = true;
    string line;
    nind = 0;
    while (getline(fh_ind, line))
    {
        nind=nind+1;
    }

    printf("  Found %d Mode Indices.\n", nind);

    memory->create(indices,nind,"mode:indices");

    fh_ind.close();

    ifstream fh_ind;
    fh_ind.open("INDICES");

    for (int n=0; n<nind; n++){
      fh_ind >> indices[n];
    }

    fh_ind.close();
   
  }
  else {
    indices_bool = false;
    printf("Unable to open INDICES file, won't output any mode quantities.\n");
    nind = 0.0;
  }

  /* PAIRS, for mode flux calculation */
  ifstream fh_pair;
  fh_pair.open("PAIRS");

  if (fh_pair.is_open()) {
    pairs_bool = true;
    string line;
    npairs= 0;
    while (getline(fh_pair, line))
    {
        npairs=npairs+1;
    }

    printf("  Found %d mode pairs for mode flux calculation.\n", npairs);

    memory->create(pairs,npairs,2,"mode:pairs");

    fh_pair.close();

    ifstream fh_pair;
    fh_pair.open("PAIRS");

    for (int n=0; n<npairs; n++){
      fh_pair >> pairs[n][0] >> pairs[n][1];
    }

    fh_pair.close();
  }
  else {
    pairs_bool = false;
    printf("Unable to open PAIRS file, won't output any mode flux quantities.\n");
    npairs = 0.0;
  }

  /*
  ifstream fh_ind("INDICES");
  string line;
  nind = 0;
  while (getline(fh_ind, line))
  {
      nind=nind+1;
  }

  printf("  Found %d Mode Indices.\n", nind);

  memory->create(indices,nind,"mode:indices");

  fh_ind.close();

  ifstream readfile6;
  readfile6.open("INDICES");

  if (!readfile6.is_open()) {
      printf("Unable to open INDICES file, won't output any mode quantities.\n");
  }

  for (int n=0; n<nind; n++){
    readfile6 >> indices[n];
  }

  readfile6.close();
  */

  /* SIDEMAP */
  /*
  ifstream readfile4;

  readfile4.open("SIDEMAP");
  //readfile2.open("ev_real.txt");

  if (!readfile4.is_open()) {
      printf("Unable to open SIDEMAP.\n");
      exit(1);
  }

  int ntypes = atom->ntypes;
  printf("    %d atom types.\n", ntypes);
  memory->create(sidemap,ntypes,"mode:sidemap");

  for (int t=0; t<ntypes; t++){
    readfile4>>sidemap[t];
    //printf(" %d\n", sidemap[t]);

  }
  */
  /* REGIONS */
  /*
  ifstream readfile5;

  readfile5.open("REGIONS");

  if (!readfile5.is_open()) {
      printf("Unable to open REGIONS.\n");
      exit(1);
  }
  memory->create(regions,natoms,"mode:regions");
  memory->create(zvals,3,"mode:zvals");
  for (int i=0; i<3; i++){
    readfile5>>zvals[i];
  }
  // assign regions to each atom based on equilibrium position
  for (int i=0; i<natoms; i++){
    if (x0[i][2]<zvals[0]) regions[i]=1;
    else if (zvals[0] <= x0[i][2] && x0[i][2] < zvals[1]) regions[i]=2;
    else if (zvals[1] <= x0[i][2] && x0[i][2] < zvals[2]) regions[i]=3;
    else if (x0[i][2] >= zvals[3]) regions[i]=4;
    else printf("EQUILIBRIUM POSITIONS OR REGIONS ARE MESSED UP.\n");

  }
  */
  /* MCC3 */
  /*
  ifstream fh_mcc3("MCC3");
  string line;
  nmcc3 = 0;
  while (getline(fh_mcc3, line))
  {
      nmcc3=nmcc3+1;
  }

  printf("  Found %d MCC3s.\n", nmcc3);

  memory->create(mcc3,nmcc3,"mode:mcc3");

  fh_mcc3.close();

  ifstream readfile6;
  readfile6.open("MCC3");

  for (int n=0; n<nmcc3; n++){
    readfile6 >> mcc3[n].i >> mcc3[n].j >> mcc3[n].k >> mcc3[n].val;
  }

  readfile6.close();

  //printf(" ****************** %d %d %d %e\n", mcc3[5].i, mcc3[5].j, mcc3[5].k, mcc3[5].val);

  /*
  readfile6.open("MCC3");

  if (!readfile6.is_open()) {
      printf("Unable to open MCC3.\n");
      exit(1);
  }
  int nmcc3;
  //readfile6>>nmcc3;
  //printf(" %d MCC3s.\n",nmcc3);
  while (getline(fh_input, line))
  {
      ncommands++;
      nleft = findNumberOf(line,"[");
      nright = findNumberOf(line,"]");
      //printf("nleft nright: %d %d\n", nleft,nright);
      if (nleft==nright) nparams = nparams + nleft;
      else popserror->exit("popsinput.cpp", "Number of left brackets not equal number of right brackets in LAMMPS input script.");

  }
  
  memory->create(mcc3,nmcc3,"mode:mcc3");
  int i,j,k;
  for (int q=0; q<nmcc3; q++){
    readfile6>>i>>j>>k>>mcc3[q];
  }
  */

  /*
  for (int q=0; q<nmcc3; q++){
    printf("%e\n", mcc3[q]);
  }
  */
  
}

/* ---------------------------------------------------------------------- */

double ComputeMode::compute_scalar()
{

  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;

  // Compute atomic displacements.
  //fprintf(fh_disp, "\n");-
  //fprintf(fh_disp, "Timestep: %d\n", update->ntimestep);
  //fprintf(fh_disp, "\n");
  int iindx;
  for (int i=0; i<nlocal; i++){
    iindx = tag[i]-1;
    //printf("%d %d %f %f %f\n", i, iindx, x[i][0], x[i][1],x[i][2]);
    if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
      //u_p[i]=0.0;
      //u[i] = 0.0; // zero the total coordinates
      for (int a=0; a<3; a++){
        u_p[i][a] = (x[i][a]-x0[tag[i]-1][a]);
        if (std::abs(u_p[i][a]) > domain->boxhi[a]/2.0 && u_p[i][a] > 0) u_p[i][a] -= domain->boxhi[a];
        else if (std::abs(u_p[i][a]) > domain->boxhi[a]/2.0 && u_p[i][a] < 0) u_p[i][a] += domain->boxhi[a];

        //u_p[i][a] = u_p[i][a]*1e-10;

        //fprintf(fh_disp,"%e ", x[i][a]-x0[i][a]);
      }
    }
  }
  // Check displacements.
  /*
  for (int i=0; i<nlocal; i++){
    printf("%f %f %f\n", u_p[i][0], u_p[i][1],u_p[i][2]);

  }
  */
  // NEED TO MAKE DISPLACEMENT VECTOR 1D if you want a total displacement vector for the system.
  //MPI_Allreduce(u_p,u,3*natoms,MPI_DOUBLE,MPI_SUM,world);

  //printf(" natoms: %d\n", natoms);
  // Mode coordinates, velocities, and forces. 
  for (int n=0; n<3*natoms; n++){
  //for (int n=7; n<8; n++){
    //printf("n: %d\n", n);
    xm_p[n]=0.0;
    xm[n] = 0.0; // zero the total coordinates
    vm_p[n]=0.0;
    vm[n] = 0.0; // zero the total velocities
    fm_p[n]=0.0;
    fm[n] = 0.0; // zero the total forces
    for (int i=0; i<nlocal; i++){
      if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
        for (int a=0; a<3; a++){
          //u[i][a] = u[i][a]*1e-10;
          double ms = mass[type[i]]; //mass[type[i]]*(1.0/6.0221409e+23)*1e-3;
          double masskg = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
          //xm_p[n] += sqrt(mass[type[i]])*emat[3*i+a][n]*u[i][a]*1e-10*(1.0/6.0221409e+23)*1e-3; // Convert to SI units.;
          //xm_p[n] += sqrt(ms)*emat[3*i+a][n]*u_p[i][a]; // sqrt(kg/kmol) * Angstrom
          xm_p[n] += sqrt(masskg)*emat[3*(tag[i]-1)+a][n]*u_p[i][a]; // sqrt(kg) * Angstrom
          fm_p[n] += sqrt(ms)*emat[3*(tag[i]-1)+a][n]*f[i][a]; // sqrt(kg/kmol) * eV/A
          //fmi_p[natoms*n+i] = sqrt(ms)*emat[3*i+a][n]*f[i][a];
          //vm_p[n] += sqrt(masskg)*emat[3*i+a][n]*v[i][a]*100; // sqrt(kg) * m/s, since 1 A/ps = 100 m/s.
          vm_p[n] += sqrt(masskg)*emat[3*(tag[i]-1)+a][n]*v[i][a]; // sqrt(kg) * A/ps
          
        }
      }
    }
  }

  MPI_Allreduce(xm_p,xm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(vm_p,vm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(fm_p,fm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  //MPI_Allreduce(fmi_p,fmi,3*natoms*natoms,MPI_DOUBLE,MPI_SUM,world);

  // Check that the basis represnts the displacements accurately.
  /*
  double ui;
  for (int i=0; i<nlocal; i++){
    double mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
    ui = 0.0;
    for (int n=0; n<3*natoms; n++){
      ui += (1.0/sqrt(mi))*emat[3*(tag[i]-1)+0][n]*xm[n];

    }
    fprintf(fh_disp, "%d %e\n", tag[i],u_p[i][0]);
  }
  */



  // Mode temperatures and energies
  double pi = 3.1415926535;
  double pem = 0.0; 
  for (int n=3; n<3*natoms; n++){
    //temp_mode[n] = ( xm[n]*xm[n]*input->mcc2[n].val )/(1.38064852e-23);
    //tm[n] = 0.5*(vm[n]*vm[n])/(1.38064852e-23);
    tm[n] = 0.5*(vm[n]*vm[n]*100*100)/(1.38064852e-23); // factors of 100 since 100 m/s = 1 A/ps
    em[n] = (0.5*vm[n]*vm[n]*100*100 + 0.5*freq[n]*freq[n]*2.0*3.1415926535*2.0*3.1415926535*1e12*1e12*xm[n]*xm[n]*1e-10*1e-10)*6.242e+18; // eV, since 6.242e+18 eV = 1 J
    pem += 0.5*freq[n]*freq[n]*2.0*3.1415926535*2.0*3.1415926535*1e12*1e12*xm[n]*xm[n]*1e-10*1e-10*6.242e+18; // eV, since 6.242e+18 eV = 1 J
    //if (n==2209){
    //    printf(" vm[n]: %e\n", vm[n]);
    //    printf(" temp_mode[n]: %e\n", temp_mode[n]);
    //}
    //kemi[n] = 0.5*vm[n]*vm[n];
    //emi[n] = (kemi[n] + pemi[n])*6.242e+18;
  }

  scalar = pem;


  // Print mode quantities.
  // Convert units to sqrt(M)*A, where M = g/mol or kg/kmol for coordinates and velocities.
  if (universe->me==0){
    //printf("%d %d %d\n", nind, indices[0], indices[1]);
    if (out_setting==0){ // Print .dat files for plotting.
      fprintf(fh_xm, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_vm, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_em, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_tm, "%f", 0.5*update->ntimestep*1e-3);
      for (int i=0; i<nind; i++){
        fprintf(fh_xm, " %e", xm[indices[i]]*sqrt(6.02214076e23*1e3));
        fprintf(fh_vm, " %e", vm[indices[i]]*sqrt(6.02214076e23*1e3));
        fprintf(fh_em, " %e", em[indices[i]]);
        fprintf(fh_tm, " %e", tm[indices[i]]);
        //fprintf(fh_mc, "%e %e %e\n", freq[indices[i]], xm[indices[i]]*sqrt(6.02214076e23*1e3),vm[indices[i]]*sqrt(6.02214076e23*1e3));      
      }
      fprintf(fh_xm, "\n");
      fprintf(fh_vm, "\n");
      fprintf(fh_em, "\n");
      fprintf(fh_tm, "\n");
    }
    else if (out_setting == 1){ // print .dat files for gifs.
      fprintf(fh_em2, "\n");
      fprintf(fh_em2, "\n");  
      
      for (int n=0; n<3*natoms; n++){
        fprintf(fh_em2, "%d %e %e\n", n,freq[n],em[n]);
      }
    }
    else if (out_setting == 2){ // print .dat files for plots and gifs.

      // For plots.
      fprintf(fh_xm, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_vm, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_em, "%f", 0.5*update->ntimestep*1e-3);
      fprintf(fh_tm, "%f", 0.5*update->ntimestep*1e-3);
      for (int i=0; i<nind; i++){
        fprintf(fh_xm, " %e", xm[indices[i]]*sqrt(6.02214076e23*1e3));
        fprintf(fh_vm, " %e", vm[indices[i]]*sqrt(6.02214076e23*1e3));
        fprintf(fh_em, " %e", em[indices[i]]);
        fprintf(fh_tm, " %e", tm[indices[i]]);
        //fprintf(fh_mc, "%e %e %e\n", freq[indices[i]], xm[indices[i]]*sqrt(6.02214076e23*1e3),vm[indices[i]]*sqrt(6.02214076e23*1e3));      
      }
      fprintf(fh_xm, "\n");
      fprintf(fh_vm, "\n");
      fprintf(fh_em, "\n");
      fprintf(fh_tm, "\n");

      // For gifs.
      fprintf(fh_em2, "\n");
      fprintf(fh_em2, "\n");  
      for (int n=0; n<3*natoms; n++){
        fprintf(fh_em2, "%d %e %e\n", n,freq[n],em[n]);
      }

    }

    else if (out_setting == 3){ // calculate and print mode fluxes in jnm.dat
      // Convert units to sqrt(M)*A, where M = g/mol or kg/kmol for coordinates and velocities.

      //fprintf(fh_jnm, "%f ", 0.5*update->ntimestep*1e-3);
      for (int i=0; i<npairs; i++){
        //fprintf(fh_jnm, "%e ", xm[pairs[i][0]]*vm[pairs[i][1]]*sqrt(6.02214076e23*1e3)*sqrt(6.02214076e23*1e3));
        // For now, make the first column zero just so we can use Kia's code, and be sure to use only one pair.
        fprintf(fh_jnm, "0 %e", xm[pairs[i][0]]*vm[pairs[i][1]]*sqrt(6.02214076e23*1e3)*sqrt(6.02214076e23*1e3));
      }
      fprintf(fh_jnm, "\n");

    }

  }


  return scalar;
}
