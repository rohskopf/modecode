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
   Compute mode F*V using MCC3s. 
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include "compute_mode_fv.h"
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

ComputeModeFv::ComputeModeFv(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute mode command");

  scalar_flag = 1;
  extscalar = 1;

  // Get output setting.
  /*
  out_setting = atoi(arg[3]);
  if (out_setting > 2){
    error->all(FLERR,"Illegal compute mode output setting > 2.");
  }
  */
  //printf("%c\n", out_setting[0]);

  //rank = MPI_COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  fh_fv = fopen("fv.dat", "w");

}

/* ---------------------------------------------------------------------- */

ComputeModeFv::~ComputeModeFv()
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
  memory->destroy(mcc3);

  fclose(fh_fv);

}

/* ---------------------------------------------------------------------- */

void ComputeModeFv::init()
{


  int natoms = atom->natoms;

  printf("    -------------- Initializing Mode F*V Computations on Proc %d. -------------------\n", universe->me);
  printf("    Reading mode inputs for %d atoms.\n", natoms);

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
      printf("Unable to open EMAT.\n");
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
      //printf("ASDF\n");
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

  /* MCC3 */

  ifstream fh_mcc3("../MCC3");
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
  readfile6.open("../MCC3");

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

double ComputeModeFv::compute_scalar()
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

  // Anharmonic forces and power transfer.
  if (universe->me==0){
    double fv;
    int i,j,k;
    fprintf(fh_fv, " %0.2f ", 0.5*update->ntimestep*1e-3);
    //printf(" %d\n", nmcc3);
    int counter = 0;
    int n_indx, m_indx, l_indx;
    int n,m,l;
    bool thismode;
    //for (int n=0; n<nmcc3; n++){
    for (int s=0; s<nmcc3; s++){

      // Calculate power transfer
      /*
      i=mcc3[n].i;
      j=mcc3[n].j;
      k=mcc3[n].k;
      //printf("%e\n", mcc3[n].val);
      //if (abs(mcc3[n].val)>1e48){
      //printf("ASDF\n");
      fv = -1.0*mcc3[s].val*xm[j]*xm[k]*1e-10*1e-10*vm[i]*100*6.242e+6; // eV/ps
      */
      
      // This code block is the recent working example as of 07-14-21
      /*
      n=mcc3[s].i;
      m=mcc3[s].j;
      l=mcc3[s].k;
      //fv = -(1.0/2.0)*mcc3[s].val*xm[m]*xm[l]*1e-10*1e-10*vm[n]*100*6.242e+6; // eV/ps
      fv = -(1.0/1.0)*mcc3[s].val*xm[m]*xm[l]*1e-10*1e-10*vm[n]*100*6.242e+6; // eV/ps
      */

      
      // This code block uses the entire MCC3 matrix
      // Find index of "n" (the mode under consideration)
      // Also assign m and l indices accordingly.
      thismode=false;
      if (mcc3[s].i==10){ 
        n_indx=0;
        m_indx=1;
        l_indx=2;
        n=mcc3[s].i;
        m=mcc3[s].j;
        l=mcc3[s].k;
        thismode=true;
      }
      else if (mcc3[s].j==10){ 
        n_indx=1;
        m_indx=0;
        l_indx=2;
        n=mcc3[s].j;
        m=mcc3[s].i;
        l=mcc3[s].k;
        thismode=true;
      }
      else if (mcc3[s].k==10){ 
        n_indx=2;
        m_indx=0;
        l_indx=1;
        n=mcc3[s].k;
        m=mcc3[s].i;
        l=mcc3[s].j;
        thismode=true;
      }
      //printf("%d %d %d\n", n,m,l);
      
      //if (mcc3[s].i==10){
      //  fv = -(1.0/1.0)*mcc3[s].val*xm[m]*xm[l]*1e-10*1e-10*vm[n]*100*6.242e+6; // eV/ps
      //}

      if (thismode){

        // Playing around with other expressions
        //fv = -(1.0/3.0)*mcc3[s].val*(  xm[m]*xm[l]*vm[n] - xm[n]*xm[l]*vm[m])*1e-10*1e-10*100*6.242e+6; // eV/ps
        fv = -(1.0/6.0)*mcc3[s].val*( xm[m]*xm[l]*vm[n] - xm[n]*xm[l]*vm[m])*1e-10*1e-10*100*6.242e+6; // eV/ps
        //fv = -(1.0/2.0)*mcc3[s].val*( xm[m]*xm[l]*vm[n])*1e-10*1e-10*100*6.242e+6; // eV/ps
        // print the power transfers.
        fprintf(fh_fv, "%e ", fv);

      }
      else{ 
        fv=0.0;
      }


      // print the power transfers.
      //fprintf(fh_fv, "%e ", fv);
      counter++;
      //}

    }
    fprintf(fh_fv, "\n");

    if (update->ntimestep==0){
      printf(" ***** Found %d MCC3s.\n", counter);
    }
  }

  return 0;
}
