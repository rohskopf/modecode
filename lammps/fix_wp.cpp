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

#include <cstring>
#include <cstdlib>
#include "fix_wp.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "comm.h"

#include <sstream>
#include <vector>
#include <fstream>
#include <string>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixWp::FixWp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), idregion(NULL), sforce(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix wp command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = nlevels_respa = 0;

  // optional args

  iregion = -1;
  idregion = NULL;

  maxatom = 1;
  memory->create(sforce,maxatom,3,"setforce:sforce");

  amp = atof(arg[3]);
  k0 = atof(arg[4]);
  xi = atof(arg[5]);
  vg = atof(arg[6]);
  z0 = atof(arg[7]);
  a0 = atof(arg[8]);


  char debug[64];
  sprintf (debug, "D_FIX_WP_%d", comm->me);

  fh_d = fopen(debug, "w");


}

/* ---------------------------------------------------------------------- */

FixWp::~FixWp()
{
  if (copymode) return;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);


  memory->destroy(x0);
  memory->destroy(emat);
  memory->destroy(freq);
  memory->destroy(tm);
  memory->destroy(vm);
  memory->destroy(vm_p);
  memory->destroy(xm);
  memory->destroy(xm_p);
  memory->destroy(up);
  memory->destroy(vp);
  memory->destroy(em);
  memory->destroy(u_p);


  fclose(fh_d);

}

/* ---------------------------------------------------------------------- */

int FixWp::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWp::init()
{

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix setforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,nlevels_respa-1);
    else ilevel_respa = nlevels_respa-1;
  }

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
    if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
    if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
    if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  }
  if (flag)
    error->all(FLERR,"Cannot use non-zero forces in an energy minimization");

  int natoms = atom->natoms;
  int nlocal = atom->nlocal;

  memory->create(emat,natoms*3,natoms*3,"wp:emat");
  memory->create(freq,natoms*3,"wp:freq");
  memory->create(x0,natoms,3,"wp:x0");
  memory->create(up,nlocal,3,"wp:up");
  for (int i=0; i<nlocal; i++){
    for (int a=0; a<3; a++) up[i][a]=0.0;
  }
  memory->create(vp,nlocal,3,"wp:vp");
  for (int i=0; i<nlocal; i++){
    for (int a=0; a<3; a++) vp[i][a]=0.0;
  }

  memory->create(xm_p,natoms*3,"mode:xm_p");
  memory->create(xm,natoms*3,"mode:xm");
  memory->create(vm_p,natoms*3,"mode:vm_p");
  memory->create(vm,natoms*3,"mode:vm");
  for (int n=0; n<3*natoms; n++){
    xm[n]=0.0;
    xm_p[n]=0.0;
    vm[n]=0.0;
    vm_p[n]=0.0;
  }
  memory->create(em,natoms*3,"mode:em");
  memory->create(u_p,natoms,3,"mode:u_p");
  for (int i=0; i<natoms; i++){
    for (int a=0; a<3; a++){
      u_p[i][a] = 0.0;
    }
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

  std::ifstream readfile2;

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

  std::ifstream readfile3;

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

  // Allocate other things.
  memory->create(vm, 3*natoms, "wp:vm");
  memory->create(xm, 3*natoms, "wp:xm");
  memory->create(tm, 3*natoms, "wp:tm");
  for (int i=0; i<3*natoms; i++){
    vm[i]=0.0;
    xm[i]=0.0;
    tm[i]=0.0;
  }


}

/* ---------------------------------------------------------------------- */

void FixWp::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixWp::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWp::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  /*
  for (int i=0; i<32; i++){
    fprintf(fh_d, "%f %f %f\n", v[i][0],v[i][1],v[i][2]);
  }
  */
  double *mass = atom->mass;
  int *type = atom->type;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;
  int me = comm->me; // If called, me is the varibale to give us the number of processor that performs the command.

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

  //foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  //force_flag = 0;

  /*
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
      if (xstyle) f[i][0] += xvalue;
      if (ystyle) f[i][1] += yvalue;
      if (zstyle) f[i][2] += zvalue;
    }

  }
  */

  // Calculate current temperature
  /*
  double t_current = 0.0;
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit){
      t_current += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
        mass[type[i]];
    }
  }
  */

  // Excite modes
  double modeamp, forceamp;
  double pi = 3.1415926535;
  double kb = 1.38064852e-23;
  double temperature = 1e-18; // 1e-17 was easily visible
  //double forceamp = 2*pi*freq[modeindx]*1e12*sqrt(kb*temperature);
  
  // When declaring modeindx, realize it starts at ZERO! (so subtract 1 from mode ID).
  //modeindx = 1906; // 1906 - interface mode, 1907 is the periodic interface
  //modeindx = 1647-1; // 10 THz, partially extended on heavy side
  /*
  modeindx = 1906-1;
  //modeindx = 100; // 100 - extended low freq mode
  modeamp = sqrt((kb*temperature)/(freq[modeindx]*freq[modeindx]*2*pi*2*pi*1e12*1e12));
  forceamp = 2*pi*freq[modeindx]*1e12*sqrt(kb*temperature);
  //printf("forceamp: %e\n", forceamp);
  time = (update->dt*update->ntimestep)*1e-12;
  
  for (int i = 0; i < nlocal; i++){

    double ms = mass[type[i]]*(1.0/6.0221409e+23)*1e-3;
    // Convert from J/m to eV/A with *1.6022e-9
    f[i][0] += (1.0/sqrt(ms))*emat[3*i+0][modeindx]*forceamp*sin(2*pi*freq[modeindx]*1e12*time)*1.6022e-9;
    f[i][1] += (1.0/sqrt(ms))*emat[3*i+1][modeindx]*forceamp*sin(2*pi*freq[modeindx]*1e12*time)*1.6022e-9;
    f[i][2] += (1.0/sqrt(ms))*emat[3*i+2][modeindx]*forceamp*sin(2*pi*freq[modeindx]*1e12*time)*1.6022e-9;
    
    // Rescale velocities
    //v[i][0] = v[i][0]*sqrt(temperature/t_current);
    //v[i][1] = v[i][1]*sqrt(temperature/t_current);
    //v[i][2] = v[i][2]*sqrt(temperature/t_current);
    
  }
  */

  // Sinusoidal driving force
  /*
  double sinterm;
  double massconv=(1.0/6.0221409e+23)*1e-3;
  double ex,ey,ez;
  double factor;
  double massfactor;
  //for (int n=0; n<3*natoms; n++){
    //if (1.54 < freq[n] && freq[n]<1.541280572309643127e+01){
      //printf("yes\n");
      //modeindx = n;
      modeindx=1962-1;
      //modeamp = sqrt((kb*temperature)/(freq[modeindx]*freq[modeindx]*2*pi*2*pi*1e12*1e12));
      //forceamp = 2*pi*freq[modeindx]*1e12*sqrt(kb*temperature);
      //printf("forceamp: %e\n", forceamp);
      time = (update->dt*update->ntimestep)*1e-12;
      //sinterm = sin(2*pi*freq[modeindx]*1e12*time)*1.6022e-9;
      factor = 2*pi*freq[modeindx]*1e12*sqrt(kb*temperature)*sin(2*pi*freq[modeindx]*1e12*time)*1.6022e-9;
      for (int i = 0; i < nlocal; i++){

        //double ms = mass[type[i]]*massconv;
        massfactor = (1.0/sqrt(mass[type[i]]*massconv));
        // Convert from J/m to eV/A with *1.6022e-9
        f[i][0] += massfactor*emat[3*i+0][modeindx]*factor;
        f[i][1] += massfactor*emat[3*i+1][modeindx]*factor;
        f[i][2] += massfactor*emat[3*i+2][modeindx]*factor;
        
        // Rescale velocities
        //v[i][0] = v[i][0]*sqrt(temperature/t_current);
        //v[i][1] = v[i][1]*sqrt(temperature/t_current);
        //v[i][2] = v[i][2]*sqrt(temperature/t_current);
        
      }
    //}
  //}
  */

  // Initialize mode temperatures.
  
  //if (me == 0){
    //printf(" me: %d\n", me);
  if (update->ntimestep==0){
      printf(" ******* Creating WP on proc %d. *******\n", me);

      //double a0=5.431;
      k0 = k0*(2.0*3.14159265359/a0);
      double eta = 1./(xi);

      double ux,uy,uz;
      double vx,vy,vz;

      double zp,v1,v2,v3,v4;
      //double z0 = 135.775;
      for (int i=0; i<nlocal; i++){

        // Calculate displacement.
        /* 
        Python script:
        zp = z-z0
        v1 = -1.0*eta*eta*zp*zp;
        v2 = k0*zp;
        v3 = 2.0*eta*eta*vg*zp;
        v4 = -1.0*k0*vg;
        u = amp*np.exp(v1)*np.cos(v2);
        */
        
        zp = x[i][0] - z0;
        v1 = -1.0*eta*eta*zp*zp;
        v2 = k0*zp;
        v3 = 2.0*eta*eta*vg*zp;
        v4 = -1.0*k0*vg;
        uz = amp*exp(v1)*cos(v2);


        // Calculate velocity
        /*
        Python script:
        zp = z-z0
        v1 = -1.0*eta*eta*zp*zp;
        v2 = k0*zp;
        v3 = 2.0*eta*eta*vg*zp;
        v4 = -1.0*k0*vg;
        v = amp*np.exp(v1)*(v3*np.cos(v2)-v4*np.sin(v2));
        */
        v1 = -1.0*eta*eta*zp*zp;
        v2 = k0*zp;
        v3 = 2.0*eta*eta*vg*zp;
        v4 = -1.0*k0*vg;
        vz = amp*exp(v1)*(v3*cos(v2)-v4*sin(v2));

        // Add to position and velocity
        //x[i][2] += uz;
        //v[i][2] = vz;
        up[i][0] = uz;
        up[i][1] = 0;
        up[i][2] = 0;
        vp[i][0] = vz;
        vp[i][1] = 0;
        vp[i][2] = 0;
        
        //if (uz>0 || vz>0){
        // fprintf(fh_d, "%e %e\n", uz,vz);
        //}

        
        
      }

      
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

      // Find mode coordinates that give wave-packet, and keep coordinates above a certain tolerance.
      for (int n=0; n<3*natoms; n++){
      //for (int n=7; n<8; n++){
        //printf("n: %d\n", n);
        xm_p[n]=0.0;
        xm[n] = 0.0; // zero the total coordinates
        vm_p[n]=0.0;
        vm[n] = 0.0; // zero the total velocities
        for (int i=0; i<nlocal; i++){
          if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
            for (int a=0; a<3; a++){
              //u[i][a] = u[i][a]*1e-10;
              double ms = mass[type[i]]; //mass[type[i]]*(1.0/6.0221409e+23)*1e-3;
              double masskg = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
              //xm_p[n] += sqrt(mass[type[i]])*emat[3*i+a][n]*u[i][a]*1e-10*(1.0/6.0221409e+23)*1e-3; // Convert to SI units.;
              //xm_p[n] += sqrt(ms)*emat[3*i+a][n]*u_p[i][a]; // sqrt(kg/kmol) * Angstrom
              xm_p[n] += sqrt(masskg)*emat[3*(tag[i]-1)+a][n]*up[i][a]; // sqrt(kg) * Angstrom
              //if (u_p[i][a] != 0) printf("%e\n", u_p[i][a]);
              //fmi_p[natoms*n+i] = sqrt(ms)*emat[3*i+a][n]*f[i][a];
              //vm_p[n] += sqrt(masskg)*emat[3*i+a][n]*v[i][a]*100; // sqrt(kg) * m/s, since 1 A/ps = 100 m/s.
              vm_p[n] += sqrt(masskg)*emat[3*(tag[i]-1)+a][n]*vp[i][a]; // sqrt(kg) * A/ps
              
            }
          }
        }
      }

      MPI_Allreduce(xm_p,xm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
      MPI_Allreduce(vm_p,vm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
      

      
      // If mode energy < 0.001 eV, keep its coordinates.

      // Mode energies
      // Only modes with energy greater than certain tolerance are considered WP modes.
      // For amplitude of 2.5e-2, that tolerance is about 1e-3 eV.
      double pi = 3.1415926535;
      double pem = 0.0; 
      for (int n=3; n<3*natoms; n++){
        em[n] = (0.5*vm[n]*vm[n]*100*100 + 0.5*freq[n]*freq[n]*2.0*3.1415926535*2.0*3.1415926535*1e12*1e12*xm[n]*xm[n]*1e-10*1e-10)*6.242e+18; // eV, since 6.242e+18 eV = 1 J
        //fprintf(fh_d, "%e\n", em[n]);
      }

      /*
      printf(" ******* Setting mode coordinates on proc %d. *******\n", me);
      //printf(" %d atoms.\n", natoms);
      int modeindx;
      double modetemp;
      double modecoor;
      string line;

      ifstream fh_mt("MODETEMP");
      ifstream fh_mc("MODECOOR");

      //printf(" Reading MODETEMP and MODECOOR.\n");

      if (!fh_mt.is_open()) {
          printf("Unable to open MODETEMP!\n");
          exit(1);
      }

      if (!fh_mc.is_open()) {
          printf("Unable to open MODECOOR!\n");
          exit(1);
      }

      // Reading mode temperatures.
      while (getline(fh_mt, line)){
          if (line.find('#') != std::string::npos){
              //printf("Ignore this line\n");
          }
          else{

              stringstream ss3(line);
              ss3 >> modeindx >> modetemp;
              //printf(" Setting mode %d to %f K.\n", modeindx, modetemp);
              //vm[modeindx]=sqrt(2*1.38064852e-23*modetemp);
              if (em[modeindx] < 0.001){
                tm[modeindx]=modetemp;
              }
              else{
                tm[modeindx]=0.0;
              }
              //if (modeindx==2209){
              //    printf(" vm[modeindx]: %e\n", vm[modeindx]);
              //}

          }
      }
      fh_mt.close();

      // Reading mode coordinates.
      while (getline(fh_mc, line)){
          if (line.find('#') != std::string::npos){
              //printf("Ignore this line\n");
          }
          else{

              stringstream ss3(line);
              ss3 >> modeindx >> modecoor;
              //printf(" Setting mode %d to %e\n", modeindx, modecoor);
              if (em[modeindx]<0.001){
                xm[modeindx]=modecoor;
              }
              else{
                xm[modeindx]=0.0;
              }
              //if (modeindx==2209){
              //    printf(" vm[modeindx]: %e\n", vm[modeindx]);
              //}

          }
      }
      fh_mc.close();

      printf("%d atoms on proc %d.\n", nlocal, me);

      // Solve for mode velocities given mode energy is kT and potential energy is known from mode coordinate.
      for (int n=0; n<3*natoms; n++){
       
        if (em[n]>0.001){
          vm[n] = sqrt( (2*kb*tm[n]) - (freq[n]*freq[n]*1e12*1e12*4.0*pi*pi*xm[n]*xm[n]) );
        }
        else{
          vm[n]=0.0;
        }

      }
      */
      // Set atomic velocities and coordinates based on these mode velocities and coordinates.
      
      double velocity_factor = 1; // Units are already A/ps
      double coord_factor = 1; // Units are already A
      for (int n=0; n<3*natoms; n++){
        if (em[n] > 1e-3){
          if (me==0) printf("%d\n", n);
          for (int i=0; i<nlocal; i++){
            double amass = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
            double invsqrtmass = 1.0/sqrt(amass);
            //printf(" amass: %e\n", amass);
            for (int a=0; a<3; a++){
              
              v[i][a] += invsqrtmass*vm[n]*emat[3*(tag[i]-1)+a][n]*velocity_factor;
              x[i][a] += invsqrtmass*xm[n]*emat[3*(tag[i]-1)+a][n]*coord_factor;
              //fprintf(fh_d, "%d %d %d %e %e\n", n,i,a,invsqrtmass*vm[n]*emat[3*(tag[i]-1)+a][n]*velocity_factor,invsqrtmass*xm[n]*emat[3*(tag[i]-1)+a][n]*coord_factor);

            }

          }
        }
      }
      
      

  } // if (update->ntimestep==0){
  
}

/* ---------------------------------------------------------------------- */

void FixWp::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // set force to desired value on requested level, 0.0 on other levels
  if (ilevel == ilevel_respa) post_force(vflag);
  else {
    Region *region = NULL;
    if (iregion >= 0) {
      region = domain->regions[iregion];
      region->prematch();
    }

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if (xstyle) f[i][0] = 0.0;
        if (ystyle) f[i][1] = 0.0;
        if (zstyle) f[i][2] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixWp::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixWp::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixWp::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*3 * sizeof(double);
  return bytes;
}
