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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_tep.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "update.h"
#include <math.h>       /* tanh, log */

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#define MAXLINE 1024

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

PairTep::PairTep(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;

  fh_debug = fopen("D_TEP", "w");
  fh_pe = fopen("pe.dat", "w");
  fh_pe3 = fopen("pe3.dat", "w");
  fh_ht = fopen("ht.dat", "w");
  fh_ht2 = fopen("ht2.dat", "w");
  fh_ht3 = fopen("ht3.dat", "w");
}

/* ---------------------------------------------------------------------- */

PairTep::~PairTep()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(morse1);
    memory->destroy(offset);

    memory->destroy(fc2);
    if (order >2) memory->destroy(fc3);
    if (order >3) memory->destroy(fc4);
    memory->destroy(phi);
    memory->destroy(emat);

    memory->destroy(u_p);
    memory->destroy(u);
    memory->destroy(xm_p);
    memory->destroy(xm);
    memory->destroy(vm_p);
    memory->destroy(vm);
    memory->destroy(fm_p);
    memory->destroy(fm);
    memory->destroy(psi);

  }

  fclose(fh_debug);
  fclose(fh_pe);
  fclose(fh_pe3);
  fclose(fh_ht);
  fclose(fh_ht2);
  fclose(fh_ht3);
}

/* ---------------------------------------------------------------------- */

void PairTep::compute(int eflag, int vflag)
{
  int i,j,ii,jj,kk,inum,jnum,itype,jtype;
  int a,b;
  int k,c;
  int itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,dr,dexp,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double fc;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  /*
  for (int i=0; i<8; i++){
    fprintf(fh_debug, "%f %f %f\n", v[i][0],v[i][1],v[i][2]);

  }
  */
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  int natoms = atom->natoms;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double boxlo_x, boxlo_y, boxlo_z, boxhi_x, boxhi_y, boxhi_z;
  boxlo_x = domain->boxlo[0];
  boxlo_y = domain->boxlo[1];
  boxlo_z = domain->boxlo[2];
  boxhi_x = domain->boxhi[0];
  boxhi_y = domain->boxhi[1];
  boxhi_z = domain->boxhi[2];

  double xi,yi,zi;
  double xj,yj,zj;
  double uix,uiy,uiz;
  double ujx,ujy,ujz;
  double ukx,uky,ukz;
  double ux,uy,uz;
  double uia, uib,ujb;
  // loop over all atoms

  //fh_umat = fopen("UMAT", "w");

  double pe_a = 0.0;
  double pe_b = 0.0;
  double pe = 0.0;
  double ht = 0.0;
  double pe3 = 0.0;
  double pe3_a = 0.0;
  double pe3_b = 0.0;
  double ht3 = 0.0;

  for (int w=0; w<nfc2; w++){

    /*
    i = fc2[w].i;
    j = fc2[w].j;
    a = fc2[w].a;
    b = fc2[w].b;
    fc = fc2[w].fc;

    itag = tag[i]-1;
    jtag = tag[j]-1;
    itype = type[itag];
    jtype = type[jtag];

    uia = x[i][a] - x0[itag][a];
    uib = x[i][b] - x0[itag][b];

    if (std::abs(uia) > boxhi_x/2.0 && uia > 0) uia -= boxhi_x;
    else if (std::abs(uia) > boxhi_x/2.0 && uia < 0) uia += boxhi_x;

    if (std::abs(uib) > boxhi_x/2.0 && uib > 0) uib -= boxhi_x;
    else if (std::abs(uib) > boxhi_x/2.0 && uib < 0) uib += boxhi_x;

    jtag = tag[j]-1;

    ujb = x[j][b] - x0[jtag][b];

    if (std::abs(ujb) > boxhi_x/2.0 && ujb > 0) ujb -= boxhi_x;
    else if (std::abs(ujb) > boxhi_x/2.0 && ujb < 0) ujb += boxhi_x;

    // TEP forces
    f[i][a] -= fc*ujb;

    // TEP energy
    pe += 0.5*fc*uia*ujb;
    if (type[itag]==1){
      pe_a += 0.5*fc*uia*ujb;
    }
    else if (type[itag]==2){
      pe_b += 0.5*fc*uia*ujb;
    }
    else{
      printf("FUCK!\n");
    }

    if ( (i==1 && j==17) || (i==17 || j==1) ){
      fprintf(fh_debug, "%d %d %d %d\n", i,a,j,b);
      fprintf(fh_debug, "  %f\n", fc);
    }
    */

    // Power transfer from A to B.
    /*
    // If i is in A, and j is in B, then compute Q_ij and add to total Q_AB.
    if (itype==1 && jtype==2){ 
        //printf("%d %d\n", tag[i],tag[j]);
        //for (int a=0; a<3; a++){
          //Q -= 0.5*( fij[a]*(v[i][a] + v[j][a]) );
          //Q += 0.5*( fij[a]*(v[i][a] + v[j][a]) );
        //}
      ht += 1.0*fc*( (v[i][a]*ujb) - (v[j][a]*uib) );
    }

    // If i is in B, and j is in A, then compute Q_ji and add to total Q_AB.
    
    else if (itype==2 && jtype==1){ 
        //printf("%d %d\n", tag[i],tag[j]);
        //for (int a=0; a<3; a++){
          //Q += 0.5*( fij[a]*(v[i][a] + v[j][a]) );
        //}
      ht -= 1.0*fc*( (v[i][a]*ujb) - (v[j][a]*uib) );
    }
    */

    // This expression integrates to zero!
    //ht += 1.0*fc*( (v[i][a]*ujb) - (v[j][a]*uib) );

    // TITEP forces
    //f[i][a] += fc*(uib-ujb);


    // TITEP forces
    //f[i][0] += 1.0*(fc2[ii][jj][0]*(uix-ujx) + fc2[ii][jj][1]*(uiy-ujy) + fc2[ii][jj][2]*(uiz-ujz));
    //f[ii][1] += 1.0*(fc2[ii][jj][3]*(uix-ujx) + fc2[ii][jj][4]*(uiy-ujy) + fc2[ii][jj][5]*(uiz-ujz));
    //f[ii][2] += 1.0*(fc2[ii][jj][6]*(uix-ujx) + fc2[ii][jj][7]*(uiy-ujy) + fc2[ii][jj][8]*(uiz-ujz));

    
  }

  //printf("%e %e %e %e\n", pe_a, pe_b, pe_a+pe_b, pe);


  double fij[3],fji[3];
  double qij=0.0;

  
  for (int w=0; w<nfc2; w++){

    i = fc2[w].i;
    j = fc2[w].j;
    a = fc2[w].a;
    b = fc2[w].b;
    fc = fc2[w].fc;

    itype = type[i];
    jtype = type[j];

    uia = x[i][a] - x0[i][a];
    uib = x[i][b] - x0[i][b];

    if (std::abs(uia) > boxhi_x/2.0 && uia > 0) uia -= boxhi_x;
    else if (std::abs(uia) > boxhi_x/2.0 && uia < 0) uia += boxhi_x;

    if (std::abs(uib) > boxhi_x/2.0 && uib > 0) uib -= boxhi_x;
    else if (std::abs(uib) > boxhi_x/2.0 && uib < 0) uib += boxhi_x;

    jtag = tag[j]-1;

    ujb = x[j][b] - x0[jtag][b];

    if (std::abs(ujb) > boxhi_x/2.0 && ujb > 0) ujb -= boxhi_x;
    else if (std::abs(ujb) > boxhi_x/2.0 && ujb < 0) ujb += boxhi_x;

    // TEP forces
    f[i][a] -= fc*ujb;
    // TITEP forces
    //f[i][a] += fc*(uib-ujb);



  } // for (int w=0; w<nfc2; w++){


  double ukc;
  if (order > 2){
    for (int w=0; w<nfc3; w++){

      i = fc3[w].i;
      j = fc3[w].j;
      a = fc3[w].a;
      b = fc3[w].b;
      k = fc3[w].k;
      c = fc3[w].c;
      fc = fc3[w].fc;

      itype = type[i];
      jtype = type[j];

      uia = x[i][a] - x0[i][a];
      uib = x[i][b] - x0[i][b];
      ukc = x[k][c] - x0[k][c];

      if (std::abs(uia) > boxhi_x/2.0 && uia > 0) uia -= boxhi_x;
      else if (std::abs(uia) > boxhi_x/2.0 && uia < 0) uia += boxhi_x;

      if (std::abs(uib) > boxhi_x/2.0 && uib > 0) uib -= boxhi_x;
      else if (std::abs(uib) > boxhi_x/2.0 && uib < 0) uib += boxhi_x;

      jtag = tag[j]-1;

      ujb = x[j][b] - x0[jtag][b];

      if (std::abs(ujb) > boxhi_x/2.0 && ujb > 0) ujb -= boxhi_x;
      else if (std::abs(ujb) > boxhi_x/2.0 && ujb < 0) ujb += boxhi_x;

      if (std::abs(ukc) > boxhi_x/2.0 && ukc > 0) ukc -= boxhi_x;
      else if (std::abs(ukc) > boxhi_x/2.0 && ukc < 0) ukc += boxhi_x;

      // TEP forces
      f[i][a] -= 0.5*fc*ujb*ukc;
      


    } // for (int w=0; w<nfc3; w++){
  } // if order > 2

  double uld;
  int l,d;
  int ktype;
  if (order > 3){
    for (int w=0; w<nfc4; w++){

      i = fc4[w].i;
      j = fc4[w].j;
      a = fc4[w].a;
      b = fc4[w].b;
      k = fc4[w].k;
      c = fc4[w].c;
      l = fc4[w].l;
      d = fc4[w].d;
      fc = fc4[w].fc;

      itype = type[i];
      jtype = type[j];
      ktype = type[k];

      uia = x[i][a] - x0[i][a];
      uib = x[i][b] - x0[i][b];
      ukc = x[k][c] - x0[k][c];
      uld = x[l][d] - x0[l][d];

      if (std::abs(uia) > boxhi_x/2.0 && uia > 0) uia -= boxhi_x;
      else if (std::abs(uia) > boxhi_x/2.0 && uia < 0) uia += boxhi_x;

      if (std::abs(uib) > boxhi_x/2.0 && uib > 0) uib -= boxhi_x;
      else if (std::abs(uib) > boxhi_x/2.0 && uib < 0) uib += boxhi_x;

      jtag = tag[j]-1;

      ujb = x[j][b] - x0[jtag][b];

      if (std::abs(ujb) > boxhi_x/2.0 && ujb > 0) ujb -= boxhi_x;
      else if (std::abs(ujb) > boxhi_x/2.0 && ujb < 0) ujb += boxhi_x;

      if (std::abs(ukc) > boxhi_x/2.0 && ukc > 0) ukc -= boxhi_x;
      else if (std::abs(ukc) > boxhi_x/2.0 && ukc < 0) ukc += boxhi_x;

      if (std::abs(uld) > boxhi_x/2.0 && uld > 0) uld -= boxhi_x;
      else if (std::abs(uld) > boxhi_x/2.0 && uld < 0) uld += boxhi_x;

      // TEP forces
      f[i][a] -= (1.0/6.0)*fc*ujb*ukc*uld;
      


    } // for (int w=0; w<nfc4; w++){
  } // if order > 3

  // Compute mode coordinates and velocities.
  double *mass = atom->mass;
  int *mask = atom->mask;

  // Compute atomic displacements.
  //fprintf(fh_disp, "\n");-
  //fprintf(fh_disp, "Timestep: %d\n", update->ntimestep);
  //fprintf(fh_disp, "\n");
  for (int i=0; i<nlocal; i++){
    //if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
      //u_p[i]=0.0;
      //u[i] = 0.0; // zero the total coordinates
      for (int a=0; a<3; a++){
        u[i][a] = (x[i][a]-x0[i][a]);
        //fprintf(fh_debug, "%f\n", u_p[i][a]);
        if (std::abs(u[i][a]) > domain->boxhi[a]/2.0 && u[i][a] > 0) u[i][a] -= domain->boxhi[a];
        else if (std::abs(u[i][a]) > domain->boxhi[a]/2.0 && u[i][a] < 0) u[i][a] += domain->boxhi[a];

        //u_p[i][a] = u_p[i][a]*1e-10;

        //fprintf(fh_disp,"%e ", x[i][a]-x0[i][a]);
      }
      //fprintf(fh_disp, "%e %e\n", x0[i][2],u_p[i][2]);
    //}
  }
  // NEED TO MAKE DISPLACEMENT VECTOR 1D!!!!!!!!
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
    for (int i=0; i<natoms; i++){
      //if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
        for (int a=0; a<3; a++){
          //u[i][a] = u[i][a]*1e-10;
          double ms = mass[type[i]]; //mass[type[i]]*(1.0/6.0221409e+23)*1e-3;
          double masskg = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
          //xm_p[n] += sqrt(mass[type[i]])*emat[3*i+a][n]*u[i][a]*1e-10*(1.0/6.0221409e+23)*1e-3; // Convert to SI units.;
          //xm_p[n] += sqrt(ms)*emat[3*i+a][n]*u_p[i][a]; // sqrt(kg/kmol) * Angstrom
          //xm_p[n] += sqrt(masskg)*emat[3*i+a][n]*u_p[i][a]; // sqrt(kg) * Angstrom
          xm_p[n] += sqrt(masskg)*emat[3*i+a][n]*u[i][a]; // sqrt(kg) * Angstrom
          //fprintf(fh_debug, "%e\n", u_p[i][a]);
          fm_p[n] += sqrt(ms)*emat[3*i+a][n]*f[i][a]; // sqrt(kg/kmol) * eV/A
          //fmi_p[natoms*n+i] = sqrt(ms)*emat[3*i+a][n]*f[i][a];
          //vm_p[n] += sqrt(masskg)*emat[3*i+a][n]*v[i][a]*100; // sqrt(kg) * m/s, since 1 A/ps = 100 m/s.
          vm_p[n] += sqrt(masskg)*emat[3*i+a][n]*v[i][a]; // sqrt(kg) * A/ps

          /*
          if (n==7){
            
            fprintf(fh_xm, "n,i,a: %d,%d,%d\n", n,i,a);
            fprintf(fh_xm, "  mass: %e\n", ms);
            fprintf(fh_xm, "  emat[3*i+a][n]: %e\n", emat[3*i+a][n]);
            fprintf(fh_xm, "  u: %e\n", u[i][a]);
            fprintf(fh_xm, "  %e\n", xm_p[n]);
            
          }
          */
          
        }
      //}
    }

  }

  MPI_Allreduce(xm_p,xm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(vm_p,vm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(fm_p,fm,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  //MPI_Allreduce(fmi_p,fmi,3*natoms*natoms,MPI_DOUBLE,MPI_SUM,world);

  /*
  for (ii = 0; ii < natoms; ii++) {

    itag = tag[ii]-1;

    itype = type[i];

    uix = x[ii][0] - x0[itag][0];
    uiy = x[ii][1] - x0[itag][1];
    uiz = x[ii][2] - x0[itag][2];

    if (std::abs(uix) > boxhi_x/2.0 && uix > 0) uix -= boxhi_x;
    else if (std::abs(uix) > boxhi_x/2.0 && uix < 0) uix += boxhi_x;

    if (std::abs(uiy) > boxhi_y/2.0 && uiy > 0) uiy -= boxhi_y;
    else if (std::abs(uiy) > boxhi_y/2.0 && uiy < 0) uiy += boxhi_y;

    if (std::abs(uiz) > boxhi_z/2.0 && uiz > 0) uiz -= boxhi_z;
    else if (std::abs(uiz) > boxhi_z/2.0 && uiz < 0) uiz += boxhi_z;

    
    
    //fprintf(fh_debug, "x[%d]: %f %f %f\n", itag, x[ii][0],x[ii][1],x[ii][2]);
    //fprintf(fh_debug, "x0[%d]: %f %f %f\n", itag, x0[itag][0],x0[itag][1],x0[itag][2]);
    //fprintf(fh_debug, "u[%d]: %f %f %f\n", itag, uix, uiy, uiz);
    //fprintf(fh_debug, "\n");
    
    

    for (jj = 0; jj < natoms; jj++) {

      jtag = tag[jj]-1;

      ujx = x[jj][0] - x0[jtag][0];
      ujy = x[jj][1] - x0[jtag][1];
      ujz = x[jj][2] - x0[jtag][2];

      if (std::abs(ujx) > boxhi_x/2.0 && ujx > 0) ujx -= boxhi_x;
      else if (std::abs(ujx) > boxhi_x/2.0 && ujx < 0) ujx += boxhi_x;

      if (std::abs(ujy) > boxhi_y/2.0 && ujy > 0) ujy -= boxhi_y;
      else if (std::abs(ujy) > boxhi_y/2.0 && ujy < 0) ujy += boxhi_y;

      if (std::abs(ujz) > boxhi_z/2.0 && ujz > 0) ujz -= boxhi_z;
      else if (std::abs(ujz) > boxhi_z/2.0 && ujz < 0) ujz += boxhi_z;



      // TEP forces
      
      //f[ii][0] += -1.0*fc2[ii][jj][0]*ujx - 1.0*fc2[ii][jj][1]*ujy - 1.0*fc2[ii][jj][2]*ujz;
      //f[ii][1] += -1.0*fc2[ii][jj][3]*ujx - 1.0*fc2[ii][jj][4]*ujy - 1.0*fc2[ii][jj][5]*ujz;
      //f[ii][2] += -1.0*fc2[ii][jj][6]*ujx - 1.0*fc2[ii][jj][7]*ujy - 1.0*fc2[ii][jj][8]*ujz;
      

      // TITEP forces
      
      f[ii][0] += 1.0*(fc2[ii][jj][0]*(uix-ujx) + fc2[ii][jj][1]*(uiy-ujy) + fc2[ii][jj][2]*(uiz-ujz));
      f[ii][1] += 1.0*(fc2[ii][jj][3]*(uix-ujx) + fc2[ii][jj][4]*(uiy-ujy) + fc2[ii][jj][5]*(uiz-ujz));
      f[ii][2] += 1.0*(fc2[ii][jj][6]*(uix-ujx) + fc2[ii][jj][7]*(uiy-ujy) + fc2[ii][jj][8]*(uiz-ujz));
      
     

      if (eflag) {
        //evdwl = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
        //evdwl *= factor_lj;
        evdwl = 0.5*(fc2[ii][jj][0]*uix*ujx + fc2[ii][jj][1]*uix*ujy + fc2[ii][jj][2]*uix*ujz 
              + fc2[ii][jj][3]*uiy*ujx + fc2[ii][jj][4]*uiy*ujy + fc2[ii][jj][5]*uiy*ujz
              + fc2[ii][jj][6]*uiz*ujx + fc2[ii][jj][7]*uiz*ujy + fc2[ii][jj][8]*uiz*ujz);
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);

      // Compute 3rd order terms
      if (order>2){

        for (kk = 0; kk < natoms; kk++) {

          ukx = x[kk][0] - x0[kk][0];
          uky = x[kk][1] - x0[kk][1];
          ukz = x[kk][2] - x0[kk][2];

          if (std::abs(ukx) > boxhi_x/2.0 && ukx > 0) ukx -= boxhi_x;
          else if (std::abs(ukx) > boxhi_x/2.0 && ukx < 0) ukx += boxhi_x;

          if (std::abs(uky) > boxhi_y/2.0 && uky > 0) uky -= boxhi_y;
          else if (std::abs(uky) > boxhi_y/2.0 && uky < 0) uky += boxhi_y;

          if (std::abs(ukz) > boxhi_z/2.0 && ukz > 0) ukz -= boxhi_z;
          else if (std::abs(ukz) > boxhi_z/2.0 && ukz < 0) ukz += boxhi_z;

          //ux = uix - 0.5*(ujx+ukx);
          //uy = uiy - 0.5*(ujy+uky);
          //uz = uiz - 0.5*(ujz+ukz);

          
          f[ii][0] += -0.5*(fc3[ii][jj][kk][0]*ujx*ukx + fc3[ii][jj][kk][1]*ujx*uky + fc3[ii][jj][kk][2]*ujx*ukz \
                          + fc3[ii][jj][kk][3]*ujy*ukx + fc3[ii][jj][kk][4]*ujy*uky + fc3[ii][jj][kk][5]*ujy*ukz \
                          + fc3[ii][jj][kk][6]*ujz*ukx + fc3[ii][jj][kk][7]*ujz*uky + fc3[ii][jj][kk][8]*ujz*ukz);

          f[ii][1] += -0.5*(fc3[ii][jj][kk][9]*ujx*ukx + fc3[ii][jj][kk][10]*ujx*uky + fc3[ii][jj][kk][11]*ujx*ukz \
                          + fc3[ii][jj][kk][12]*ujy*ukx + fc3[ii][jj][kk][13]*ujy*uky + fc3[ii][jj][kk][14]*ujy*ukz \
                          + fc3[ii][jj][kk][15]*ujz*ukx + fc3[ii][jj][kk][16]*ujz*uky + fc3[ii][jj][kk][17]*ujz*ukz);

          f[ii][2] += -0.5*(fc3[ii][jj][kk][18]*ujx*ukx + fc3[ii][jj][kk][19]*ujx*uky + fc3[ii][jj][kk][20]*ujx*ukz \
                          + fc3[ii][jj][kk][21]*ujy*ukx + fc3[ii][jj][kk][22]*ujy*uky + fc3[ii][jj][kk][23]*ujy*ukz \
                          + fc3[ii][jj][kk][24]*ujz*ukx + fc3[ii][jj][kk][25]*ujz*uky + fc3[ii][jj][kk][26]*ujz*ukz);
          
          


          
        }

      }

    }
  }

  // Total forces on cell
  double fx, fy, fz;
  fx =0;
  fy = 0;
  fz=0;
  for (ii=0; ii<natoms; ii++){
    fx += f[ii][0];
    fy += f[ii][1];
    fz += f[ii][2];
  }

  //printf("%f %f %f\n", fx,fy,fz);
  */

  if (vflag_fdotr) virial_fdotr_compute();

  //fclose(fh_umat);
  //printf("ASDF~_~_~_~_~_~\n");

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTep::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(morse1,n+1,n+1,"pair:morse1");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTep::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  //order = force->numeric(FLERR,arg[0]);
  order = atoi(arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = 10.0;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTep::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  //force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  //force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_one = 10.0;
  //if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");



  int natoms = atom->natoms;
  /*
  memory->create(fc2,natoms,natoms,9,"pair:fc2");
  for (int i=0; i<natoms; i++){
    for (int j=0; j<natoms; j++){
      for (int indx=0; indx<9; indx++){
        fc2[i][j][indx] = 0.0;
      }
    }
  }

  if (order > 2){
    memory->create(fc3,natoms,natoms,natoms,27,"pair:fc3");
    for (int i=0; i<natoms; i++){
      for (int j=0; j<natoms; j++){
        for (int k=0; k<natoms; k++){
          for (int indx=0; indx<27; indx++){
            fc3[i][j][k][indx] = 0.0;
          }
        }
      }
    }
  }
  */
  memory->create(x0,natoms,3,"pair:x0");
  for (int i = 0; i < natoms; i++) {
    for (int j = 0; j < 3; j++) {
      x0[i][j] = 0.0;
    }
  }

  read_fc2();
  nfc3=0;
  nfc4=0;
  if (order > 2){
    read_fc3();
  }
  if (order > 3){
    read_fc4();
  }
  //printf("ASDF------------\n");
  read_equil();
  //printf("ASDF------\n");
  /* EMAT */

  memory->create(emat,natoms*3,natoms*3,"mode:emat");

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
          //fprintf(fh_debug, "%e\n", emat[i][j]);
      }
  }  

  // Allocate mode arrays
  memory->create(u_p,natoms,3,"mode:u_p");
  memory->create(u,natoms,3,"mode:u");
  memory->create(xm_p,natoms*3,"mode:xm_p");
  memory->create(xm,natoms*3,"mode:xm");
  memory->create(vm_p,natoms*3,"mode:vm_p");
  memory->create(vm,natoms*3,"mode:vm");
  memory->create(fm_p,natoms*3,"mode:fm_p");
  memory->create(fm,natoms*3,"mode:fm");

}

/* ----------------------------------------------------------------------
   read fc2
------------------------------------------------------------------------- */

void PairTep::read_fc2()
{


  /*
  char line[MAXLINE],*ptr;
  FILE *fp;
  int n,nwords,nuc_poscar;
  char **words = new char*[MAXLINE];
  int info = 0;
  int i,j,a,b;
  double fc;

  fp = fopen("FC2", "r");

  if (fp == NULL) {
    char str[MAXLINE];
    sprintf(str,"Cannot open the file %s", "FC2");
    error->all(FLERR,str);
  }

  int eof = 0;
  
  // extract the fc2s
  
  nfc2=-1; // minus 1 since first line is ignored
  while (!eof) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    //MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    //MPI_Bcast(&n,1,MPI_INT,0,world);
    //MPI_Bcast(line,n,MPI_CHAR,0,world);

    //std::cout << line << std::endl;
    nfc2++;
  }

  printf(" Found %d FC2s.\n", nfc2);
  
  memory->create(fc2,nfc2,"pair:fc2");
  
  
  fp = fopen("FC2", "r");

  if (fp == NULL) {
    char str[MAXLINE];
    sprintf(str,"Cannot open the file %s", "FC2");
    error->all(FLERR,str);
  }

  // extract the fc2s
  
  eof=0;
  int indx=-1; 
  while (!eof) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
        printf("ASDF\n");
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    //std::cout << line << std::endl;
    // strip comment, skip line if blank
    
    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    //nwords = atom->count_words(line);
    nwords=5;
    //fprintf(fh_debug, "nwords: %d\n", nwords);
    if (nwords == 0) {
      continue;
    } else {

      std::cout << line << std::endl;
      nwords = 0;
      words[nwords++] = strtok(line," \t\n\r\f");
      indx++;
      //printf("ADSFASDFASDF\n");
      //printf(" indx: %d\n", indx);
      while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;
      
      // Convert (a,b) coordinate to index
      i = atoi(words[0])-1;
      j = atoi(words[2])-1;
      a = atoi(words[1])-1;
      b = atoi(words[3])-1;
      fc = atof(words[4])*13.605698066*(1./0.529177249)*(1./0.529177249); // Convert Ryd/Bohr^2 to eV/A^2;;
      //fc = atof(words[4]);
      
      fc2[indx].i=i;
      fc2[indx].j=j;
      fc2[indx].a=a;
      fc2[indx].b=b;
      fc2[indx].fc=fc;

      printf("ASDFADF\n");
      
    }
  }
  */

  string line;
  ifstream fh_fc2("FC2");
  //string line;
  nfc2 = -1; // Skip first line.
  while (getline(fh_fc2, line))
  {
      nfc2=nfc2+1;
  }

  printf("  Found %d FC2s.\n", nfc2);

  memory->create(fc2,nfc2,"mode:fc2");

  fh_fc2.close();

  ifstream readfile8;
  readfile8.open("FC2");
  string junk;
  readfile8 >> junk;

  for (int n=0; n<nfc2; n++){
    readfile8 >> fc2[n].i >> fc2[n].a >> fc2[n].j >> fc2[n].b >> fc2[n].fc;
  }

  readfile8.close();

  int i,j,k,a,b,c;
  for (int w=0; w<nfc2; w++){
    fc2[w].i = fc2[w].i-1;
    fc2[w].a = fc2[w].a-1;
    fc2[w].j = fc2[w].j-1;
    fc2[w].b = fc2[w].b-1;
    fc2[w].fc = fc2[w].fc*13.605698066*(1./0.529177249)*(1./0.529177249); // Convert Ryd/bohr^2 to eV/A^2
  }

  // Make an FC2 array
  //double fc2_arr[32][32][9][9];
  memory->create(phi,8,8,3,3,"pair:cutsq");
  for (int i=0; i<8; i++){
    for (int j=0; j<8; j++){
      for (int a=0; a<3; a++){
        for (int b=0; b<3; b++){
          phi[i][j][a][b] = 0.0;
        }
      }
    }
  }

  double fc;
  for (int w=0; w<nfc2; w++){

    i = fc2[w].i;
    j = fc2[w].j;
    a = fc2[w].a;
    b = fc2[w].b;
    fc = fc2[w].fc;
    //printf("%d %d %d %d %f\n", i,j,a,b,fc);
    phi[i][j][a][b] = fc;
  }
  
  
}

/* ----------------------------------------------------------------------
   read fc3
------------------------------------------------------------------------- */

void PairTep::read_fc3()
{

  string line;
  ifstream fh_fc3("FC3");
  //string line;
  nfc3 = -1; // Skip first line.
  while (getline(fh_fc3, line))
  {
      nfc3=nfc3+1;
  }

  printf("  Found %d FC3s.\n", nfc3);

  memory->create(fc3,nfc3,"mode:fc3");

  fh_fc3.close();

  ifstream readfile8;
  readfile8.open("FC3");
  string junk;
  readfile8 >> junk;

  for (int n=0; n<nfc3; n++){
    readfile8 >> fc3[n].i >> fc3[n].a >> fc3[n].j >> fc3[n].b >> fc3[n].k >> fc3[n].c >> fc3[n].fc;
  }

  readfile8.close();

  int i,j,k,a,b,c;
  for (int w=0; w<nfc3; w++){
    fc3[w].i = fc3[w].i-1;
    fc3[w].a = fc3[w].a-1;
    fc3[w].j = fc3[w].j-1;
    fc3[w].b = fc3[w].b-1;
    fc3[w].k = fc3[w].k-1;
    fc3[w].c = fc3[w].c-1;
    fc3[w].fc = fc3[w].fc*13.605698066*(1./0.529177249)*(1./0.529177249)*(1./0.529177249); // Convert Ryd/bohr^3 to eV/A^3
  }

  //printf("%d ATOMS ----------------------------------------\n", atom->natoms);
  int natoms = atom->natoms;
  memory->create(psi,3*natoms,3*natoms,3*natoms,"mode:psi");
  //printf("ADSFASDFASDFAFD\n");
  for (int ii=0; ii<3*natoms; ii++){
    //printf("%d\n", ii);
    for (int jj=0; jj<3*natoms; jj++){
      //printf("%d\n", ii);
      for (int kk=0; kk<3*natoms; kk++){
        //printf("%d\n", kk);
        psi[ii][jj][kk] = 0.0;
      }
    }
  }
  //printf("ADSFASDFASDFAFD\n");
  //int i,j,k,a,b,c;
  double fc;
  int ii,jj,kk;
  for (int w=0; w<nfc3; w++){
    i = fc3[w].i;
    j = fc3[w].j;
    k = fc3[w].k;
    a = fc3[w].a;
    b = fc3[w].b;
    c = fc3[w].c;
    fc = fc3[w].fc;
    ii = 3*i+a;
    jj = 3*j+b;
    kk = 3*k+c;
    //printf("%d %d %d\n", ii,jj,kk);
    psi[ii][jj][kk] = fc;
  }
  
  /*
  for (int i=0; i<natoms; i++){
    for (int j=0; j<natoms; j++){
      for (int k=0; k<natoms; k++){
        for (int a=0; a<3; a++){
          for (int b=0; b<3; b++){
            for (int c=0; c<3; c++){

            }
          }
        }
      }
    }
  }
  */

}

/* ----------------------------------------------------------------------
   read fc4
------------------------------------------------------------------------- */

void PairTep::read_fc4()
{

  string line;
  ifstream fh_fc4("FC4");
  //string line;
  nfc4 = -1; // Skip first line.
  while (getline(fh_fc4, line))
  {
      nfc4=nfc4+1;
  }

  printf("  Found %d FC4s.\n", nfc4);

  memory->create(fc4,nfc4,"mode:fc4");

  fh_fc4.close();

  ifstream readfile8;
  readfile8.open("FC4");
  string junk;
  readfile8 >> junk;

  for (int n=0; n<nfc4; n++){
    readfile8 >> fc4[n].i >> fc4[n].a >> fc4[n].j >> fc4[n].b >> fc4[n].k >> fc4[n].c >> fc4[n].l >> fc4[n].d >> fc4[n].fc;
  }

  readfile8.close();

  int i,j,k,a,b,c;
  for (int w=0; w<nfc4; w++){
    fc4[w].i = fc4[w].i-1;
    fc4[w].a = fc4[w].a-1;
    fc4[w].j = fc4[w].j-1;
    fc4[w].b = fc4[w].b-1;
    fc4[w].k = fc4[w].k-1;
    fc4[w].c = fc4[w].c-1;
    fc4[w].l = fc4[w].l-1;
    fc4[w].d = fc4[w].d-1;
    fc4[w].fc = fc4[w].fc*13.605698066*(1./0.529177249)*(1./0.529177249)*(1./0.529177249)*(1./0.529177249); // Convert Ryd/bohr^4 to eV/A^4
  }

  /*
  //printf("%d ATOMS ----------------------------------------\n", atom->natoms);
  int natoms = atom->natoms;
  memory->create(psi,3*natoms,3*natoms,3*natoms,"mode:psi");
  //printf("ADSFASDFASDFAFD\n");
  for (int ii=0; ii<3*natoms; ii++){
    //printf("%d\n", ii);
    for (int jj=0; jj<3*natoms; jj++){
      //printf("%d\n", ii);
      for (int kk=0; kk<3*natoms; kk++){
        //printf("%d\n", kk);
        psi[ii][jj][kk] = 0.0;
      }
    }
  }
  //printf("ADSFASDFASDFAFD\n");
  //int i,j,k,a,b,c;
  double fc;
  int ii,jj,kk;
  for (int w=0; w<nfc3; w++){
    i = fc3[w].i;
    j = fc3[w].j;
    k = fc3[w].k;
    a = fc3[w].a;
    b = fc3[w].b;
    c = fc3[w].c;
    fc = fc3[w].fc;
    ii = 3*i+a;
    jj = 3*j+b;
    kk = 3*k+c;
    //printf("%d %d %d\n", ii,jj,kk);
    psi[ii][jj][kk] = fc;
  }
  */
  /*
  for (int i=0; i<natoms; i++){
    for (int j=0; j<natoms; j++){
      for (int k=0; k<natoms; k++){
        for (int a=0; a<3; a++){
          for (int b=0; b<3; b++){
            for (int c=0; c<3; c++){

            }
          }
        }
      }
    }
  }
  */

}

/* ----------------------------------------------------------------------
   read EQUIL
------------------------------------------------------------------------- */

void PairTep::read_equil()
{

  char line[MAXLINE],*ptr;
  FILE *fp;
  int n,nwords,nuc_poscar;
  char **words = new char*[MAXLINE];
  int info = 0;


  fp = fopen("EQUIL", "r");

  if (fp == NULL) {
    char str[MAXLINE];
    sprintf(str,"Cannot open the equilibrium position file %s", "EQUIL");
    error->all(FLERR,str);
  }

  int eof = 0;
  int atom_count = 0;
  // extract the equilibrium position
  while (!eof) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    //std::cout << line << std::endl;

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    //nwords = atom->count_words(line);
    nwords=3;
    //fprintf(fh_debug, "nwords: %d\n", nwords);
    if (nwords == 0) {
      continue;
    } else {
      nwords = 0;
      words[nwords++] = strtok(line," \t\n\r\f");
      while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;
      x0[atom_count][0] = atof(words[0]);
      x0[atom_count][1] = atof(words[1]);
      x0[atom_count][2] = atof(words[2]);
      //fprintf(fh_debug, "%f %f %f\n", x0[atom_count][0],x0[atom_count][1],x0[atom_count][2]);
      atom_count++;
    }
  }


}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTep::init_one(int i, int j)
{

  /*
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  morse1[i][j] = 2.0*d0[i][j]*alpha[i][j];

  if (offset_flag) {
    double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    offset[i][j] = d0[i][j] * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr));
  } else offset[i][j] = 0.0;

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  morse1[j][i] = morse1[i][j];
  offset[j][i] = offset[i][j];
  */

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTep::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTep::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&d0[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTep::write_restart_settings(FILE *fp)
{
  fwrite(&order,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTep::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&order,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&order,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairTep::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,d0[i][i],alpha[i][i],r0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairTep::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],r0[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairTep::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r,dr,dexp,phi;

  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  dexp = exp(-alpha[itype][jtype] * dr);
  fforce = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;

  phi = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairTep::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  return NULL;
}
