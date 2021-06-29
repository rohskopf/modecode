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
#include "compute_mode_heatflux.h"
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

ComputeModeHeatflux::ComputeModeHeatflux(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_ke(NULL), id_pe(NULL), id_stress(NULL)
{
  if (narg != 6) error->all(FLERR,"Illegal compute mode command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 1;

  // store ke/atom, pe/atom, stress/atom IDs used by heat flux computation
  // insure they are valid for these computations

  int n = strlen(arg[3]) + 1;
  id_ke = new char[n];
  strcpy(id_ke,arg[3]);

  n = strlen(arg[4]) + 1;
  id_pe = new char[n];
  strcpy(id_pe,arg[4]);

  n = strlen(arg[5]) + 1;
  id_stress = new char[n];
  strcpy(id_stress,arg[5]);

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");
  if (strcmp(modify->compute[ike]->style,"ke/atom") != 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute ke/atom");
  if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Compute heat/flux compute ID does not compute pe/atom");
  if (modify->compute[istress]->pressatomflag == 0)
    error->all(FLERR,
               "Compute heat/flux compute ID does not compute stress/atom");

  vector = new double[6];
  iflux_arr = new double[1];

  //rank = MPI_COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  fh_disp = fopen("disp.dat", "w");
  fh_xm = fopen("xm.dat", "w");
  fh_fm = fopen("fm.dat", "w");
  fh_tm = fopen("tm.dat", "w");
  fh_miflux = fopen("miflux.dat", "w");
  fh_iflux = fopen("iflux.dat", "w");
  fh_fv = fopen("fv.dat", "w");
  fh_fvtot = fopen("ipt.dat", "w");
  fh_em = fopen("em.dat", "w");
  fh_qtep = fopen("q_tep.dat", "w");
  fh_etep = fopen("tep_energy.dat", "w");

  fh_debug = fopen("D_MODE_HEATFLUX", "w");

}

/* ---------------------------------------------------------------------- */

ComputeModeHeatflux::~ComputeModeHeatflux()
{
  delete [] id_ke;
  delete [] id_pe;
  delete [] id_stress;
  delete [] vector;
  delete [] iflux_arr;

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

  memory->destroy(sidemap);
  memory->destroy(regions);
  memory->destroy(zvals);

  fclose(fh_disp);
  fclose(fh_xm);
  fclose(fh_fm);
  fclose(fh_miflux);
  fclose(fh_iflux);
  fclose(fh_fv);
  fclose(fh_fvtot);
  fclose(fh_em);
  fclose(fh_qtep);
  fclose(fh_etep);


  fclose(fh_debug);

}

/* ---------------------------------------------------------------------- */

void ComputeModeHeatflux::init()
{

  // error checks

  int ike = modify->find_compute(id_ke);
  int ipe = modify->find_compute(id_pe);
  int istress = modify->find_compute(id_stress);
  if (ike < 0 || ipe < 0 || istress < 0)
    error->all(FLERR,"Could not find compute heat/flux compute ID");

  c_ke = modify->compute[ike];
  c_pe = modify->compute[ipe];
  c_stress = modify->compute[istress];

  int natoms = atom->natoms;

  if (universe->me == 0) {
    printf("--------------------------- %d -----------------------------\n", universe->me);
  }

  printf("    -------------- Initializing Mode Computations on Proc %d. -------------------\n", universe->me);
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

  readfile2.open("EMAT");
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

  /* SIDEMAP */

  ifstream readfile4;

  readfile4.open("SIDEMAP");
  //readfile2.open("ev_real.txt");

  if (readfile4.is_open()) {

    int ntypes = atom->ntypes;
    printf("    %d atom types.\n", ntypes);
    memory->create(sidemap,ntypes,"mode:sidemap");

    for (int t=0; t<ntypes; t++){
      readfile4>>sidemap[t];
      //printf(" %d\n", sidemap[t]);

    }

  }

  /* REGIONS */

  ifstream readfile5;

  readfile5.open("REGIONS");

  if (readfile5.is_open()) {
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
  }

  /* MCC2 */

  ifstream fh_mcc2("MCC2");
  string line;
  nmcc2 = 0;
  while (getline(fh_mcc2, line))
  {
      nmcc2=nmcc2+1;
  }

  printf("  Found %d MCC2s.\n", nmcc2);

  memory->create(mcc2,nmcc2,"mode:mcc2");

  fh_mcc2.close();

  ifstream readfile7;
  readfile7.open("MCC2");

  for (int n=0; n<nmcc2; n++){
    readfile7 >> mcc2[n].i >> mcc2[n].j >> mcc2[n].val;
  }

  readfile7.close();
  

  /* MCC3 */

  
  ifstream fh_mcc3("MCC3");
  //string line;
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


  /* FC2*/

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
    readfile8 >> fc2[n].i >> fc2[n].a >> fc2[n].j >> fc2[n].b >> fc2[n].val;
  }

  readfile8.close();

  int i,j,a,b;
  for (int w=0; w<nfc2; w++){
    fc2[w].i = fc2[w].i-1;
    fc2[w].a = fc2[w].a-1;
    fc2[w].j = fc2[w].j-1;
    fc2[w].b = fc2[w].b-1;
    //printf("%d %d %d %d %e\n", i,a,j,b,fc2[w].val);
    fc2[w].val = fc2[w].val*13.605698066*(1./0.529177249)*(1./0.529177249); // Convert Ryd/bohr^2 to eV/A^2
    //printf("%d %d %d %d %e\n", i,a,j,b,fc2[w].val);
  }
  

  /* FC3*/

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

  ifstream readfile9;
  readfile9.open("FC3");
  junk;
  readfile9 >> junk;

  for (int n=0; n<nfc3; n++){
    readfile9 >> fc3[n].i >> fc3[n].a >> fc3[n].j >> fc3[n].b >> fc3[n].k >> fc3[n].c >> fc3[n].val;
  }

  readfile9.close();

  int k,c;
  for (int w=0; w<nfc3; w++){
    fc3[w].i = fc3[w].i-1;
    fc3[w].a = fc3[w].a-1;
    fc3[w].j = fc3[w].j-1;
    fc3[w].b = fc3[w].b-1;
    fc3[w].k = fc3[w].k-1;
    fc3[w].c = fc3[w].c-1;
    //printf("%d %d %d %d %e\n", i,a,j,b,fc2[w].val);
    fc3[w].val = fc3[w].val*13.605698066*(1./0.529177249)*(1./0.529177249)*(1./0.529177249); // Convert Ryd/bohr^3 to eV/A^3
    //printf("%d %d %d %d %e\n", i,a,j,b,fc2[w].val);
  }

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

void ComputeModeHeatflux::compute_vector()
{

  invoked_vector = update->ntimestep;

  // invoke 3 computes if they haven't been already

  if (!(c_ke->invoked_flag & INVOKED_PERATOM)) {
    c_ke->compute_peratom();
    c_ke->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
    c_pe->compute_peratom();
    c_pe->invoked_flag |= INVOKED_PERATOM;
  }
  if (!(c_stress->invoked_flag & INVOKED_PERATOM)) {
    c_stress->compute_peratom();
    c_stress->invoked_flag |= INVOKED_PERATOM;
  }

  // heat flux vector = jc[3] + jv[3]
  // jc[3] = convective portion of heat flux = sum_i (ke_i + pe_i) v_i[3]
  // jv[3] = virial portion of heat flux = sum_i (stress_tensor_i . v_i[3])
  // normalization by volume is not included

  double *ke = c_ke->vector_atom;
  double *pe = c_pe->vector_atom;
  double **stress = c_stress->array_atom;

  double **x = atom->x;
  double **v = atom->v;
  /*
  for (int i=0; i<32; i++){
    fprintf(fh_debug, "%f %f %f\n", v[i][0],v[i][1],v[i][2]);

  }
  */
  double **f = atom->f;
  int *tag = atom->tag;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int natoms = atom->natoms;

  double jc[3] = {0.0,0.0,0.0};
  double jv[3] = {0.0,0.0,0.0};
  double eng;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      eng = pe[i] + ke[i];
      jc[0] += eng*v[i][0];
      jc[1] += eng*v[i][1];
      jc[2] += eng*v[i][2];
      jv[0] -= stress[i][0]*v[i][0] + stress[i][3]*v[i][1] +
        stress[i][4]*v[i][2];
      jv[1] -= stress[i][3]*v[i][0] + stress[i][1]*v[i][1] +
        stress[i][5]*v[i][2];
      jv[2] -= stress[i][4]*v[i][0] + stress[i][5]*v[i][1] +
        stress[i][2]*v[i][2];
    }
  }

  // convert jv from stress*volume to energy units via nktv2p factor

  double nktv2p = force->nktv2p;
  jv[0] /= nktv2p;
  jv[1] /= nktv2p;
  jv[2] /= nktv2p;

  // sum across all procs
  // 1st 3 terms are total heat flux
  // 2nd 3 terms are just conductive portion

  double data[6] = {jc[0]+jv[0],jc[1]+jv[1],jc[2]+jv[2],jc[0],jc[1],jc[2]};
  MPI_Allreduce(data,vector,6,MPI_DOUBLE,MPI_SUM,world);


  // Displace all atoms along n=3.
  /*
  for (int i=0; i<natoms; i++){
    for (int a=0; a<3; a++){
      x[i][a] += emat[3*i+a][3];
    }
    //fprintf(fh_disp, "%f %f %f\n", x[i][0],x[i][1],x[i][2]);
  }
  */

  // Compute atomic displacements.
  //fprintf(fh_disp, "\n");-
  //fprintf(fh_disp, "Timestep: %d\n", update->ntimestep);
  //fprintf(fh_disp, "\n");
  //printf("x0[%d]: %f %f %f\n", 17-1, x0[17-1][0], x0[17-1][1], x0[17-1][2]);
  int iindx;
  for (int i=0; i<nlocal; i++){
    iindx = tag[i]-1;
    //printf("%d %d %d\n", i, tag[i]-1,universe->me);
    //printf("%d %d %f %f %f\n", i, iindx, x[i][0], x[i][1],x[i][2]);
    if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
      //u_p[i]=0.0;
      //u[i] = 0.0; // zero the total coordinates
      //printf("x[%d]: %f %f %f\n", i, x[i][0],x[i][1],x[i][2]);
      //printf("x0[%d]: %f %f %f\n", tag[i]-1, x0[tag[i]-1][0], x0[tag[i]-1][1], x0[tag[i]-1][2]);
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
 

  // Convert back to cartesian coordinates just to check.
  /*
  double uia;
  double via;
  double mi;
  double u2[32][3],v2[32][3]; // u and v from mode caluclations
  //fprintf(fh_debug,"---------------------- TIMESTEP %d\n", update->ntimestep);
  for (int i=0; i<natoms; i++){
  //for (int i=0; i<1; i++){
    mi = mass[type[i]]*(1.0/(6.02214076e23*1e3)); // atom mass in kg
    for (int a=0; a<3; a++){
    //for (int a=2; a<3; a++){
      uia = 0.0;
      via = 0.0;
      v2[i][a] = 0.0;
      u2[i][a] = 0.0;
      for (int n=0; n<3*natoms; n++){
        via += (1.0/sqrt(mi))*vm[n]*emat[3*i+a][n];
        v2[i][a] += (1.0/sqrt(mi))*vm[n]*emat[3*i+a][n];
        uia += (1.0/sqrt(mi))*xm[n]*emat[3*i+a][n];
        //fprintf(fh_disp, "%d %d %d %e\n", i,a,n,emat[3*i+a][n]);
      }
      //fprintf(fh_debug, "%d %d %e %e\n", i,a,v[i][a],v2[i][a]);
     
      
      if ( abs(u[i][a]-uia)>1e-6 || abs(v[i][a]-via)>1e-6){
        //fprintf(fh_debug, "%d %d %e %e\n", i,a,u_p[i][a] - uia,v[i][a] - v2[i][a]);
        fprintf(fh_disp,"TIMESTEP %d\n", update->ntimestep);
        fprintf(fh_disp, "  %e %e\n", u[i][a],uia);
      }
      

    }
      
  }
  */

  // Interfacial flux
  /*
  int side_i; // 1 - left side of interface
              // 2 - right side of interface
  int side_j; // 1 - left side of interface
              // 2 - right side of interface
  iflux_p=0.0;
  iflux_arr[0]=0.0; // this is the MPI reduced iflux
  for (int n=0; n<3*natoms; n++){
    miflux_p[n]=0.0;
    miflux[n]=0.0;
  }

  int regi; // region i
  int regj; // region j
  */
  /*
  for (int i=0; i<nlocal; i++){
    if (mask[i] & groupbit){ // This keeps the atoms in a particular group.
      regi=regions[i];
      //printf(" type, side: %d %d\n", type[i],side);
      if (regi == 2){
        for (int j=0; j<nlocal; j++){
          if (mask[j] & groupbit){ // This keeps the atoms in a particular group.
            regj = regions[j];
            if (regj==3){
              double mj = mass[type[j]]; //mass[type[j]]*(1.0/6.0221409e+23)*1e-3;
              for (int n=0; n<3*natoms; n++){
                for (int a=0; a<3; a++){
                  iflux_p -= 0.5*(1.0/sqrt(mj))*fmi[natoms*n+i]*emat[3*j+a][n]*(v[i][a]+v[j][a]);
                  miflux_p[n] -= 0.5*(1.0/sqrt(mj))*fmi[natoms*n+i]*emat[3*j+a][n]*(v[i][a]+v[j][a]);
                }
              }
            }
          }
        }
      }
    }
  }
  */
  /*
  MPI_Allreduce(miflux_p,miflux,3*natoms,MPI_DOUBLE,MPI_SUM,world);
  double data_iflux[1] = {iflux_p};
  MPI_Allreduce(data_iflux,iflux_arr,1,MPI_DOUBLE,MPI_SUM,world);
  iflux = iflux_arr[0];
  double iflux2 = 0.0;
  for (int n=0;n<3*natoms;n++){
    iflux2 += miflux[n];
  }
  */

  double iflux_test = 0.0;
  double sqmr; //sqrt mass ratio
  double edotf;
  double edotv;
  double edotfi,edotfj;
  double psij[3],psji[3];
  double vdotp_ij;
  double vdotp_ji;

  /*
  for (int n=0; n<3*natoms; n++){
    for (int i=0; i<natoms; i++){
      for (int j=0; j<natoms; j++){
        regi=regions[i];
        regj=regions[j];
        if (regi==2 && regj==3){
          sqmr = sqrt(mass[type[i]]/mass[type[j]]);
          edotfi = emat[3*i+0][n]*f[i][0] + emat[3*i+1][n]*f[i][1] + emat[3*i+2][n]*f[i][2];
          edotfj = emat[3*j+0][n]*f[j][0] + emat[3*j+1][n]*f[j][1] + emat[3*j+2][n]*f[j][2];

          psij[0] = emat[3*j+0][n]*edotfi;
          psij[1] = emat[3*j+1][n]*edotfi;
          psij[2] = emat[3*j+2][n]*edotfi;

          psji[0] = emat[3*i+0][n]*edotfj;
          psji[1] = emat[3*i+1][n]*edotfj;
          psji[2] = emat[3*i+2][n]*edotfj;

          vdotp_ij = v[i][0]*psij[0] + v[i][1]*psij[1] + v[i][2]*psij[2];
          vdotp_ji = v[j][0]*psji[0] + v[j][1]*psji[1] + v[j][2]*psji[2];

          iflux_test -= sqmr*vdotp_ij - (1.0/sqmr)*vdotp_ji;

        } // if (regi==2 && regj==3){
      }
    }
  }
  */


  // Mode temperatures and energies
  double pi = 3.1415926535;
  for (int n=3; n<3*natoms; n++){
    //temp_mode[n] = ( xm[n]*xm[n]*input->mcc2[n].val )/(1.38064852e-23);
    //tm[n] = 0.5*(vm[n]*vm[n])/(1.38064852e-23);
    tm[n] = 0.5*(vm[n]*vm[n]*100*100)/(1.38064852e-23); // factors of 100 since 100 m/s = 1 A/ps
    em[n] = (0.5*vm[n]*vm[n]*100*100 + 0.5*freq[n]*freq[n]*2.0*3.1415926535*2.0*3.1415926535*1e12*1e12*xm[n]*xm[n]*1e-10*1e-10)*6.242e+18; // eV, since 6.242e+18 eV = 1 J
    //if (n==2209){
    //    printf(" vm[n]: %e\n", vm[n]);
    //    printf(" temp_mode[n]: %e\n", temp_mode[n]);
    //}
    //kemi[n] = 0.5*vm[n]*vm[n];
    //emi[n] = (kemi[n] + pemi[n])*6.242e+18;
  }

  // Print stuff
  
  //fprintf(fh_xm, "Timestep: %d\n", update->ntimestep);
  
  
  //fprintf(fh_xm, "\n");
  //fprintf(fh_xm, "\n");

  //fprintf(fh_fm, "\n");
  //fprintf(fh_fm, "\n");

  //fprintf(fh_tm, "\n");
  //fprintf(fh_tm, "\n");
  

  /*
  fprintf(fh_em, "\n");
  fprintf(fh_em, "\n");  
  
  for (int n=3; n<3*natoms; n++){
    //fprintf(fh_xm, "%d %e %e\n", n,freq[n],xm[n]);
    //fprintf(fh_fm, "%d %e %e\n", n,freq[n],fm[n]);
    //fprintf(fh_tm, "%d %e %e\n", n,freq[n],tm[n]);
    fprintf(fh_em, "%d %e %e\n", n,freq[n],em[n]);
  }
  */

  // SINGLE PLOT DAT FILES
  
  //fprintf(fh_xm, "%.2f %e\n",  0.5*update->ntimestep*1e-3, xm[1186]);
  //fprintf(fh_tm, "%.2f %e\n",  0.5*update->ntimestep*1e-3, tm[1906]);
  //fprintf(fh_fm, "%.2f %e\n",  0.5*update->ntimestep*1e-3, fm[605]);
  //fprintf(fh_em, "%.2f %e %e\n",  0.5*update->ntimestep*1e-3, em[605], em[1906]);
  //double fvtot = fm[931]*vm[931]*1e-3; // Originally units are (eV/A)*(A/ps) = eV/ps. Multiply by 1e-3 to get eV/fs
  //fprintf(fh_fvtot, "%.2f %e\n",  0.5*update->ntimestep*1e-3, fvtot);


  // TEP heat transfer
  /*
  int i,j,k;
  int a,b,c;
  double fc;
  double ht = 0.0;
  double ft[16][3];
  for (int n=0; n<natoms; n++){
    ft[n][0] = 0.0;
    ft[n][1] = 0.0;
    ft[n][2] = 0.0;
  }
  for (int s=0; s<nfc2; s++){
    i = fc2[s].i-1;
    j = fc2[s].j-1;
    a = fc2[s].a-1;
    b = fc2[s].b-1;
    fc = fc2[s].val; // eV/A^2
    ft[i][a] -= fc*u_p[j][b];
  }

  for (int i=0; i<natoms; i++){
    //fprintf(fh_qtep, "%f %f %f\n", ft[i][0],ft[i][1],ft[i][2]);
    fprintf(fh_qtep, "u[%d]: %f %f %f\n", i, u_p[i][0],u_p[i][1],u_p[i][2]);
    fprintf(fh_qtep, "ft[%d]: %f %f %f\n", i, ft[i][0],ft[i][1],ft[i][2]);
  }
  //fprintf(fh_qtep, "%f %e\n", 0.5*update->ntimestep*1e-3, ht);
  */

  /*
  Interface modes:
  (indices starting at 0)
  73
  74
  76
  77
  */

  // Mode interaction heat transfer
  if (universe->me == 0){

    double ht = 0.0;
    int n1,n2,n3;
    double xn1,xn2,xn3;
    double vn1,vn2,vn3;
    double mcc;
    for (int s=0; s<nmcc2; s++){
      n1 = mcc2[s].i;
      n2 = mcc2[s].j;
      xn1 = xm[n1]*1e-10; // sqrt(kg) m since 1e-10 m = 1 A.
      xn2 = xm[n2]*1e-10; // sqrt(kg) m since 1e-10 m = 1 A.
      vn1 = vm[n1]*100; // sqrt(kg) m/s since 100 m/s = 1 A/ps.
      vn2 = vm[n2]*100; // sqrt(kg) m/s since 100 m/s = 1 A/ps.
      mcc = mcc2[s].val; // J/(kg*m^2)
      // Multiply by 1.0 when comparing to TEP.
      //ht += 1.0*mcc*xn2*vn1*6.242e+6; // eV/ps since 1 J/s = 6.242e+6 eV/ps
      // Multiple by 2.0 when comparing to others?
      if (abs(mcc)>1e24){
      //if (n1!=73 && n1!=74 && n1!=76 && n1!=77 && n2!=73 && n2!=74 && n2!=76 && n2!=77){
        ht += 1.0*mcc*xn2*vn1*6.242e+6; // eV/ps since 1 J/s = 6.242e+6 eV/ps
      //}
      }

      if (update->ntimestep==200){
        //fprintf(fh_debug, "%d %d %e %e %e\n", n1,n2,mcc,xm[n2],vm[n1]);
      }

      /*
      if (abs(mcc)>1e24){
        ht += mcc*xn2*vn1*6.241509e18*1e-12;
      }
      */

    }

    
    for (int s=0; s<nmcc3; s++){
      n1 = mcc3[s].i;
      n2 = mcc3[s].j;
      n3 = mcc3[s].k;
      xn1 = xm[n1]*1e-10; // sqrt(kg) m since 1e-10 m = 1 A.
      xn2 = xm[n2]*1e-10; // sqrt(kg) m since 1e-10 m = 1 A.
      xn3 = xm[n3]*1e-10; // sqrt(kg) m since 1e-10 m = 1 A.
      vn1 = vm[n1]*100; // sqrt(kg) m/s since 100 m/s = 1 A/ps.
      vn2 = vm[n2]*100; // sqrt(kg) m/s since 100 m/s = 1 A/ps.
      vn3 = vm[n3]*100; // sqrt(kg) m/s since 100 m/s = 1 A/ps.
      mcc = mcc3[s].val; // J/(kg^(3/2)*m^3
      if (abs(mcc)>1e46){
      //if (n1!=73 && n1!=74 && n1!=76 && n1!=77 && n2!=73 && n2!=74 && n2!=76 && n2!=77 && n3!=73 && n3!=74 && n3!=76 && n3!=77){
        ht += (1.0/3.0)*mcc*xn3*xn2*vn1*6.242e+6; // eV/ps since 1 J/s = 6.242e+6 eV/ps
      //}
      }
    }
    

    fprintf(fh_fvtot, "%f %e\n", 0.5*update->ntimestep*1e-3, ht);

  }
  // TEP energy
  /*
  double fc;
  double uia,uib,ujb;
  double pe_tot,pe_a,pe_b;
  pe_tot = 0.0;
  pe_a = 0.0;
  pe_b = 0.0;

  int i,j,k,a,b,c;
  for (int w=0; w<nfc2; w++){

    //printf("%d\n", w);
    i = fc2[w].i;
    j = fc2[w].j;
    a = fc2[w].a;
    b = fc2[w].b;
    fc = fc2[w].val;
    //printf("  %d %d %d %d.\n",i,a,j,b);

    uia = u[i][a];
    ujb = u[j][b];

    //f2[i][a] -= fc*ujb;

    // TEP potential energy
    pe_tot += 0.5*fc*uia*ujb;
    if (type[i]==1){
      pe_a += 0.5*fc*uia*ujb;
    }
    else if (type[i]==2){
      pe_b += 0.5*fc*uia*ujb;
    }
    else{
      printf("FUCK!\n");
    }
    
  }
  */

  /*
  double ukc;
  for (int w=0; w<nfc3; w++){

    //printf("%d\n", w);
    i = fc3[w].i;
    j = fc3[w].j;
    k = fc3[w].k;
    a = fc3[w].a;
    b = fc3[w].b;
    c = fc3[w].c;
    fc = fc3[w].val;
    //printf("  %d %d %d %d.\n",i,a,j,b);

    uia = u[i][a];
    ujb = u[j][b];
    ukc = u[k][c];

    //f2[i][a] -= fc*ujb;

    // TEP potential energy
    pe_tot += (1.0/6.0)*fc*uia*ujb*ukc;
    if (type[i]==1){
      pe_a += (1.0/6.0)*fc*uia*ujb*ukc;
    }
    else if (type[i]==2){
      pe_b += (1.0/6.0)*fc*uia*ujb*ukc;
    }
    else{
      //printf("FUCK!\n");
    }
    
  }
  */
  //fprintf(fh_etep, "%.5f %f %f %f\n", 0.5*update->ntimestep*1e-3,pe_a,pe_b,pe_tot);

  /*
  for (int i=0; i<natoms; i++){
    printf("%e %e %e\n", f2[i][0],f2[i][1],f2[i][3]);
  }
  */
  /*
  double fc;
  for (int i=0; i<natoms; i++){
    for (int j=0; j<natoms; j++){
      for (int a=0; a<3; a++){
        for (int b=0; b<3; b++){

          uia = x[i][a] - x0[i][a];
          uib = x[i][b] - x0[i][b];
          ujb = x[j][b] - x0[j][b];

          if (std::abs(uia) > boxhi_x/2.0 && uia > 0) uia -= boxhi_x;
          else if (std::abs(uia) > boxhi_x/2.0 && uia < 0) uia += boxhi_x;

          if (std::abs(uib) > boxhi_x/2.0 && uib > 0) uib -= boxhi_x;
          else if (std::abs(uib) > boxhi_x/2.0 && uib < 0) uib += boxhi_x;

          if (std::abs(ujb) > boxhi_x/2.0 && ujb > 0) ujb -= boxhi_x;
          else if (std::abs(ujb) > boxhi_x/2.0 && ujb < 0) ujb += boxhi_x;

          // TEP energy
          pe += 0.5*fc*uia*ujb;
          if (itype==1){
            pe_a += 0.5*fc*uia*ujb;
          }
          else if (itype==2){
            pe_b += 0.5*fc*uia*ujb;
          }


        }

      }
    }

  }
  */

  // Anharmonic forces and power transfer.
  /*
  double fv;
  int i,j,k;
  fprintf(fh_fv, " %0.2f ", 0.5*update->ntimestep*1e-3);
  printf(" %d\n", nmcc3);
  int counter = 0;
  //for (int n=0; n<nmcc3; n++){
  for (int n=0; n<nmcc3; n++){
    i=mcc3[n].i;
    j=mcc3[n].j;
    k=mcc3[n].k;
    //printf("%e\n", mcc3[n].val);
    //if (abs(mcc3[n].val)>1e48){
    //printf("ASDF\n");
    fv = -1.0*mcc3[n].val*xm[j]*xm[k]*1e-10*1e-10*vm[i]*100*6.242e+6; // eV/ps
    fprintf(fh_fv, "%e ", fv);
    counter++;
    //}

  }
  fprintf(fh_fv, "\n");

  if (update->ntimestep==0){
    printf(" ***** Found %d MCC3s above bound *****\n", counter);
  }
  */
  // For Mode 605.
  /*
  double f1,f2,f3;
  double fv1,fv2,fv3;
  f1 = -1.0*mcc3[0]*xm[1186]*xm[1186]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f2 = -1.0*mcc3[1]*xm[1186]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f3 = -1.0*mcc3[2]*xm[1906]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  fv1 = f1*vm[605]*100*6.242e+6; // N/sqrt(kg) * sqrt(kg)*m/s = N*m/s. Multiply by 6.242e+6 to convert to eV/ps, since 1 J/s = 6.242e+6 eV/ps. Factor of 100 converts A/ps to m/s for the velocity.
  fv2 = f2*vm[605]*100*6.242e+6; // eV/ps
  fv3 = f3*vm[605]*100*6.242e+6; // eV/ps
  fprintf(fh_fv, " %0.2f %e %e %e\n", 0.5*update->ntimestep*1e-3, fv1,fv2,fv3);
  */

  /*
  1186 604 604 -4.29913242442938561606e+35
  1186 604 605 -8.61022846845004687996e+35
  1186 604 1186 -4.42875156879688460168e+48
  1186 604 1906 -6.66703612060836330866e+47
  1186 605 604 -8.61012071814902748046e+35
  1186 605 605 -2.71580837403027022113e+33
  1186 605 1186 5.77235791965297911004e+48
  1186 605 1906 -1.95764722287588829226e+48
  1186 1186 604 -4.42875156879688330360e+48
  1186 1186 605 5.77235791965297911004e+48
  1186 1186 1186 -6.07988888457994467499e+37
  1186 1186 1906 1.68620102361842350857e+38
  1186 1906 604 -6.66703612060838764755e+47
  1186 1906 605 -1.95764722287588764323e+48
  1186 1906 1186 1.68612289128063035904e+38
  1186 1906 1906 -4.26931440218365803418e+37
  */
  /*
  double f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16;
  double fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8,fv9,fv10,fv11,fv12,fv13,fv14,fv15,fv16;
  f1 = -1.0*mcc3[0]*xm[604]*xm[604]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f2 = -1.0*mcc3[1]*xm[604]*xm[605]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f3 = -1.0*mcc3[2]*xm[604]*xm[1186]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f4 = -1.0*mcc3[3]*xm[604]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f5 = -1.0*mcc3[4]*xm[605]*xm[604]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f6 = -1.0*mcc3[5]*xm[605]*xm[605]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f7 = -1.0*mcc3[6]*xm[605]*xm[1186]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f8 = -1.0*mcc3[7]*xm[605]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f9 = -1.0*mcc3[8]*xm[1186]*xm[604]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f10 = -1.0*mcc3[9]*xm[1186]*xm[605]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f11 = -1.0*mcc3[10]*xm[1186]*xm[1186]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f12 = -1.0*mcc3[11]*xm[1186]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f13 = -1.0*mcc3[12]*xm[1906]*xm[604]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f14 = -1.0*mcc3[13]*xm[1906]*xm[605]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f15 = -1.0*mcc3[14]*xm[1906]*xm[1186]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  f16 = -1.0*mcc3[15]*xm[1906]*xm[1906]*1e-10*1e-10; // (1/sqrt(kg^3))*J/m * (sqrt(kg)*m)^2 = (1/sqrt(kg))*N
  fv1 = f1*vm[1186]*100*6.242e+6; // N/sqrt(kg) * sqrt(kg)*m/s = N*m/s. Multiply by 6.242e+6 to convert to eV/ps, since 1 J/s = 6.242e+6 eV/ps. Factor of 100 converts A/ps to m/s for the velocity.
  fv2 = f2*vm[1186]*100*6.242e+6; // eV/ps
  fv3 = f3*vm[1186]*100*6.242e+6; // eV/ps
  fv4 = f4*vm[1186]*100*6.242e+6; // eV/ps
  fv5 = f5*vm[1186]*100*6.242e+6; // eV/ps
  fv6 = f6*vm[1186]*100*6.242e+6; // eV/ps
  fv7 = f7*vm[1186]*100*6.242e+6; // eV/ps
  fv8 = f8*vm[1186]*100*6.242e+6; // eV/ps
  fv9 = f9*vm[1186]*100*6.242e+6; // eV/ps
  fv10 = f10*vm[1186]*100*6.242e+6; // eV/ps
  fv11 = f11*vm[1186]*100*6.242e+6; // eV/ps
  fv12 = f12*vm[1186]*100*6.242e+6; // eV/ps
  fv13 = f13*vm[1186]*100*6.242e+6; // eV/ps
  fv14 = f14*vm[1186]*100*6.242e+6; // eV/ps
  fv15 = f15*vm[1186]*100*6.242e+6; // eV/ps
  fv16 = f16*vm[1186]*100*6.242e+6; // eV/ps
  fprintf(fh_fv, " %0.2f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 0.5*update->ntimestep*1e-3,fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8,fv9,fv10,fv11,fv12,fv13,fv14,fv15,fv16);
  */
  
  


  //fprintf(fh_tm, "%.2f %e %e %e %e %e\n", 0.5*update->ntimestep, tm[1498], tm[1500], tm[1504], tm[1506], tm[1508]); 
  

  /*
  fprintf(fh_miflux, "\n");
  //fprintf(fh_xm, "Timestep: %d\n", update->ntimestep);
  fprintf(fh_miflux, "\n");
  for (int n=3; n<3*natoms; n++){
    fprintf(fh_miflux, "%d %e %e\n", n+1,freq[n],miflux[n]);
    //fprintf(fh_xm, "%d %e\n", n+1,xm[n]);
  }
  */

  //fprintf(fh_iflux, "\n");
  //fprintf(fh_xm, "Timestep: %d\n", update->ntimestep);
  //fprintf(fh_iflux, "\n");
  //fprintf(fh_iflux, "%e %e %e\n", iflux, iflux2, iflux_test);


}
