/*
 mode.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.

 The Mode class stores all per-mode arrays, and information pertaining to modes.
 This class also stores per-atom arrays (position, velocities, etc.).
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "mpi.h"
#include <math.h>       /* sqrt */
#include <random>

#include "mode.h"
#include "memory.h"
#include "input.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

Mode::Mode(EM3 *em3) : Pointers(em3) {


}

Mode::~Mode() 
{
    memory->deallocate(xm);
    memory->deallocate(xmta);
    memory->deallocate(vm);
    memory->deallocate(fm);
    memory->deallocate(fm2);
    memory->deallocate(fm3);
    memory->deallocate(fm4);
    memory->deallocate(em);

    memory->deallocate(xa);
    memory->deallocate(va);
    memory->deallocate(fa);
    memory->deallocate(fa2);
    memory->deallocate(fa3);
    memory->deallocate(fa4);

    memory->deallocate(fmta);
    memory->deallocate(fmta2);
    memory->deallocate(fmta3);
    memory->deallocate(fmta4);


    fclose(fh_em);
    fclose(fh_fv);
    

};

void Mode::initialize()
{

    natoms = input->natoms;
    nmodes = 3*natoms;
    //printf("%d\n", natoms);
    memory->allocate(xm,nmodes);
    memory->allocate(xmta,natoms,3);
    memory->allocate(vm,nmodes);
    memory->allocate(fm,nmodes);
    memory->allocate(fm2,nmodes);
    memory->allocate(fm3,nmodes);
    memory->allocate(fm4,nmodes);
    memory->allocate(em, nmodes);

    memory->allocate(xa,natoms,3);
    memory->allocate(va,natoms,3);
    memory->allocate(fa,natoms,3);
    memory->allocate(fa2,natoms,3);
    memory->allocate(fa3,natoms,3);
    memory->allocate(fa4,natoms,3);

    memory->allocate(fmta,natoms,3);
    memory->allocate(fmta2,natoms,3);
    memory->allocate(fmta3,natoms,3);
    memory->allocate(fmta4,natoms,3);

    for (int i=0; i<nmodes; i++){
        xm[i] = input->xm[i];
        //printf(" %e\n", q[i]);
        vm[i] = input->vm[i];
        fm[i] = 0.0; // Initialize mode forces to be zero... they may change when potential is calculated.
        fm2[i]=0.0;
        fm3[i]=0.0;
        fm4[i]=0.0;
    }

    //mode2Atom(); // Convert mode coordinates to atom positions

    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            xa[i][a] = input->xa[i][a];
            va[i][a] = input->va[i][a];
            fa[i][a] = 0.0; // Initialize atom forces to be zero... they may change when potential is calculated.
            fa2[i][a]=0.0;
            fa3[i][a]=0.0;
            fa4[i][a]=0.0;
            fmta[i][a] = 0.0;
            fmta2[i][a] = 0.0;
            fmta3[i][a] = 0.0;
            fmta4[i][a] = 0.0;
        }
    }



    fh_em = fopen("em.dat", "w");
    fh_fv = fopen("fv.dat", "w");

}

/*
Convert mode quantities to per-atom quantities
*/
void Mode::mode2Atom()
{
    
    
    double **xa0 = input->xa0;
    double *mass = input->mass;
    int *type = input->type;
    double **emat = input->emat;
    for (int i=0; i<natoms; i++){
        xmta[i][0] = xa0[i][0];
        xmta[i][1] = xa0[i][1];
        xmta[i][2] = xa0[i][2];
        for (int n=3; n<nmodes; n++){
            xmta[i][0] += (1.0/sqrt(mass[type[i]-1]))*emat[0+i*3][n]*xm[n];
            xmta[i][1] += (1.0/sqrt(mass[type[i]-1]))*emat[1+i*3][n]*xm[n];
            xmta[i][2] += (1.0/sqrt(mass[type[i]-1]))*emat[2+i*3][n]*xm[n];
        }
    }
    
    

    for (int i=0; i<natoms; i++){
        for (int a=0; a<3; a++){
            fmta[i][a]=0.0;
            fmta2[i][a]=0.0;
            fmta3[i][a]=0.0;
            fmta4[i][a]=0.0;
            for (int n=3; n<nmodes; n++){
                fmta[i][a] += fm[n]*sqrt(mass[type[i]-1])*input->emat[3*i+a][n];
                fmta2[i][a] += fm2[n]*sqrt(mass[type[i]-1])*input->emat[3*i+a][n];
                fmta3[i][a] += fm3[n]*sqrt(mass[type[i]-1])*input->emat[3*i+a][n];
                fmta4[i][a] += fm4[n]*sqrt(mass[type[i]-1])*input->emat[3*i+a][n];
            }
        }
    }   
    
    /*
    for (int n=0; n<nmodes; n++){
        fm2a[n] = 0.0;
        for (int i=0; i<natoms; i++){
            for (int a=0; a<3; a++){
                fm2a[n] += (fa[i][a]*input->emat[3*i+a][n])/(sqrt(input->mass));
            }
        }
    }
    */
}

/*
Convert atom quantities to mode quantities
*/
void Mode::atom2Mode()
{

    double **xa0 = input->xa0;
    double *mass = input->mass;
    int *type = input->type;
    double **emat = input->emat;

    for (int n=3; n<nmodes; n++){
        xm[n]=0.0;
        for (int i=0; i<natoms; i++){
            for (int a=0; a<3; a++){
                xm[n] += sqrt(mass[type[i]-1])*emat[3*i+a][n]*(xa[i][a]-xa0[i][a]);
            }
        }
    }

}

/*
Calculate mode energies.
*/
void Mode::calcModeEnergy()
{

  for (int n=3; n<nmodes; n++){
    em[n] = (0.5*input->mcc2[n].val*xm[n]*xm[n] + 0.5*vm[n]*vm[n])*6.242e+18; // eV
  }

  int modeindx = 10;
  int i,j,k,l;
  int n,m,o,p;
  double mcc;
  fv = 0.0;
  bool thismode;
  for (int w=0; w<input->nmcc3; w++){
      i=input->mcc3[w].i;
      j=input->mcc3[w].j;
      k=input->mcc3[w].k;
      mcc = input->mcc3[w].val;
      thismode=false;
      if (i==modeindx){
        n=i;
        m=j;
        o=k;
        thismode=true;
        //thismode=false;
      }
      else if (j==modeindx){
        n=j;
        m=i;
        o=k;
        thismode=false;
        //thismode=true;
      }
      else if (k==modeindx){
        n=k; 
        m=i;
        o=j;
        thismode=false;
        //thismode=true;
      }

      if (thismode){
        //fv += -(1.0/2.0)*mcc*( xm[m]*xm[l]*vm[n] - xm[n]*xm[l]*vm[m])*6.242e+18*1e-12; // eV/ps
        fv += -(1.0/2.0)*mcc*( xm[m]*xm[o]*vm[n])*6.242e+18*1e-12; // eV/ps THIS WORKS!
        //fv += (1.0/2.0)*mcc*( xm[n]*xm[m]*vm[o])*6.242e+18*1e-12; // eV/ps THIS WORKS!

      }
      
      /*
      if (i==j && j==k){ // iii
          f[i] -= input->mcc3[w].val*q[i]*q[i];
          f3[i] -= input->mcc3[w].val*q[i]*q[i];
      }
      
      if (i==j && i!=k){ //iij
          f[i] -= input->mcc3[w].val*q[j]*q[k];
          f[k] -= input->mcc3[w].val*q[i]*q[j];
          f3[i] -= input->mcc3[w].val*q[j]*q[k];
          f3[k] -= input->mcc3[w].val*q[i]*q[j];
      }
      
      if (i!=j && j==k){ // ijj
          f[i] -= input->mcc3[w].val*q[j]*q[k];
          f[j] -= input->mcc3[w].val*q[i]*q[k];
          f3[i] -= input->mcc3[w].val*q[j]*q[k];
          f3[j] -= input->mcc3[w].val*q[i]*q[k];
      }    
      
      if (i!=j && j!=k){ // ijk
          f[i] -= input->mcc3[w].val*q[j]*q[k];
          f[j] -= input->mcc3[w].val*q[i]*q[k];
          f[k] -= input->mcc3[w].val*q[i]*q[j];
          f3[i] -= input->mcc3[w].val*q[j]*q[k];
          f3[j] -= input->mcc3[w].val*q[i]*q[k];
          f3[k] -= input->mcc3[w].val*q[i]*q[j];
      }
      */
      //f[i] -= input->mcc3[f].val*q[
  }

  for (int w=0; w<input->nmcc4; w++){
      i=input->mcc4[w].i;
      j=input->mcc4[w].j;
      k=input->mcc4[w].k;
      l=input->mcc4[w].l;
      mcc = input->mcc4[w].val;
      thismode=false;
      if (i==modeindx){
        n=i;
        m=j;
        o=k;
        p=l;
        thismode=true;
        //thismode=false;
      }
      else if (j==modeindx){
        n=j;
        m=i;
        o=k;
        p=l;
        thismode=false;
        //thismode=true;
      }
      else if (k==modeindx){
        n=k; 
        m=i;
        o=j;
        p=l;
        thismode=false;
        //thismode=true;
      }
      else if (l==modeindx){
        n=l; 
        m=i;
        o=j;
        p=k;
        thismode=false;
        //thismode=true;
      }

      if (thismode){
        //fv += -(1.0/2.0)*mcc*( xm[m]*xm[l]*vm[n] - xm[n]*xm[l]*vm[m])*6.242e+18*1e-12; // eV/ps
        fv += -(1.0/6.0)*mcc*( xm[m]*xm[o]*xm[p]*vm[n])*6.242e+18*1e-12; // eV/ps THIS WORKS!
        //fv += -(1.0/6.0)*mcc*( xm[n]*xm[m]*xm[o]*vm[p])*6.242e+18*1e-12; // eV/ps THIS WORKS!
      }
      
  }

  // Print quantities
  fprintf(fh_em, "%f %.10e\n", em3->t*0.5*0.001, em[10]);
  fprintf(fh_fv, "%f %.10e\n", em3->t*0.5*0.001, fv);
}
