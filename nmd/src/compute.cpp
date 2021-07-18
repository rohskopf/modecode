/*
 compute.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "mpi.h"
#include <math.h>       /* sqrt */
#include <random>

#include "compute.h"
#include "output.h"
#include "potential.h"
#include "input.h"
#include "update.h"
#include "memory.h"
#include "mode.h"

using namespace std;

using namespace EM3_NS;

Compute::Compute(EM3 *em3) : Pointers(em3) {

    //fh_debug = fopen("DEBUG_compute", "w");

}

Compute::~Compute() 
{
    //fclose(fh_debug);

};

void Compute::compute_ke()
{

    
    int natoms = input->natoms;
    //double **v = update->v;

    // Convert to atom coordinates
    
    /*

    ke = 0.0;
    for (int i=3; i<mode->nmodes; i++){
        ke += freqs[i]*freqs[i]*v[i]*v[i];
    }

    // Kinetic energy is half times the sum of squared speeds in LJ units

    ke=0.0;
    for (int i=0; i<natoms; i++){

        ke += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])/2.0;

    }

    etot = ke + potential->pe;

    temp = ke*(2.0/3.0)/natoms;
    */

    
    
}

void Compute::compute_pressure()
{

    /*
    int natoms = input->natoms;

    double p_xx = potential->p_xx;
    double p_yy = potential->p_yy;
    double p_zz = potential->p_zz;

    pressure = natoms*temp + (1.0/3.0)*(p_xx + p_yy + p_zz);
    pressure = pressure/(input->volume);
    */

}

void Compute::compute_msd()
{

    // Simply compute MSD for current positions against original positions

    /*
    double **x0 = input->x0;
    double **x = neighbor->x;
    int natoms = input->natoms;

    double dxi,dyi,dzi;
    double dist2; // distance squared
    double sum=0.0;
    for (int i=0; i<natoms; i++){

        // Compute distance moved for atom i

        dxi = x[i][0] - x0[i][0];
        dyi = x[i][1] - x0[i][1];
        dzi = x[i][2] - x0[i][2];

        dist2 = dxi*dxi + dyi*dyi + dzi*dzi;

        sum += dist2;
    }

    msd = sum/input->natoms;
    */

}

void Compute::compute_com()
{

    

}

void Compute::compute_rmsd()
{

    
    

}
