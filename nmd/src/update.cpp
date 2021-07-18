/*
 update.cpp

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

#include "update.h"
#include "potential.h"
#include "memory.h"
#include "input.h"
#include "compute.h"
#include "mode.h"

using namespace std;

using namespace EM3_NS;

Update::Update(EM3 *em3) : Pointers(em3) {


    fh_debug = fopen("DEBUG_update", "w");
}

Update::~Update() 
{

    memory->deallocate(aa);

};

void Update::initialize()
{
    memory->allocate(aa, input->natoms,3);

}

void Update::integrate()
{

    double dt = input->dt;

    /* Real space integration */

    if (input->space == 0 || input->space == 2){
        double **fa = mode->fa;
        double **xa = mode->xa;
        double **va = mode->va;
        double *mass = input->mass;
        int *type = input->type;
        //int nmodes = mode->nmodes;

        for (int i=0; i<input->natoms; i++){
            for (int b=0; b<3; b++){
                aa[i][b] = fa[i][b]/mass[type[i]-1];
                //fprintf(fh_debug,"type[i]-1: %d\n", type[i]-1);
                //fprintf(fh_debug,"mass[type[i]-1]: %e\n", mass[type[i]-1]);
                //fprintf(fh_debug,"aa[%d][%d]: %e\n", i,b,aa[i][b]);
            }
        }

        for (int i=0; i<input->natoms; i++){

            va[i][0] += 0.5*dt*aa[i][0];
            va[i][1] += 0.5*dt*aa[i][1];
            va[i][2] += 0.5*dt*aa[i][2];

            xa[i][0] += dt*va[i][0];
            xa[i][1] += dt*va[i][1];
            xa[i][2] += dt*va[i][2];

        }
        // Calculate forces with these new positions

        potential->calculate();

        fa = mode->fa;

        for (int i=0; i<input->natoms; i++){
            for (int b=0; b<3; b++){
                aa[i][b] = fa[i][b]/mass[type[i]-1];
            }
        }

        for (int i=0; i<input->natoms; i++){
            va[i][0] += + 0.5*dt*aa[i][0];
            va[i][1] += + 0.5*dt*aa[i][1];
            va[i][2] += + 0.5*dt*aa[i][2];
        }

    } // real space integration

    /* Mode space integration */

    if (input->space == 1 || input->space == 2){ 

        // In mode space, the mode force is the acceleration
        double *am = mode->fm;
        double *xm = mode->xm;
        double *vm = mode->vm;
        int nmodes = mode->nmodes;

        /*
        for (int i=3; i<nmodes; i++){
            printf(" v[%d]: %e\n", i,v[i]);
        }
        */

        for (int i=3; i<nmodes; i++){
            vm[i] += 0.5*dt*am[i];
            xm[i] += dt*vm[i];
            //printf(" vm[%d]: %e\n", i,vm[i]);
            //printf(" xm[%d]: %e\n", i,xm[i]);
        }
        
        // Calculate forces with these new positions

        potential->calculate();

        am = mode->fm;

        // Perform final step of Verlet algorithm
        for (int i=3; i<nmodes; i++){
            vm[i] += 0.5*dt*am[i];
        }



    } // mode space integration


}
