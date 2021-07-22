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
    memory->deallocate(xa_old);
    memory->deallocate(xm_old);

};

void Update::initialize()
{
    memory->allocate(aa, input->natoms,3);
    memory->allocate(xa_old, input->natoms,3);
    memory->allocate(xm_old, input->natoms*3);

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

            //fprintf(fh_debug, "%e %e %e\n", xa[i][0], xa[i][1], xa[i][2]);

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


void Update::integrate2()
{

    double dt = input->dt;

    // Calculate forces

    potential->calculate();

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

        if (em3->t==1){
            for (int i=0; i<input->natoms; i++){

                xa_old[i][0] = xa[i][0];
                xa_old[i][1] = xa[i][1];
                xa_old[i][2] = xa[i][2];

                //fprintf(fh_debug, "%e %e %e\n", xa_old[i][0], xa_old[i][1], xa_old[i][2]);

                xa[i][0] += va[i][0]*dt + 0.5*aa[i][0]*dt*dt;
                xa[i][1] += va[i][1]*dt + 0.5*aa[i][1]*dt*dt;
                xa[i][2] += va[i][2]*dt + 0.5*aa[i][2]*dt*dt;

                //fprintf(fh_debug, "%e %e %e\n", xa[i][0], xa[i][1], xa[i][2]);

                //fprintf(fh_debug, "x0: %e %e %e\n", xa_old[i][0], xa_old[i][1], xa_old[i][2]);
                //fprintf(fh_debug, "x1: %e %e %e\n", xa[i][0], xa[i][1], xa[i][2]);
                //fprintf(fh_debug, "\n");

            }

            //fprintf(fh_debug, "END OF T=1\n");

        } // if em3->t == 0

        else{
            for (int i=0; i<input->natoms; i++){

                //xa_old[i][0] = xa[i][0];
                //xa_old[i][1] = xa[i][1];
                //xa_old[i][2] = xa[i][2];

                //fprintf(fh_debug, "%e %e %e\n", xa_old[i][0], xa_old[i][1], xa_old[i][2]);

                //fprintf(fh_debug, "x0: %e %e %e\n", xa_old[i][0], xa_old[i][1], xa_old[i][2]);
                //fprintf(fh_debug, "x1: %e %e %e\n", xa[i][0], xa[i][1], xa[i][2]);

                xa[i][0] = 2.0*xa[i][0] - xa_old[i][0] + aa[i][0]*dt*dt;
                xa[i][1] = 2.0*xa[i][1] - xa_old[i][1] + aa[i][1]*dt*dt;
                xa[i][2] = 2.0*xa[i][2] - xa_old[i][2] + aa[i][2]*dt*dt;

                // Get new x_old
                xa_old[i][0] = 0.5*(xa[i][0] + xa_old[i][0] - aa[i][0]*dt*dt);
                xa_old[i][1] = 0.5*(xa[i][1] + xa_old[i][1] - aa[i][1]*dt*dt);
                xa_old[i][2] = 0.5*(xa[i][2] + xa_old[i][2] - aa[i][2]*dt*dt);

                //fprintf(fh_debug, "x1: %e %e %e\n", xa_old[i][0], xa_old[i][1], xa_old[i][2]);
                //fprintf(fh_debug, "x2: %e %e %e\n", xa[i][0], xa[i][1], xa[i][2]);
                //fprintf(fh_debug, "\n");

                va[i][0] = (xa[i][0]-xa_old[i][0])/dt;
                va[i][1] = (xa[i][1]-xa_old[i][1])/dt;
                va[i][2] = (xa[i][2]-xa_old[i][2])/dt;


            }
        } // else

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

        if (em3->t==1){
            for (int i=3; i<nmodes; i++){

                xm_old[i] = xm[i];
                xm[i] += vm[i]*dt + 0.5*am[i]*dt*dt;

            }
        } // if em3->t == 0

        else{

            for (int i=3; i<nmodes; i++){

                xm[i] = 2.0*xm[i] - xm_old[i] + am[i]*dt*dt;

                // Get new x_old
                xm_old[i] = 0.5*(xm[i] + xm_old[i] - am[i]*dt*dt);

                vm[i] = (xm[i]-xm_old[i])/dt;

            }

        } // else

    } // mode space integration

    // Calculate forces

    potential->calculate();



}
