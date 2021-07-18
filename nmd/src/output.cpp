/*
 output.cpp

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

#include "output.h"
#include "potential.h"
#include "input.h"
#include "memory.h"
#include "mode.h"

using namespace std;

using namespace EM3_NS;

Output::Output(EM3 *em3) : Pointers(em3) {

    fh_xyz_a = fopen("dump_atom.xyz", "w");
    fh_xyz_m = fopen("dump_mode.xyz", "w");
    fh_forces = fopen("dump.forces", "w");
    fh_forces2 = fopen("dump.forces2", "w");
    fh_forces3 = fopen("dump.forces3", "w");
    fh_forces4 = fopen("dump.forces4", "w");
    fh_forces_mta = fopen("dump.forces_mta", "w");
    fh_forces_mta2 = fopen("dump.forces_mta2", "w");
    fh_forces_mta3 = fopen("dump.forces_mta3", "w");
    fh_forces_mta4 = fopen("dump.forces_mta4", "w");
    fh_dfile = fopen("DFILE", "w");
    fh_ffile = fopen("FFILE", "w");


}

Output::~Output() 
{
    fclose(fh_xyz_a);
    fclose(fh_xyz_m);
    fclose(fh_forces);
    fclose(fh_forces2);
    fclose(fh_forces3);
    fclose(fh_forces4);
    fclose(fh_forces_mta);
    fclose(fh_forces_mta2);
    fclose(fh_forces_mta3);
    fclose(fh_forces_mta4);
    fclose(fh_dfile);
    fclose(fh_ffile);

};

void Output::write_xyz()
{

    /* Real space xyz */
    if (input->space == 0 || input->space==2){
        //printf("ASDF\n");
        /* XYZ file */
        fprintf(fh_xyz_a, "%d\n", input->natoms);
        fprintf(fh_xyz_a, "Atoms. Timestep: %d\n", em3->t);
        for (int i=0; i<input->natoms; i++){
            fprintf(fh_xyz_a, "%d %e %e %e\n", input->type[i],mode->xa[i][0]*1e10, mode->xa[i][1]*1e10, mode->xa[i][2]*1e10);
        }

        /* Alamode DFILE */

        if (em3->t != 0){
            double ux,uy,uz;
            for (int i=0; i<input->natoms; i++){
                ux = mode->xa[i][0] - input->xa0[i][0];
                uy = mode->xa[i][1] - input->xa0[i][1];
                uz = mode->xa[i][2] - input->xa0[i][2];
                // Convert from m to Bohr
                ux = ux*1.89e+10;
                uy = uy*1.89e+10;
                uz = uz*1.89e+10;
                fprintf(fh_dfile, "%f %f %f\n", ux,uy,uz);
            }
        }

    }

    /* Mode space xyz */
    if (input->space == 1 || input->space==2){
        fprintf(fh_xyz_m, "%d\n", input->natoms);
        fprintf(fh_xyz_m, "Atoms. Timestep: %d\n", em3->t);
        for (int i=0; i<input->natoms; i++){
            fprintf(fh_xyz_m, "%d %e %e %e\n", input->type[i],mode->xmta[i][0]*1e10, mode->xmta[i][1]*1e10, mode->xmta[i][2]*1e10);
        }
    }
    


}

void Output::write_forces()
{
    /*
    fprintf(fh_forces, "Modes. Timestep: %d\n", em3->t);
    for (int i=0; i<mode->nmodes; i++){
        fprintf(fh_forces, "%f\n", mode->f[i]);
    }
    */

    

    /* Real space */
    if (input->space == 0 || input->space==2){

        /* Atom forces */
        fprintf(fh_forces, "Timestep: %d\n", em3->t);
        for (int i=0; i<input->natoms; i++){
            fprintf(fh_forces, "%e %e %e\n", mode->fa[i][0],mode->fa[i][1],mode->fa[i][2]);
        }
        if (input->order>=2){
            /* Atom forces (2nd order) */
            fprintf(fh_forces2, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces2, "%e %e %e\n", mode->fa2[i][0],mode->fa2[i][1],mode->fa2[i][2]);
            }
        }
        if (input->order>=3){
            /* Atom forces (3rd order) */
            fprintf(fh_forces3, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces3, "%e %e %e\n", mode->fa3[i][0],mode->fa3[i][1],mode->fa3[i][2]);
            }
        }
        if (input->order>=4){
            /* Atom forces (4th order) */
            fprintf(fh_forces4, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces4, "%e %e %e\n", mode->fa4[i][0],mode->fa4[i][1],mode->fa4[i][2]);
            }
        }
        /* Alamode FFILE */
        if (em3->t != 0){
            double fx,fy,fz;
            for (int i=0; i<input->natoms; i++){
                fx = mode->fa[i][0];
                fy = mode->fa[i][1];
                fz = mode->fa[i][2];
                // Convert from N/m to Ryd/Bohr
                fx = fx*4.587420897e+17*(1.0/1.89e+10);
                fy = fy*4.587420897e+17*(1.0/1.89e+10);
                fz = fz*4.587420897e+17*(1.0/1.89e+10);
                fprintf(fh_ffile, "%f %f %f\n", fx,fy,fz);
            }
        }

    }

    /* Mode space forces converted to atomic forces */
    if (input->space == 1 || input->space==2){
    
        /* Atom forces converted from mode forces */
        fprintf(fh_forces_mta, "Timestep: %d\n", em3->t);
        for (int i=0; i<input->natoms; i++){
            fprintf(fh_forces_mta, "%e %e %e\n", mode->fmta[i][0],mode->fmta[i][1],mode->fmta[i][2]);
        }
        if (input->order>=2){
            /* Atom forces converted from mode forces (2nd order) */
            fprintf(fh_forces_mta2, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces_mta2, "%e %e %e\n", mode->fmta2[i][0],mode->fmta2[i][1],mode->fmta2[i][2]);
            }
        }
        if (input->order>=3){
            /* Atom forces converted from mode forces (3rd order)*/
            fprintf(fh_forces_mta3, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces_mta3, "%e %e %e\n", mode->fmta3[i][0],mode->fmta3[i][1],mode->fmta3[i][2]);
            }
        }
        if (input->order>=4){
            /* Atom forces converted from mode forces (4th order)*/
            fprintf(fh_forces_mta4, "Timestep: %d\n", em3->t);
            for (int i=0; i<input->natoms; i++){
                fprintf(fh_forces_mta4, "%e %e %e\n", mode->fmta4[i][0],mode->fmta4[i][1],mode->fmta4[i][2]);
            }
        }
    } 
    

}
