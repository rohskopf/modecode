/*
 em3.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi.h"

#include "memory.h"
#include "timer.h"
#include "input.h"
#include "potential.h"
#include "update.h"
#include "output.h"
#include "compute.h"
#include "mode.h"

using namespace std;

using namespace EM3_NS;

EM3::EM3(int narg, char **arg)
{

    //fh_debug = fopen("DEBUG_em3", "w");

    /************************** Set up MPI settings **************************/

    int color,key,global,local;
    MPI_Comm comm;

    // Split the communicators so that multiple instances can be run
    MPI_Comm_rank(MPI_COMM_WORLD, &global);
    color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
    key = global; 
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    MPI_Comm_rank(comm,&local);

    //  Get the number of processes.
    procs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
    //  Get the individual process ID.
    rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

    /************************** Initial Screen Output **************************/
    if (rank == 0)
    {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                            EM3 0.0                              +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " Running on " << procs << " procs" << std::endl;
    }
  
    timer = new Timer(this);

    if (rank == 0) std::cout << " Job started at " << timer->DateAndTime() << std::endl;

    /************************** Proceed with Program **************************/

    // Dynamically allocated pointers

    create();

    // Initialize system and mode quantities

    initialize();

    int natoms = input->natoms;
    // Calculate potential and forces

    potential->calculate();
    if (input->space == 1 || input->space==2) mode->mode2Atom();

    // Perform MD

    printf(" Running for %d steps.\n", input->nsteps);
    //printf("Step PE KE E T P MSD AverageNeigh\n");
    if (input->space!=2) printf("Step PE\n");
    if (input->space==2) printf("Step AtomPE ModePE\n");
    //compute->compute_ke();
    //compute->compute_pressure();
    //compute->compute_msd();
    //compute->compute_com();

    //printf("%d %f %f %f %f %f %f %f\n", 0, potential->pe/natoms, compute->ke/natoms, compute->etot/natoms, compute->temp,\
            compute->pressure, compute->msd, neighbor->neighavg);
    //output->write_xyz();
    //output->write_forces();

    t=0;
    //printf(" 0 %e\n", potential->pe*6.242e+18);
    if (input->space==0) printf(" %d %f \n", t, potential->pea*6.242e+18);
    if (input->space==1) printf(" %d %f \n", t, potential->pem*6.242e+18);
    if (input->space==2) printf(" %d %f %f\n", t, potential->pem*6.242e+18, potential->pea*6.242e+18);
    output->write_xyz();
    output->write_forces();
    for (t=1; t < input->nsteps+1;t++){

        //mode->atom2Mode(); // Convert atom coordinates to mode coordinates
        //mode->mode2Atom(); // Convert mode forces to atomic forces

        // Compute KE for temperature rescaling

        compute->compute_ke();
        mode->calcModeEnergy();

        // Update the timestep

        update->integrate();

        // Write output files

        if (t % input->nout == 0){

            //compute->compute_pressure();
            //compute->compute_msd();
            //compute->compute_com();

            //printf("%d %f %f %f %f %f %f %f\n", t, potential->pe/natoms, compute->ke/natoms, compute->etot/natoms, compute->temp,\
                    compute->pressure, compute->msd, neighbor->neighavg);

            if (input->space==0) printf(" %d %f \n", t, potential->pea*6.242e+18);
            if (input->space==1) printf(" %d %f \n", t, potential->pem*6.242e+18);
            if (input->space==2) printf(" %d %f %f\n", t, potential->pem*6.242e+18, potential->pea*6.242e+18);
            if (input->space == 1 || input->space==2) mode->mode2Atom();
            output->write_xyz();
            output->write_forces();
        }
    }
    

    // Delete dynamically allocated pointers

    finalize();

    if (rank == 0) std::cout << std::endl << " Job finished at " 
        << timer->DateAndTime() << std::endl;
    if (rank == 0) timer->print_elapsed();
    //printf("ASDFASDF\n");
}

void EM3::create()
{
    input = new Input(this);
    memory = new Memory(this);
    potential = new Potential(this);
    update = new Update(this);
    compute = new Compute(this);
    mode = new Mode(this);
    output = new Output(this);

}

void EM3::initialize()
{
    input->readconfig();
    input->readinput();
    //input->readconfig();
    input->readParams();
    input->initialize();
    mode->initialize();
    update->initialize();
}

EM3::~EM3()
{
    delete timer;
    //printf("BEfore input\n");
    delete input;
    delete mode;
    //printf("AFter input\n");
    delete output;

    //fclose(fh_debug);

}

void EM3::finalize()
{

    delete memory;
    delete potential;
    delete update;
    delete compute;


}

