#include <iostream>
#include <iomanip>
#include <vector>

/* ModeCode (MC) include files */
#include "mem.h"
//#include "poptimer.h"
#include "in.h"
#include "verify.h"
#include "ifc2mcc.h"
#include "asr.h"
#include "compute.h"
#include "tic.h"
#include "qhgk.h"
//#include "config.h"

/* MPI include file */
#include "mpi.h"

using namespace std;

using namespace MC_NS;

MC::MC(int narg, char **arg)
{

    /************************** Set up MPI settings **************************/

    int color,key,global,local;
    MPI_Comm comm;

    // Split the communicators so that multiple instances of LAMMPS can be run
    MPI_Comm_rank(MPI_COMM_WORLD, &global);
    color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
    key = global; 
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    MPI_Comm_rank(comm,&local);

    //  Get the number of processes.
    nprocs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
    //  Get the individual process ID.
    rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

    /************************** Initial Screen Output **************************/
    if (rank == 0)
    {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                        IFC Calculator 0.0                       +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        //std::cout << " Running on " << nprocs << " procs" << std::endl;
    }
  
    //poptimer = new PopTimer(this);

    //if (rank == 0) std::cout << " Job started at " << poptimer->DateAndTime() << std::endl;

    /************************** Create LAMMPS Pointer **************************/
    char *args1[] = {
        (char *) "lmp",
        (char *) "-screen",
        (char*) "none",
    0};

    lmp = new LAMMPS_NS::LAMMPS(3,args1,comm);

    /************************** Proceed **************************/

    // Extract desired task
    task = std::string(arg[1]);
    //std::cout << task << std::endl;

    // Dynamically allocated pointers
    create();

    // Grab user input to initialize settings
    initialize();
    
    // Run the code and desired task
    run(narg,arg);

    // Delete dynamically allocated pointers
    finalize();

    /*
    if (rank == 0) std::cout << std::endl << " Job finished at " 
        << poptimer->DateAndTime() << std::endl;
    */
}

void MC::create()
{
    mem = new Mem(this);
    in = new In(this);

    if (task=="fd"){
    }
    if (task=="verify"){
        verify = new Verify(this);
    }
    if (task=="ifc2mcc"){
        ifc2mcc = new Ifc2mcc(this);
    }
    if (task=="asr"){
        asr = new Asr(this);
    }
    if (task=="compute"){
        compute = new Compute(this);
    }
    if (task=="tic"){
        tic = new Tic(this);
    }
    if (task=="qhgk"){
        qhgk = new Qhgk(this);
    }
    //else{
    //    printf(" INVALID TASK.\n");
    //}
    

}

void MC::initialize()
{
    if (rank==0) printf(" Reading INPUT.\n");
    in->readInput();

}

void MC::run(int nargs, char **args)
{

  if (task=="fd"){
    if (nargs!=5){
        printf(" Need 5 arguments for FD! fd delta cutoff order\n");
        exit(1);
    }
    delta = atof(args[2]);
    cutoff = atof(args[3]);
    order = atoi(args[4]);
    if (rank==0){
      printf(" Performing %d order FD on %d procs.\n",order,nprocs);
      printf(" Cutoff: %f\n", cutoff);
    }
    in->calcFC2();
  }

  if (task=="verify"){
    string verification = std::string(args[2]);
    string mcc3_filename = std::string(args[6]);
    verify->readMCC(mcc3_filename.c_str());
    if (verification=="random"){
        printf(" Verifying an existing TEP with random atomic displacements.\n");
        int nconfigs = atoi(args[3]); // number of random configurations to compare against
        double mag = atof(args[4]); // absolute value of random mode amplitude magnitudes
        int seed = atof(args[5]); // absolute value of random mode amplitude magnitudes
        verify->compareRandom(nconfigs,mag,seed);
    }
    else printf(" INVALID VERIFICATION.\n");
  }

  if (task=="ifc2mcc"){
      ifc2mcc->task = atoi(args[2]);
      // Convert IFC 2 MCC
      if (ifc2mcc->task==0){
          order = atoi(args[3]);
          ifc2mcc->go();  
      }
      // Extract a range of MCC3s into a single MCC3 file.
      else if (ifc2mcc->task==1){
          ifc2mcc->lowmode = atoi(args[3]);
          ifc2mcc->highmode = atoi(args[4]);
          printf(" Extracting MCC3s for modes %d - %d.\n", ifc2mcc->lowmode, ifc2mcc->highmode);
          ifc2mcc->extract();
      }
      // For a list of modes "i", take the average of all MCC3s with modes j and k.
      else if (ifc2mcc->task==2){
          //printf(" Extracting MCC3s for a list of modes.\n");
          int startfile = atoi(args[3]);
          int endfile = atoi(args[4]);
          ifc2mcc->average(startfile, endfile);
      }

      // For a list of modes "i", extract all MCC3s associated with triplets of these modes.
      else if (ifc2mcc->task==3){
          //printf(" Extracting MCC3s for a list of modes.\n");
          int startfile = atoi(args[3]);
          int endfile = atoi(args[4]);
          ifc2mcc->extractFew(startfile,endfile);
      }

      // For a list of modes "i", calculate all MCC3s associated with triplets of these modes.
      else if (ifc2mcc->task==4){
          //printf(" Extracting MCC3s for a list of modes.\n");
          order = atoi(args[3]);
          ifc2mcc->calcFew();
      }

      // Convert IFC 2 MCC within spatial regions (based on atom types 1 and 2). 
      if (ifc2mcc->task==5){
          order = atoi(args[3]);
          ifc2mcc->go_spatial();  
      }

      // Convert IFC 2 MCC3 for a particular mode n1. 
      if (ifc2mcc->task==6){
          //order = atoi(args[3]);
          ifc2mcc->go_n1(atoi(args[3]));  
      }

      // Extract all SMCC2s into a single SMCC2 file.
      else if (ifc2mcc->task==7){
          printf(" Extracting SMCC2s.\n");
          ifc2mcc->extract_smcc2(atoi(args[3]));
      }

      // Calculate generalized velocities. 
      if (ifc2mcc->task==8){
          order = atoi(args[3]);
          ifc2mcc->go_gv();  
      }

  }

  if (task=="asr"){
      if (nargs != 3){
        printf("Need 3 args for ASR!\n");
        exit(1);
      }
      order = atoi(args[2]);
      asr->go(order);
  }

  if (task=="compute"){

      //printf(" Compute.\n");
      // Extract compute task.
      string task_compute = std::string(args[2]);
      printf(" Compute task: %s\n", task_compute.c_str());

      if (task_compute=="pr"){
          if (nargs != 4){
            printf("Need 4 args for PR calculation! compute pr natoms\n");
            exit(1);
          }
          compute->natoms=atoi(args[3]);
          compute->pr_call = true;
          printf(" Calculating participation ratio based on mode eigenvectors.\n");
          compute->readEmat();
          compute->participationRatio();
      }

      if (task_compute=="esp"){
          if (nargs != 3){
            printf("Need 3 args for ESP calculation! compute esp\n");
            exit(1);
          }
          //compute->natoms=atoi(args[3]);
          compute->esp_call = true;
          //printf(" Reading EMAT\n");
          //compute->readEmat();
          printf(" Calculating ESP.\n");
          compute->eigenvectorSpatialParameter();
      }

  }

  if (task=="tic"){
      // Compute thermal interface conductance.
      tic->go();

  }

  if (task=="qhgk"){
      // Compute thermal interface conductance.
      // First argument - 0 or 1 for classical or quantum, respectively.
      // Second argument - temperature in K.
      qhgk->go(atoi(args[2]),atof(args[3]));

  }

  
}

void MC::finalize()
{

    delete in;
    delete mem;

    if (task=="fd"){
    }
    if (task=="verify"){
        delete verify;
    }
    if (task=="ifc2mcc"){
        delete ifc2mcc;
    }
    if (task=="asr"){
        delete asr;
    }
    if (task=="compute"){
        delete compute;
    }
    if (task=="tic"){
        delete tic;
    }
    if (task=="qhgk"){
        delete qhgk;
    }

}

MC::~MC()
{
    delete lmp;

};
