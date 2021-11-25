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
#include "visualize.h"
#include "postproc.h"
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
        std::cout << " +                          ModeCode 0.0                           +" << std::endl;
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
    if (task=="visualize"){
        visualize = new Visualize(this);
    }
    if (task=="postproc"){
        postproc = new Postproc(this);
    }
    //else{
    //    printf(" INVALID TASK.\n");
    //}
    

}

void MC::initialize()
{
    if (rank==0) printf(" Reading INPUT.\n");
    in->readInput();
    if (rank==0) printf(" Done reading INPUT.\n");

}

void MC::run(int nargs, char **args)
{

  if (task=="fd"){
    if (nargs!=6){
        printf(" Need 6 arguments for FD! fd delta cutoff tolerance order\n");
        exit(1);
    }
    delta = atof(args[2]);
    cutoff = atof(args[3]);
    order = atoi(args[5]);
    if (order == 2) in->tol2 = atof(args[4]);
    if (order == 3) in->tol3 = atof(args[4]);
    if (order == 4) in->tol4 = atof(args[4]);
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
          ifc2mcc->go(atof(args[4]));  
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
          //order = atoi(args[3]);
          ifc2mcc->go_gv(atoi(args[3]));  
      }

  }

  if (task=="asr"){
      if (nargs != 4){
        printf("Need 3 args for ASR! modecode asr order tolerance\n");
        exit(1);
      }
      order = atoi(args[2]);
      double tolerance = atof(args[3]);
      asr->go(order, tolerance);
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

  if (task=="visualize"){

      if (nargs != 8){
        printf("Need 8 args for visualization! visualize temperature n_indx timestep scale_factor largest_gv_setting additional_time\n");
        printf("You gave %d args.\n", nargs);
        exit(1);
      }
      visualize->temperature = atof(args[2]);
      visualize->n_indx = atoi(args[3]);
      visualize->timestep = atof(args[4]);
      visualize->scale_factor = atof(args[5]);
      visualize->largest_gv_setting = atoi(args[6]); // 0 - take largest GV within +/- 10 modes, 1 - take largest GV across all modes.
      visualize->additional_time = atof(args[7]);

      visualize->initialize();
      visualize->readEmat();
      visualize->readGV();

      visualize->calcInitialState();
      visualize->calcTimeDependence();


  }
  
  if (task=="postproc"){
  
    if (rank==0) printf(" Post processing simulation data.\n");
    //if (rank==0) printf(" Reading EMAT.\n");
    
    //postproc->readEmat();
    //if (rank==0) printf(" Finished reading EMAT.\n");
    
    int task_postproc = atoi(args[2]); // 1 - loop over ensembles, calculate FFT and |FFT|^2 (power spectrum) of all mode amplitudes and velocities, calculate DOS overlap for all pairs, then ensemble average.
    postproc->task=task_postproc;
    
    postproc->initialize();
    
    postproc->ensemble_dirname = std::string(args[3]); // name of ensemble directories, not including the ensemble number.
                                                        // E.g. "e_100_" is the general name, but "e_100_1" is the dirname of a particular ensemble directory.
    //int task_postproc = atoi(args[2]);
    //printf(" Postproc task: %d\n", task_postproc);
    
    if (task_postproc==1){
      postproc->ntimesteps = atoi(args[4]);
      postproc->sampling_interval = atof(args[5]);
      postproc->nens = atoi(args[6]);
      postproc->overlap_output_tag = atoi(args[7]);
      postproc->task1();
    }
    // Calculate auto-correlation of mode action for certain pairs.
    //postproc->calc1();
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
    if (task=="visualize"){
        delete visualize;
    }
    if (task=="postproc"){
        delete postproc;
    }

}

MC::~MC()
{
    delete lmp;

};
