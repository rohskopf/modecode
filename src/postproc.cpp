/*
This class post processes simulation data to analyze the modes.
Current functions include:

*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"

#include "postproc.h"
#include "mem.h"
#include "in.h"
//#include "config.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "pair.h"
#include "domain.h"

#include <ctime>

using namespace std;

using namespace MC_NS;

Postproc::Postproc(MC *mc) : Ptrs(mc) {
    //fh_debug = fopen("DEBUG_ASR","w");
    //fh_fc2 = fopen("FC2_ASR","w");
    
    rank = mc->rank;

}

Postproc::~Postproc() 
{

    //fclose(fh_debug);
    //fclose(fh_fc2);

    //if (order >= 3) mem->deallocate(fc3);
    //if (order >= 4) mem->deallocate(fc4);
    mem->deallocate(emat);

};

/*
Initialize
*/

void Postproc::initialize()
{

  if (rank==0) printf(" Initializing.\n");
  natoms = lmp->atom->natoms;

}

void Postproc::calc1()
{

  if (rank==0) printf(" Calculating autocorrelation of mode action pairs.\n");
  double **test;
  double **test2;
  
  long long int n1 = 400000;
  long long int n2 = 3000;
  
  mem->allocate(test,n1,n2);
  mem->allocate(test2,n1,n2);
  size_t size_of_one = sizeof(test[0]);
  printf(" %ld\n", size_of_one);

  long long int n = n1*n2*size_of_one;
  std::cout << n << std::endl;
  n = n/1e6;
  printf(" MB: %lld\n", n);
  
  mem->deallocate(test);
  mem->deallocate(test2);
}

void Postproc::readEmat()
{

    if (rank==0) printf(" Reading EMAT.\n");
    if (rank==0) printf(" Reading eigenvector matrix for %d atoms.\n",natoms);

    ifstream readfile;

    mem->allocate(emat,natoms*3,natoms*3);

    for (int i = 0; i < natoms*3; i++) {
        for (int j = 0; j < natoms*3; j++) {
            emat[i][j] = 0.0;
        }
    }

    readfile.open("EMAT");
    //readfile2.open("ev_real.txt");
  
    if (!readfile.is_open()) {
        cout<<"Unable to open the file!"<<endl;
        exit(1);
    }

    //printf("natoms: %d\n",  natoms);
    for (int i=0;i<3*natoms;i++){
        for (int j=0;j<3*natoms;j++){
            readfile>>emat[i][j];
        }
    }
}
