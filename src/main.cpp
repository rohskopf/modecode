#include <stdlib.h>
#include <iostream>
#include "mc.h"
#include "mpi.h"

using namespace MC_NS;

int main(int argc, char **argv)
{

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* Begin a MC instance */
  MC *mc = new MC(argc, argv);

  /* Delete the memory */
  delete mc;

  /* Close MPI */
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

