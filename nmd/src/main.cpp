/*
 main.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <stdlib.h>
#include <iostream>
#include "em3.h"
#include "mpi.h"

using namespace EM3_NS;

int main(int argc, char **argv)
{

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* Begin a EM3 instance */
  EM3 *em3 = new EM3(argc, argv);

  /* Delete the memory */

  delete em3;

  /* Close MPI */
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

