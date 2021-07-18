/*
 em3.h

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>
#include <vector>
#include "mpi.h"


namespace EM3_NS
{
    class EM3
    {
    public:

        class Memory *memory;
        class Timer *timer;
        class Input *input;
        class Potential *potential;
        class Update *update;
        class Output *output;
        class Compute *compute;
        class Mode *mode;
        EM3(int, char **);
        ~EM3();

        FILE * fh_debug; // Debug file handle        

        void create(); // Function to create class instances
        void initialize(); // Initializes the program, reads input
        void finalize(); // Deallocates class instances
        int procs; // number of processes
        int rank; // rank of process
        int t; // timestep

    };
}

