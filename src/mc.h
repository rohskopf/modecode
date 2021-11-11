
/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>
#include <vector>
#include "mpi.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"


namespace MC_NS
{
    class MC
    {
    public:


        LAMMPS_NS::LAMMPS *lmp;
        class Mem *mem;
        class In *in;
        class Verify *verify;
        class Ifc2mcc *ifc2mcc;
        class Asr *asr;
        class Compute *compute;
        class Tic *tic;
        class Qhgk *qhgk;
        class Visualize *visualize;
        class Postproc *postproc;
        //class PopTimer *poptimer;
        //class Config;
        MC(int, char **);
        ~MC();
        void create();
        void initialize();
        void run(int,char **);
        void finalize();
        int nprocs;
        int rank;

        int seed;

        double delta;
        double cutoff;
        int order;

        std::string task; // task to perform

    };
}

