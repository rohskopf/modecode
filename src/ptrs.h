

#pragma once

#include "mpi.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

#include "mc.h"

namespace MC_NS
{
    class Ptrs
    {
    public:
        Ptrs(MC *ptr) :
            mc(ptr),
            mem(ptr->mem),
            //poptimer(ptr->poptimer),
            //popinput(ptr->popinput),
            //config(ptr->config),
            lmp(ptr->lmp)
            {}

        virtual ~Ptrs() {}

    protected:
        MC *mc;
        Mem *&mem;
        //PopTimer *&poptimer;
        //PopInput *&popinput;
        //Config *&config;
        LAMMPS_NS::LAMMPS *&lmp;
    };
}

