/*
 pointers.h

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "mpi.h"

#include "em3.h"

namespace EM3_NS
{
    class Pointers
    {
    public:
        Pointers(EM3 *ptr) :
            em3(ptr),
            memory(ptr->memory),
            timer(ptr->timer),
            input(ptr->input),
            potential(ptr->potential),
            update(ptr->update),
            output(ptr->output),
            compute(ptr->compute),
            mode(ptr->mode)
            {}

        virtual ~Pointers() {}

    protected:
        EM3 *em3;
        Memory *&memory;
        Timer *&timer;
        Input *&input;
        Potential *&potential;
        Update *&update;
        Output *&output;
        Compute *&compute;
        Mode *&mode;
    };
}

