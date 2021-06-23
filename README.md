# ModeCode 0.0

ModeCode is a massively parallel and modular program that aids in the study of vibrational modes.

![Alt Text](https://github.com/rohskopf/modecode/extended_mode.gif)

## Installation

Before we install ModeCode, we must first install LAMMPS as a shared library.

### Building LAMMPS as a shared library.

Go into the src/ directory of your LAMMPS installation and do:

    make mode=shlib g++_openmpi

or any other setting (aside from g++_openmpi) according to the instructions on the LAMMPS website 
(http://gensoft.pasteur.fr/docs/lammps/12Dec2018/Python_shlib.html).

Now we are ready to install ModeCode, by linking to the LAMMPS shared library.

### Installing ModeCode.

First download ModeCode from github:

    git clone https://github.com/rohskopf/modecode

Then edit src/Makefile to have the appropriate paths pointing towards your LAMMPS shared library
installation.

Go into src/ and install with:

    make

That's it! We just made a ModeCode executable called `mc`.

## Using ModeCode

The general format for running ModeCode is to do:

    mpirun -np P mc task setting1 setting2 setting3 ...

where 

- `P` is the number of processes.
- `task` is calculation/task we are performing, which typically
referes to the C++ module/class being used. 
- `settings` further specify which function in the 
`task` is used, or input values needed to run the function.

There may be few or many `settings`, depending on which `task` we use, so let's break it down.

***

### Finite difference (`fd`) task.

This task uses finite difference to extract the 2nd, 3rd, or 4th order interatomic force constants
(IFCs). 

The general use of this task is:

    mpirun -np P mc fd delta cutoff order

where 

- `P` is the number of processes to split the IFC calculations over.
- `fd` refers to the finite difference task.
- `delta` is the finite difference step size, anywhere from 0.0001 to 0.01 depending on the 
  stiffness of your material. 
- `cutoff` is the interatomic interaction cutoff for force constants, in whatever units declared by
  LAMMPS.
- `order` is the order of finite difference, e.g. 2, 3, 4.

#### Outputs.

- FC2, FC3, or FC4 depending on the `order` parameter.

***

### Acoustic sum rule (`asr`) task.

This task calculates the self-interaction 2nd order IFCs, given the file FC2 which does not contain
self terms. 

The general use of this task is:

    mc asr order

where 

- `asr` refers to the ASR task.
- `order` refers to the IFC order, although for now ModeCode only supports 2nd order ASR.

#### Outputs.

- FC2_ASR
  - Same as FC2, except includes self-terms.

***

### Diagonalization (`diag`) task.

This task diagonalizes the dynamical matrix to get the mode frequencies and eigenvectors.

Currently this is not implemented in the C++ code; we use a simple Python script instead. 

See tools/calc_eig3.py.
This takes FC2_ASR, DATA (LAMMPS data file), and then you must edit the molecular weights to 
get the desired dynamical matrix. Simply do:

    python calc_eig3.py

#### Outputs.

- FREQUENCIES
- EMAT

***

### Compute (`compute`) task.

This task computes various quantities associated with the modes. 

The general use of this task is:

    mc compute compute_task setting1 setting2 ...

where
- `compute` refers to the compute task.
- `compute_task` is the sub-task which refers to the quantity we are computing.
- `settings` determine the necessary inputs for the calculation.

#### Compute tasks (`compute_task`).

There are a few separate sub-tasks which will be explained here.

##### Participation ratio (`pr`).

Run with:

    mc compute pr natoms

where
- `pr` refers to the participation ratio sub-task.
- `natoms` is the number of atoms in the system.

##### Eigenvector spatial parameter (`esp`).

Run with:

    mc compute esp

where
- `esp` refers to the ESP task.


***

### IFC to MCC conversion (`ifc2mcc`) task.

This task converts IFCs to mode coupling constants (MCCs), in various different ways.

Documentation coming soon.

#### Outputs.

- Coming soon.

***

## Tutorial: Crystalline silicon.

This is nice to make sure everything works.

Go into examples/Si_8atoms.

### 2nd order finite difference.

    mpirun -np 2 mc fd 0.01 2

This creates the FC2 file.

### Acoustic sum rule correction.

    mc asr 2

This creates the FC2_ASR file, same as FC2 but with self-terms. 

### Get the modes.

There is a custom calc_eig3.py script in this directory. The inputs and molecular weights are
set so that the correct dynamical matrix is made. Run with:

    python calc_eig3.py
