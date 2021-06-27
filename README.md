# ModeCode 0.0

Vibrational modes are collective movements in a system. In atomic systems, these collective motions are composed of atoms and they vibrate at specific frequencies. All the motion in an atomic system is determined by Newton's 2nd law, and this motion may be decomposed into the individual modes. This is illustrated below for atoms in a diamond crystal with seemingly random thermal vibrations at finite temperature.

![Alt Text](https://github.com/rohskopf/modecode/blob/main/mode_decomposition.gif)

These motions influence all phenomena that we observe in materials, including thermal and electrical transport, mass diffusion, and chemical reactions. For example consider the vibrational mode of atoms in a superlattice below, modes like these aid in the transfer of heat through superlattice thermoelectric devices, and understanding their behavior will help us engineer better thermoelectric devices.

![Alt Text](https://github.com/rohskopf/modecode/blob/main/extended_mode.gif)

ModeCode is a massively parallel and modular program that aids in the study of vibrational modes. This program is the first of its kind, as all other open-source programs only deal with vibrational modes/phonons in ideal crystals. By interfacing with a large library of interatomic potentials and alleviating the assumption of periodicity, ModeCode allows for the study of vibrational modes in all materials. 

# Installation

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

    make clean
    make

That's it! We just made a ModeCode executable called `modecode`.

# Using ModeCode

ModeCode is simple - you declare a calculation task and settings and must have an INPUT file in the directory that you run it in. These details are explained below.

The general format for running ModeCode is to do:

    mpirun -np P modecode task setting1 setting2 setting3 ...

where 

- `P` is the number of processes.
- `task` is calculation/task we are performing, which typically
referes to the C++ module/class being used. 
- `settings` further specify which function in the 
`task` is used, or input values needed to run the function.

There may be few or many `settings`, depending on which `task` we use. Before we get into the multiple tasks, it's first important to understand the input files required to use ModeCode.

### INPUT file.

This file is composed of LAMMPS commands that declare your system geometry, potential, neighborlist settings, etc. Please refer to the LAMMPS documentation to declare desired settings for your system. It's important to declare LAMMPS settings here, because some tasks in ModeCode such as finite difference will evaluate the potential at many different geometries; we therefore need to declare neighborlist settings, a system geometry via LAMMPS data files, atom styles, and pair styles. 

Other than the INPUT file, other files required depend on whatever your INPUT file uses. For example if you use the LAMMPS `read_data` command in your INPUT file, you need to also include the data file in your directory. If your LAMMPS pair style has a file it reads parameters from, then that file must also be included in the directory. Any file used by your INPUT file must also be included in the directory. 

Now that we understand the different inputs used by ModeCode, let's consider the different tasks or calculations that are possible.

### Finite difference (`fd`) task.

This task uses finite difference to extract the 2nd, 3rd, or 4th order interatomic force constants
(IFCs). 

The general use of this task is:

    mpirun -np P modecode fd delta cutoff tolerance order

where 

- `P` is the number of processes to split the IFC calculations over.
- `fd` refers to the finite difference task.
- `delta` is the finite difference step size, a number depending on the 
  stiffness of your material with the same units declared by LAMMPS in the INPUT file. 
- `cutoff` is the interatomic interaction cutoff for force constants, in whatever units declared by
  LAMMPS.
- `tolerance` tells the program to ignore force constants below this absolute value. Units are determined by LAMMPS.
- `order` is the order of finite difference, e.g. 2, 3, 4.

#### Outputs.

- FC2, FC3, or FC4 depending on the `order` parameter.

***

### Acoustic sum rule (`asr`) task.

This task calculates the self-interaction 2nd order IFCs, given the file FC2 which does not contain
self terms. 

The general use of this task is:

    modecode asr order

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

    modecode compute compute_task setting1 setting2 ...

where
- `compute` refers to the compute task.
- `compute_task` is the sub-task which refers to the quantity we are computing.
- `settings` determine the necessary inputs for the calculation.

#### Compute tasks (`compute_task`).

There are a few separate sub-tasks which will be explained here.

##### Participation ratio (`pr`).

Run with:

    modecode compute pr natoms

where
- `pr` refers to the participation ratio sub-task.
- `natoms` is the number of atoms in the system.

##### Eigenvector spatial parameter (`esp`).

Run with:

    modecode compute esp

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

    mpirun -np 2 modecode fd 0.01 2

This creates the FC2 file.

### Acoustic sum rule correction.

    modecode asr 2

This creates the FC2_ASR file, same as FC2 but with self-terms. 

### Get the modes.

There is a custom calc_eig3.py script in this directory. The inputs and molecular weights are
set so that the correct dynamical matrix is made. Run with:

    python calc_eig3.py
