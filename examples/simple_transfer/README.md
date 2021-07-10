# Simple observation of mode-mode energy transfer.

This example is the simplest situation of energy transfer between modes; one mode starts with finite
energy and transfers energy to two other modes. The system is a Si-Ge superlattice, and the initial
excited mode is an interface mode. 

First you must calculate the modes for this Si/Ge superlattice to obtain the EMAT file, which is not
included here due to the large file size. To do this,
follow the same procedure in `modecode/examples/Si_8atoms` to get the eigenvectors (EMAT) and 
frequencies.

ModeCode may
also be used to get the MCC3s, but they are also supplied here for the modes of interest. The EMAT
file should be placed in this directory, and the LAMMPS simulations are run in the `files/` 
directory, because the LAMMPS `compute mode` and `fix mode` read `../EMAT`. 

Go into the `files/` directory to begin.

Run the LAMMPS input script to perform the MD simulation:

    lmp_mpi < in.run

By using `compute_mode` in LAMMPS, this input script will output a file `em.dat` which contains mode
energies at each timestep according to the order in the `INDICES` file. To observe the mode energies
simply run gnuplot:

    gnuplot < gnuplot_em_time

which will plot the mode energies as a function of time during the simulation. The simulation also 
calculated the power transfers at each timestep, and stored them in a file called `fv.dat`. This
file must be time-integrated to get the energy transfers, and this is achieved by running

    python integrate2.py

which creates integrated_power.dat. Or to integrate the power transfer components, do

    python integrate5.py

which creates integrated_fv_sum.dat. This integrated power may be viewed in gnuplot by doing

    gnuplot < gnuplot_em_compare
