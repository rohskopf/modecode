# Simple observation of mode-mode energy transfer.

This example is the simplest situation of energy transfer between modes; one mode starts with finite
energy and transfers energy to two other modes. The system is a Si-Ge superlattice, and the initial
excited mode is an interface mode. 

First, run the LAMMPS input script to perform the MD simulation:

    lmp_mpi < in.run

By using `compute_mode` in LAMMPS, this input script will output a file `em.dat` which contains mode
energies at each timestep according to the order in the `INDICES` file. To observe the mode energies
simply run gnuplot:

    gnuplot < gnuplot_em_time

which will plot the mode energies as a function of time during the simulation. The simulation also 
calculated the power transfers at each timestep, and stored them in a file called `fv.dat`. This
file must be time-integrated to get the energy transfers. To integrate the power transfer
components* (e.g. iij, iik, ijj, ijk, etc.), do

    python integrate3.py

Calculate the sum of integrated powers with

    python integrate5.py

which creates integrated_power.dat. The integrated power may be viewed in gnuplot by doing

    gnuplot < gnuplot_em_compare
