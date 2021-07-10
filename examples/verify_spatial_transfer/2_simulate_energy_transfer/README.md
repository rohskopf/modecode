To run all the scripts here, simply do:

    chmod +x run.sh
    ./run.sh

The explanations of these scripts are below.

First go into the files/ directory and put the pair_tep files in your lammps/src directory. Recompile lammps and proceed.

We'll do the harmonic calculation first. Go inside `in.run` and ensure that `compute mode/hf 2` is declared so that we are doing a 2nd order calculation.

First run the MD simulation with harmonic heat flow components:

    lmp_mpi < in.run2

Time-integrate the power transfer calculated by SMCC2s:

    python integrate3.py
    mv integrated_power_mode.dat integrated_power_harmonic.dat

Now run the simulation with 2nd+3rd order heat flow components:

    lmp_mpi < in.run3

Calculate energy change of side B by TEP:

    python extract_de.py

Time-integrate the power transfer calculated by SMCC2s + SMCC3s:

    python integrate3.py
    mv integrated_power_mode.dat integrated_power_anharmonic.dat

Compare the plots:

    gnuplot < gnuplot_echange_compare
