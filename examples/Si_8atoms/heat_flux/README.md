# Simple heat flux example.

This example involves calculating the heat flux using compute_mode_heatflux.cpp in LAMMPS. This example is mainly for debugging and testing purposes, when developing new functionalities to compute mode/heatflux.

    chmod +x run_ensembles.sh
    ./run_ensembles.sh

See the outputs in qtot_house, and verify that it agrees with one or many procs.

qtot is the sum of all mode-mode contributions to the heat flux, and we calculate it as a function of time here. 

So qtot.dat is the total harmonic heat flux as a function of time.
