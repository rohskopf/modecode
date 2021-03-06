# Simple heat flux example.

This example involves calculating the heat flux using compute_mode_heatflux.cpp in LAMMPS. This example is mainly for debugging and testing purposes, when developing new functionalities to compute mode/heatflux.

    chmod +x run_ensembles.sh
    ./run_ensembles.sh
    
This example does a 100 K simulation, where we first equilibrate the system at a temperature of 100 K.

*Make sure to edit in.run_template to have the correct variables like temperature and number of steps.*

See the outputs in qtot_house, and verify that it agrees with one or many procs.

qtot is the sum of all mode-mode contributions to the heat flux, and we calculate it as a function of time here. 

So qtot.dat is the total harmonic heat flux as a function of time.

This code also outputs qnm.dat, which is the time average of all mode heat fluxes, <Qnm> listed by their frequencies. This is used for plotting in GNUplot, and is best used with NEMD simulations.

*Keep in mind when using this code for development/debugging purposes, that qnm.dat has values output for which n2>n1, so be aware of that.*
