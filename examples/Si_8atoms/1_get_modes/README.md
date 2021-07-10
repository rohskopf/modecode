First calculate the 2nd order force constants:

    mpirun -np 4 modecode fd 0.001 3.0 1e-8 2 > outfile_ifc2

Perform acoustic sum rule correction:

    modecode asr 2

Calculate the modes:

    python calc_eig.py

Calculate the 3rd order force constants:

    mpirun -np 4 modecode fd 0.001 3.0 1e-10 3 > outfile_ifc3

Calculate MCC3s:

    mkdir mcc3
    mpirun -np 4 modecode ifc2mcc 0 3

Extract all MCC3s (a symmetric matrix) for modes in a range:

    #modecode ifc2mcc 1 3 23 # if you wanna get all of them.
    modecode ifc2mcc 1 10 10 # to get MCC3s for mode 10.

Rename MCC3 file for next step:

    mv MCC3 MCC3_10

Get the unique MCC3s (e.g. ijk is the same as ikj), for use in MD calculations:

    python get_unique.py

This creates a new file called MCC3, which will be used in the LAMMPS calculations next.

