# Getting the modes and coupling constants.

First calculate the 2nd order force constants:

    mpirun -np 4 modecode fd 0.001 3.0 1e-8 2 > outfile_ifc2

Perform acoustic sum rule correction:

    modecode asr 2 1e-8

Calculate the modes:

    python calc_eig.py

Calculate spatial MCC2s:

    mkdir smcc2
    mpirun -np 4 modecode ifc2mcc 5 2

Extract all SMCC2s into a single SMCC2 file

    modecode ifc2mcc 7 4 # 7 is the task, and 4 is the number of files we generated since used 4 procs

Move the files EMAT, FREQUENCIES, FC2_ASR, SMCC2.

    mv SMCC2 ../SMCC2
    mv EMAT ../EMAT
    mv FREQUENCIES ../2_simulate_energy_transfer/FREQUENCIES
    mv FC2_ASR ../2_simulate_energy_transfer/FC2

Now proceed to the next step of simulating the energy transfer.

    cd ../2_simulate_energy_transfer
