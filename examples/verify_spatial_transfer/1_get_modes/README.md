# Getting the modes and coupling constants.

First calculate the 2nd order force constants:

    mpirun -np 4 modecode fd 0.001 3.0 1e-8 2 > outfile_ifc2

Perform acoustic sum rule correction:

    modecode asr 2 1e-8

Calculate the modes:

    python calc_eig.py

Calculate the 3rd order force constants:

    mpirun -np 4 modecode fd 0.001 3.0 1e-10 3 > outfile_ifc3

Calculate spatial MCC2s:

    mkdir smcc2
    mpirun -np 4 modecode ifc2mcc 5 2

Extract all SMCC2s into a single SMCC2 file

    modecode ifc2mcc 7 4 # 7 is the task, and 4 is the number of files we generated since used 4 procs

Calculate spatial MCC3s and move file to SMCC3:

    mkdir mcc3
    mpirun -np 1 modecode ifc2mcc 5 3
    mv mcc3/MCC3_0 SMCC3

Move the files EMAT, FREQUENCIES, FC2_ASR, FC3, SMCC2, SMCC3.

    mv SMCC2 ../SMCC2
    mv SMCC3 ../SMCC3
    mv EMAT ../EMAT
    mv FREQUENCIES ../2_simulate_energy_transfer/FREQUENCIES
    mv FC2_ASR ../2_simulate_energy_transfer/FC2
    mv FC3 ../2_simulate_energy_transfer/FC3

Now proceed to the next step of simulating the energy transfer.

    cd ../2_simulate_energy_transfer
