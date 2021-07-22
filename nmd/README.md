# Normal mode dynamics (NMD) simulator.

Go into src/ and compile with

    make

Then go into Si_8atoms to do simulations that compare real-space and mode-space integration. The em3.cpp file in src/ calls the calcModeEnergy() function in mode.cpp, which is used for comparing the validity of the mode power transfer expressions. Simply comment out this line if you want to skip verifying the power transfer.
