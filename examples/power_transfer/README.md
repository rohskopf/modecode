# Simple energy transfer example.

This example involves energy transfer from a single excited mode to other modes in crystalline silicon. The purpose here is to show that the ensemble average of mode energy agrees with the time-integrated average of the mode power transfer. To run this example and produce a plot, first install the files in lammps_files/ in your LAMMPS src/ directory. Then do

    chmod +x run_ensembles.sh
    ./run_ensembles.sh

This produces a plot in e1/em_compare.png, showing the agreement between the two quantities. 

The inputs in this simulation, like the FC files and MCC files, are obtained from the Si_8atom example.
