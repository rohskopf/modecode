# Simple mode flux example.

This directory houses a workflow example for calculating the ensemble average of the time-integral of the mode flux correlation, using the 8-atom c-Si example.

    chmod +x run_ensembles.sh
    ./run_ensembles.sh > output

See the simulation outputs (jnm.dat files) in jnm_house, which houses Jnm(t) for the simulation.

*Units of mode flux are in (g/mol)A^2/ps*

To post process the data and view the averages, do:

    cd jnm_house
    ./postproc.sh

## Files

`cor_fft.cpp` calculates cross correlation of jnm$i.dat columns, where each column is calculated against the last column, and then the final column is calculated against itself (autocorrelation). The columns in jnm$i.dat are organized so that the first column is the first pair in files/PAIRS, second column is second pair, and so forth.
`cor_fft2.cpp` calculates autocorrelation of each column.
`integrate_ave` time-integrates the outputs of cor_fft codes, also produces a running average plot (smoother to look at).
`average_ensembles.py` takes the average of time-integrated correlations.
`gnuplot_script` plots the time-integrated correlation functions.
