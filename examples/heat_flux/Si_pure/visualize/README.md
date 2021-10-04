# Heat flux visualization example.

Here we "visualize" heat flux by observing the time dependent displacements associated with combinations of modes, where each atomic displacement is weighted by the mode coupling constant magnitude.

The `visualize` feature works with the following command:

    mpirun -np nprocs modecode visualize temperature nindx timestep scalefactor
    
where 

- `temperature` is the temperature in Kelvin, which controls displacements.
- `nindx` is the index, starting from 0, of the mode you want to visualize heat flux for.
- `timestep` is the timestep to view the visualization in ps.
- `scalefactor` is scales the atomic displacements for visualization purposes.

To try this example, first place EMAT and GV in the previous directory (these can be downloaded here: https://www.dropbox.com/sh/nq1owq4wmh285qy/AACNYH8kVLq2xQU91AfctHv9a?dl=0).

Then do:

    modecode visualize 300 21 10 1e-1
    
which visualizes mode 21 (a transverse acoustic mode) at a temperature of 300 K, with a timestep of 10 ps, and a 1e-1 scale factor on the displacements.

Then watch the movie of atomic displacements with:

    vmd disp.xyz

The scale factor and timestep should be played with to get a coherent movie.
