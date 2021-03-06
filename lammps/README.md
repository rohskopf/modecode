# ModeCode LAMMPS classes.

Put these files into the /src directory of LAMMPS and then install LAMMPS.

## `compute mode`
Compute mode computes and outputs mode amplitudes, velocities, temperature, and energies. In LAMMPS, 
it's used like so:

  compute modeComp all mode 0 # Setting value: 0 - plots, 1 - gifs, 2 - plots and gifs.

where the final integer is an output setting that determines whether to output data for plots or
gifs.

For a setting of 0, data for plots are output where the first column represents the time, and 
the other columns represent the mode indices (starting from zero) for which a quantity (amplitude,
velocity, energy, etc.) is being output. 

For a setting of 1, data for gifs are output where the x-axis of the gif are mode frequencies, and
y-axis are whatever quantity you want to view as a function of time. This option allows one to view
all mode quantities as a function of time.

For a setting of 2, both plot and gif data are output. 

If you use setting 0 or 2, you must supply an INDICES file which labels the mode indices (starting
from zero) that you wish you output data for. 

For a setting of 3, the "mode flux" Jnm = Xn*Vn is output at each time according to the mode pairs given in the PAIRS file.

## `fix mode`
Fix mode sets the initial mode temperatures and coordinates. 

These are set via the MODETEMP and MODECOOR files. For any modes that are not found in these files,
they'll be given a value of zero. 

## `compute mode/fv`
Calculates modal F*V according to the MCC3s in MCC3. The output is fv.dat, where each row is all
F*V ordered in the same way that MCC3 is ordered. 

## `compute mode/hf`
Calculates mode heat flow between two regions of space, which are defined inside the mode coupling
constants or "spatial mode coupling constants".

## `compute mode/heatflux`
Calculates Hardy's heat flux in mode coordinates, with generalized velocities as the inputs. Use like so:

    compute modeHeatflux all mode/heatflux ${nthermo} ${nsteps} ${modeheatflux_setting}
    
where output files are:
- `qtot.dat` - total harmonic heat flux as a function of time.
- `qnm.dat` - time averaged <Qnm> for all mode pairs.

This evalues the expression (Knm)(xn)(vm)

## `compute mode/heatflux2`
Same as compute mode/heatflux, except calculates the heat transfer between two regions of space. Use like so:

    compute modeHeatflux all mode/heatflux2 ${nthermo} ${nsteps} ${modeheatflux_setting}
    
This evaluates the expression (Knm)(xm)(vn) *Note that the xm and vn are switched between here and compute mode/heatflux, this is something in my derivations for QAB.*

## `fix wp`
Creates a longitudinal wave-packet at finite temperature. 

    fix fixname all wp amplitude kmagnitude width velocity z0 a0 direction
    
E.g. for crystalline silicon,

    fix fixWp all wp 5e-2 0.91 30.0 67.0 50..0 5.431 2

Launches the wave-packet in the z-direction (hence the "2" at the last argument). 
