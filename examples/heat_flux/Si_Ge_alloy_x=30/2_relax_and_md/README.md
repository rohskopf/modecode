#Heat flux of Si with a Ge impurity

Relax with 

    lmp_mpi < in.relax

This creates DATA_RELAXED. Run this script to sort the atom ids in DATA_RELAXED:

    python sort_positions_id.py

Which makes DATA_RELAXED_SORTED. To check that the forces are zero, run a single MD step with

    lmp_mpi < in.run
