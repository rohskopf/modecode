#Heat flux of Si with a Ge impurity

Relax with 

    lmp_mpi < in.relax

Then copy and paste the atoms from DATA_RELAXED into POSITIONS and sort by id with

    python sort_positions_id.py

This creates POSITIONS_RELAXED_SORTED, which can then be copied and pasted into DATA_RELAXED_SORTED.
