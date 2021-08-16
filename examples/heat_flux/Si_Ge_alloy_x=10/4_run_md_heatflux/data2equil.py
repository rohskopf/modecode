"""
Convert LAMMPS data file to EQUIL file.
"""

import numpy as np

data = np.loadtxt("DATA_RELAXED_SORTED", skiprows=17)
data = data[:,2:]
np.savetxt("EQUIL", data)
