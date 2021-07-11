"""
Sort positions by ID
"""

import numpy as np

dat = np.loadtxt("POSITIONS_RELAXED")

ids = dat[:,0]

indices = np.argsort(ids)

sorted_dat = dat[indices,:5]

print(np.shape(sorted_dat))

np.savetxt("POSITIONS_RELAXED_SORTED", sorted_dat, fmt='%d %d %.15f %.15f %.15f')
