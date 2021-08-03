"""
Sort positions by ID
"""

import numpy as np

fh = open("DATA_RELAXED", 'r')
line = fh.readline()
line = fh.readline()
line = fh.readline()
line_split = line.split()
natoms = int(line_split[0])
print("%d atoms" % (natoms))

dat = np.loadtxt("DATA_RELAXED", skiprows=16, max_rows=natoms)

ids = dat[:,0]

indices = np.argsort(ids)

sorted_dat = dat[indices,:5]

print(np.shape(sorted_dat))

np.savetxt("POSITIONS_RELAXED_SORTED", sorted_dat, fmt='%d %d %.15f %.15f %.15f')
