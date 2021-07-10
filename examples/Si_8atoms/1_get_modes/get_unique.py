"""
Open MCC3 file obtained by extract(), and get the unique MCC3s (n' < n'').
The MCC3 file obtained by extract() is used for plotting, and the unique values will be used
for mode/fv in LAMMPS.
"""

import numpy as np

mcc3 = np.loadtxt("MCC3_10")
print(np.shape(mcc3))

boole = mcc3[:,1] <= mcc3[:,2]
indices = np.where(boole)[0]
print(np.shape(indices))

unique = mcc3[indices,:]
unique[:,0:3] = unique[:,0:3].astype(int)
#print(unique)
print(np.shape(unique))
np.savetxt("MCC3",unique, fmt = "%d %d %d %e")
