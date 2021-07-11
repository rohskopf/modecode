"""
Sort the mode-mode heat fluxes by magnitude so that the largest get plotted last in gnuplot.
"""

import numpy as np

dat = np.loadtxt("fv.dat")

indices = np.argsort(dat[:,2])
#sort = np.sort(dat[:,2])
freq1 = dat[indices,0]
freq2 = dat[indices,1]
vals = dat[indices,2]

freq1 = np.array([freq1]).T
freq2 = np.array([freq2]).T
vals = np.array([vals]).T

cat = np.concatenate((freq1,freq2,vals), axis=1)

print(cat)

np.savetxt("fv_sorted.dat", cat)
