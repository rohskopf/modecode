"""
Calculate trajectory MPE between two files
"""

import numpy as np
import matplotlib.pyplot as plt

file1 = "dump_atom.xyz" # Calculate MPE wrt this file
#file2 = "dump_mode.xyz" # Calculate MPE of this file
#file2 = "dump_atom_integrate0.xyz" # Calculate MPE of this file
file2 = "dump_atom_harmonic.xyz"
natoms = 8
ndat = 10000 # number timesteps that data was dumped

fh1 = open(file1, "r")
fh2 = open(file2, "r")

x1 = []
timestep = []
for t in range(0,ndat):
  line = fh1.readline()
  line = fh1.readline()
  line_split = line.split()
  timestep.append(int(line_split[2]))
  x1t = []
  for i in range(0,natoms):
    line = fh1.readline()
    nums = [float(x) for x in line.split()]
    x1t.append([nums[1],nums[2],nums[3]])
  x1.append(x1t)
timestep = np.array(timestep)
timestep = 0.5*timestep
timestep = timestep/1000.
x1 = np.array(x1)

x2 = []
for t in range(0,ndat):
  line = fh2.readline()
  line = fh2.readline()
  x2t = []
  for i in range(0,natoms):
    line = fh2.readline()
    nums = [float(x) for x in line.split()]
    x2t.append([nums[1],nums[2],nums[3]])
  x2.append(x2t)
x2 = np.array(x2)

# Calculate MPE per timestep
mpe = []
for t in range(0,ndat):
  mpet = 0.0
  for n in range(0,natoms):
    diff = x2[t][n]-x1[t][n]
    diffnorm = np.linalg.norm(diff)
    norm = np.linalg.norm(x1[t][n])
    mpet = mpet + diffnorm/norm
  mpet = mpet/natoms
  mpe.append(mpet)
mpe = np.array(mpe)
mpe = mpe*100.

# print MPE data
fh = open("MPE " + file2, "w")
for t in range(0,ndat):
  fh.write("%f\n" % (mpe[t]))

timestep = np.array([timestep]).T
mpe = np.array([mpe]).T
dat = np.concatenate((timestep,mpe),axis=1)
print(np.shape(dat))

np.savetxt("traj_err.dat", dat)

#plt.plot(timestep,mpe, '-')
#plt.show()


fh1.close()
fh2.close()
