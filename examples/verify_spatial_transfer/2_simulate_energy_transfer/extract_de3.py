"""
Extract energy change (final - initial) for all timesteps in LAMMPS log file.
"""

import numpy as np

#data = np.loadtxt('md_energy.dat')

fh = open('log.lammps')
line = fh.readline()
while ("Step c_pe_A" not in line):
  line = fh.readline()

dat = []
counter = 0
ke_arr = []
while ("Loop" not in line):
  line = fh.readline()
  if ("Loop" not in line):
    line_split = line.split()
    nums = [float(x) for x in line_split]
    time = 0.5*nums[0]*1e-3
    #pe = nums[3]
    ke_b = nums[4]
    ke = nums[5]
    ke_arr.append(ke)
    
    #if (time==0.):
    #  e0 = pe+ke

    #etot = pe+ke - e0

    dat.append([time,ke_b,0.,0.]) # 2nd column is ke_b
    counter = counter + 1
#print(line)
fh.close()
dat = np.array(dat)
ke_arr = np.array(ke_arr)

print(dat)

#timesteps = data[:,0]
#pe_B = data[:,3]
#ke_B = data[:,4]
#e_B = pe_B+ke_B

#time = np.array([0.5*timesteps*1e-3])
#de = np.array([e_B - e_B[0]])

#print(dat)

#concat = np.concatenate((time.T,de.T),axis=1)

pe_dat = np.loadtxt("pe3.dat")
pe_b = pe_dat[:,2]
pe_arr = pe_dat[:,3]
nrows = np.shape(dat)[0]
ncol = np.shape(dat)[1]
print(nrows)
for i in range(0,nrows):
  dat[i,2]=pe_b[i] # 3rd column is pe_b
  dat[i,3]=pe_b[i]+dat[i,1] # 4th column is total e_b

dat = dat-dat[0,:] # Find the energy change of side B. 
np.savetxt('energies3.dat', dat)

etot = ke_arr+pe_arr
np.savetxt('total_energies.dat',etot)
