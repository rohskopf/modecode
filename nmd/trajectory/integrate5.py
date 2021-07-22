"""
Integate data in fv.dat
"""

import numpy as np

#dx = 0.01
#nmcc3 = 3 # This number comes from output of running in.run.
data = np.loadtxt('fv.dat')
em = np.loadtxt('em.dat')[0,1]

nsteps = np.shape(data)[0]
print(nsteps)

integrated = []
#nint = 10 # Number of integrations to perform that splits time interval.
nint = np.linspace(1,nsteps,20)
nint = nint.astype(int)
#fh_t = open("times.dat",'w')
#fh_sum = open("integrated_fv.dat", 'w')
integrated = []
for t in nint:
  indx = int(nsteps/t)
  indx = t
  print(indx)
  print("%d %f" % (indx-1,data[indx-1,0]))
  #print([data[indx-1,0], np.trapz(data[0:indx,1],x=data[0:indx,0])])
  integrated.append([data[indx-1,0], em + np.trapz(data[0:indx,1],x=data[0:indx,0]) ] )

integrated = np.array(integrated)
#print(integrated)
np.savetxt("integrated_fv.dat", integrated)

print("Done integrating!")
    

