"""
Integate data in fv.dat
"""

import numpy as np

dx = 0.01
data = np.loadtxt('fv.dat')
mcc3 = np.loadtxt('../MCC3')
nmcc3 = np.shape(mcc3)[0]

nsteps = np.shape(data)[0]
print(nsteps)

# Get the initial energy e0
emdat = np.loadtxt("em.dat")
e0 = emdat[0,1]
print(e0)

integrated = []
#nint = 10 # Number of integrations to perform that splits time interval.
nint = np.linspace(1,nsteps,20)
nint = nint.astype(int)
fh_t = open("times.dat",'w')
fh_sum = open("integrated_fv_sum.dat", 'w')
for t in nint:
  integrated = []
  indx = int(nsteps/t)
  indx = t
  print("%d %f" % (indx-1,data[indx-1,0]))
  fh_t.write("%f\n" % (data[indx-1,0]))
  for s in range(0,nmcc3):
    integrated.append( np.trapz(data[0:indx,s+1],x=data[0:indx,0]) )
  integrated = np.array(integrated)
  sum_integrated = np.sum(integrated)
  np.savetxt("integrated_fv" + str(data[indx-1,0]) + ".dat", integrated)
  fh_sum.write("%f %f\n" % (data[indx-1,0], e0 + (1/2)*sum_integrated))

print("Done integrating!")
fh_t.close()
fh_sum.close()
    

