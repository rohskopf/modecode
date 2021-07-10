"""
Integate data in q.dat
"""

import numpy as np

dx = 0.0005
counter = 0
integral1 = 0.0
integral2 = 0.0
integral3 = 0.0
data = np.loadtxt('ht.dat')
x = data[:,0]
y = data[:,1]
test = np.trapz(y,x)
print("Numpy: %e" % (test))

lendata = len(data)
print("%d timesteps." % (lendata))
integrals = []
for l in range(0,lendata):
  #print(l)
  if (l % 1 == 0):
    print(l)
    intervals = [0.0,data[l,0]]
    integral = 0.0
    for l2 in range(0,lendata):
      time = data[l2,0]

      f = data[l2,1]
      if ( (intervals[0] <= time) and (time <= intervals[1]) ):
        #print(time)
        if ( (time == intervals[0]) or (time == intervals[1]) ):
          coeff = 1.0
        else:
          coeff = 2.0

        integral = integral + 0.5*dx*coeff*f
        if ( time == intervals[1]):
          break
    integrals.append([data[l,0], integral])
integrals = np.array(integrals)
np.savetxt('integrated_power.dat', integrals, delimiter=' ')
    
"""
with open('fv_test.dat') as fh:
    for line in fh:
        line_split = line.split()
        nums = [float(x) for x in line_split]
        #print(nums)
        time = nums[0]
        f1 = nums[1]
        f2 = nums[2]
        f3 = nums[3]
        if ( (intervals[0] <= time) and (time <= intervals[1]) ):
          #print(time)
          if ( (time == intervals[0]) or (time == intervals[1]) ):
            coeff = 1.0
          else:
            coeff = 2.0

          integral1 = integral1 + 0.5*dx*coeff*f1
          integral2 = integral2 + 0.5*dx*coeff*f2
          integral3 = integral3 + 0.5*dx*coeff*f2
          if ( time == intervals[1]):
            break

print(" %e eV" % (integral2))
"""
