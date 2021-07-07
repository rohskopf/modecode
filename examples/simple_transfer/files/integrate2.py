"""
Integate data in fv.dat
"""

import numpy as np

dx = 0.01
counter = 0
intervals = [0.0,350]
integral1 = 0.0
integral2 = 0.0
integral3 = 0.0
data = np.loadtxt('fv.dat')
lendata = len(data)
print("%d timesteps." % (lendata))
integrals = []
for l in range(0,lendata):
  print(l)
  intervals = [0.0,data[l,0]]
  integral1 = 0.0
  integral2 = 0.0
  integral3 = 0.0
  for l2 in range(0,lendata):
    time = data[l2,0]

    f1 = data[l2,1]
    f2 = data[l2,2]
    f3 = data[l2,3]
    if ( (intervals[0] <= time) and (time <= intervals[1]) ):
      #print(time)
      if ( (time == intervals[0]) or (time == intervals[1]) ):
        coeff = 1.0
      else:
        coeff = 2.0

      integral1 = integral1 + 0.5*dx*coeff*f1
      integral2 = integral2 + 0.5*dx*coeff*f2
      integral3 = integral3 + 0.5*dx*coeff*f3
      if ( time == intervals[1]):
        break
  integrals.append([data[l,0], integral1, integral2, integral3])
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
