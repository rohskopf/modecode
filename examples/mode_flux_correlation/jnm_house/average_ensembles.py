import numpy as np

nens = 5 # number of ensembles

i = 1
filename = "area_evol_CCs%d.txt" % (i)
filename_ave = "running_ave%d.txt" % (i)
data = np.loadtxt(filename)
data_ave = np.loadtxt(filename_ave)

shape = np.shape(data)
shape_ave = np.shape(data_ave)
#print(shape)

nrows = shape[0]
nrows_ave = shape_ave[0]

summ = np.zeros(shape)
summ_ave = np.zeros(shape_ave)

for i in range(1, nens+1):

  filename = "area_evol_CCs%d.txt" % (i)
  filename_ave = "running_ave%d.txt" % (i)
  data = np.loadtxt(filename)
  data_ave = np.loadtxt(filename_ave)

  summ = summ + data
  summ_ave = summ_ave + data_ave

average = np.divide(summ, nens)
average_ave = np.divide(summ_ave, nens)

np.savetxt("ensemble_average.dat", average)
np.savetxt("ensemble_average_ave.dat", average_ave)




