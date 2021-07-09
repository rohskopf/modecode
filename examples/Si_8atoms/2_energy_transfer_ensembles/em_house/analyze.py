import numpy as np

em_list = []
summ_list = []
for i in range(1,10+1):

  em_name = "em%d.dat" % (i)
  summ_name = "summ%d.dat" % (i)
  em = np.loadtxt(em_name)
  summ = np.loadtxt(summ_name)
  em_list.append(em)
  summ_list.append(summ)

em_list=np.array(em_list)
summ_list=np.array(summ_list)

em_avg = np.average(em_list,axis=0)
summ_avg = np.average(summ_list,axis=0)

np.savetxt("em_avg.dat", em_avg)
np.savetxt("summ_avg.dat", summ_avg)
  
