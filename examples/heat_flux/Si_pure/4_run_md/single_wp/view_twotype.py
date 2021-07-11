"""
View wave packet
"""
import numpy as np

"""
Read original (not displaced) positions.
"""
data_o = "DATA" # original data (not displaced)
mag = 1. # magnitude factor to view

fh_o = open(data_o,'r')
# Get to 17th line
for i in range(0,16):
    line = fh_o.readline()
    if(i==2):
        line_split = line.split()
        natoms = int(line_split[0])

x0 = []
for i in range(0,natoms):
    print(i)
    line = fh_o.readline()
    #print(line)
    line_split = line.split()
    x = float(line_split[2])
    y = float(line_split[3])
    z = float(line_split[4])
    x0.append([x,y,z])
x0 = np.array(x0)
fh_o.close()

"""
Read dump file.
"""
fh = open('dump.xyz', 'r')
fh_w = open('dump_magnified.xyz','w')
fh_dat = open('disp.dat','w')
images = 100
for t in range(0,images):
    line = fh.readline()
    fh_w.write(line)
    line = fh.readline()
    fh_w.write(line)
    fh_dat.write('\n')
    fh_dat.write('\n')
    for i in range(0,natoms):
        line = fh.readline()
        line_split = line.split()
        typ = int(line_split[0])
        x = float(line_split[1])
        y = float(line_split[2])
        z = float(line_split[3])
        ux = (x - x0[i][0])*mag
        uy = (y - x0[i][1])*mag
        uz = (z - x0[i][2])*mag
        #fh_w.write('%d %.10f %.10f %.10f\n' % (typ,x0[i][0]+ux,x0[i][1]+uy,x0[i][2]+uz))
        fh_dat.write('%e %e\n' % (x0[i][2],uz))
            
fh.close()
fh_w.close()
fh_dat.close()
            

