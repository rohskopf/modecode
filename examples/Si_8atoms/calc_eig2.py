import numpy as np
from numpy import linalg as LA
from scipy.linalg import eigh
import matplotlib.pyplot as plt

def dd(X):
    D = np.diag(np.abs(X)) # Find diagonal coefficients
    S = np.sum(np.abs(X), axis=1) - D # Find row sum without diagonal
    if np.all(D > S):
        print('matrix is diagonally dominant')
    else:
        print('NOT diagonally dominant')
    return

def is_diagonally_dominant(x):
    abs_x = np.abs(x)
    return np.all( 2*np.diag(abs_x) >= np.sum(abs_x, axis=1) )


#fh = open('FC2_FD, 'r')

natoms = 8
#mass = 4.6637066e-26
mws1 = 28.0855 # Si
mws2 = 15.999 # O
mws = np.zeros(natoms)
mws[0:natoms]=mws1
#mws[184:328]=mws2
# Convert mws [kg/kmol] to kg
nav = 6.0221409e+23 # number/mol
mass = mws/(nav*1e3)
#print(mass)


#line = fh.readline()

dm = np.zeros([natoms*3,natoms*3]) # dynamical matrix

nfc2 = -1
with open('FC2') as fh:
    for line in fh:
        nfc2 = nfc2+1
        if (nfc2 != 0):
            #print line,  # The comma to suppress the extra new line char
            line_split = line.split()
            i = int(line_split[0])-1
            a = int(line_split[1])-1
            j = int(line_split[2])-1
            b = int(line_split[3])-1
            fc = float(line_split[4])*(2.179874099E-18)*1.89e+10*1.89e+10 # Convert Ryd/Bohr^2 to J/m^2
            #print(fc)
            dm[3*i+a][3*j+b]=fc/(np.sqrt(mass[i]*mass[j]))

w, v = LA.eigh(dm)

#w = np.lib.scimath.sqrt(w)
#w = w/(2*np.pi)

print(w)

#np.savetxt('MCC2', np.real(w), delimiter='\n')

#plt.hist(w, bins=50)  # arguments are passed to np.histogram
#plt.show()

#print(v)

vsq = v**2
vsqsum = np.sum(vsq, axis=1) # should sum to 1
#print(" Norm of eigenvectors: %f" % (vsqsum))
print(vsqsum)    

np.savetxt('EMAT', v, delimiter=' ')

# Check orthogonality
"""
dots = []
for m1 in range(0,natoms*3):
    print('m1: %d' % (m1))
    for m2 in range(0,natoms*3):
        eindx = 0
        dotsum = 0.
        for n in range(0,natoms):
            v1 = []
            v2 = []
            for a in range(0,3):
                v1.append(v[eindx+a][m1])
                v2.append(v[eindx+a][m2])
            #print(v1)
            #print(v2)
            dot = np.dot(v1,v2)
            #print(dot)
            #print('----------')
            eindx +=3
            dotsum += dot
            
        #print(dotsum)
        if (m1 == m2 and np.round(dotsum,6) != 1.):
            print('%d %d %f' % (m1,m2,dotsum))
        if (m1 != m2 and np.round(dotsum,6) != 0.):
            print('%d %d %f' % (m1,m2,dotsum))

        if (m1 != m2):
            dots.append(dotsum)

dots = np.array(dots)
print(np.mean(dots))
print(np.max(dots))
"""

# Calculate MCC2 manually
"""
fh = open("MCC2_PYTHON", 'w')
e=v
for n1 in range(3,natoms*3):
    for n2 in range(3,natoms*3):
        mcc2 = 0.0
        for i in range(0,natoms):
            for j in range(0,natoms):
                for a in range(0,3):
                    for b in range(0,3):
                        mcc2 = mcc2 + dm[3*i+a][3*j+b]*e[3*i+a][n1]*e[3*j+b][n2]
        fh.write("%d %d %e\n" % (n1,n2,np.sqrt(mcc2/(2*np.pi))))
fh.close()
"""
e_inv = np.linalg.inv(v)
e = v
first = np.matmul(e_inv,dm)
omega = np.matmul(first,e)
#omega = np.sqrt(omega/(2.*np.pi))
#np.savetxt("OMEGA", omega)
fh = open("MCC2_PYTHON", 'w')
for i in range(3,3*natoms):
    for j in range(3,3*natoms):
        if (i==j):
            fh.write("%d %d %e\n" % (i,j,omega[i][j]))
fh.close()
        
                        

x = np.arange(0,natoms*3)
w = np.lib.scimath.sqrt(w)
w = w/(2*np.pi)
w = w/1e12
plt.scatter(x, w)
plt.ylabel("Frequency (THz)")
plt.xlabel("Mode")
plt.show()
