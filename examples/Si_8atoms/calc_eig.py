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


fh = open('HESSIAN', 'r')

natoms = 8

line = fh.readline()

# Read ALM hessian
nav = 6.02e23

mws1 = 28.0855 
mws2 = 16.
mws = np.zeros(natoms)
mws[0:]=mws1
#mws[183:]=mws2
masses = (mws/nav)*1e-3
#print(masses)

#prefactor = ((2.1798741*10**-18)/(5.29177*10**-11)**2)/(28./(6.02*10**23)) # SI

j_over_r = 2.179874099e-18
b_over_m = 1/5.29177e-11
nav = 6.02e23

#convert = j_over_r*(b_over_m)**2 * 1e3 * nav * mw
convert = j_over_r*(b_over_m)**2 # (J/m^2)=(Ryd/bohr^2)*(J/Ryd)*(bohr/m)^2

#print(prefactor)

hessian_alm = np.zeros([3*natoms,3*natoms,3,3])
for i in range(0,natoms):
    for a in range(0,3):
        for j in range(0,natoms):
            for b in range(0,3):
                line = fh.readline()
                line_split = line.split()
                ii = int(line_split[0])-1
                aa = int(line_split[1])-1
                jj = int(line_split[2])-1
                bb = int(line_split[3])-1
                fc = float(line_split[4])
                #print("%d %d %d %d\n" % (ii,aa,jj,bb))
                hessian_alm[ii][jj][aa][bb] = fc*convert/(np.sqrt(masses[i]*masses[j]))

#hessian_alm = hessian_alm * convert


# Make properly indexed Hessian

hessian = np.zeros([3*natoms, 3*natoms])
ii = 0
jj = 0
for i in range(0,natoms):
    ii = 0
    for j in range(0,natoms):
        for a in range(0,3):
            for b in range(0,3):

                #print('ii+a: %d, jj+b: %d' % (ii+a,jj+b))
                hessian[ii+a][jj+b] = hessian_alm[i][j][a][b]

        ii += 3
    jj += 3
    
fh.close()

#print(is_diagonally_dominant(hessian))

#np.savetxt('dynmat.csv', hessian, delimiter=',')

#w, v = LA.eig(hessian)
w, v = LA.eigh(hessian)

w = np.lib.scimath.sqrt(w)
w = w/(2*np.pi)

print(w)

#w = np.sqrt(w)/(1e12*2.*np.pi)
#w = np.nan_to_num(w)

np.savetxt('frequencies.txt', np.real(w), delimiter='\n')

#plt.hist(w, bins=50)  # arguments are passed to np.histogram
#plt.show()

#print(v)

vsq = v**2
vsqsum = np.sum(vsq, axis=1) # should sum to 1
#print(" Norm of eigenvectors: %f" % (vsqsum))
print(vsqsum)


# Check orthogonality of all modes

"""
dots = []
for m1 in range(0,216*3):
    print('m1: %d' % (m1))
    for m2 in range(0,216*3):
        eindx = 0
        dotsum = 0.
        for n in range(0,216):
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

        
                

np.savetxt('eigenvectors.csv', v, delimiter=',')
np.savetxt('eigenvectors.txt', v, delimiter=' ')

"""
#print(np.shape(v))

ncount = 0
for n in range(0,natoms):
    for a in range(0,3):
        print(v[0][ncount+a])
    ncount += 3
    print('\n')
"""

x = np.arange(0,natoms*3)
w = w/1e12
plt.scatter(x, w)
plt.ylabel("Frequency (THz)")
plt.xlabel("Mode")
plt.show()
