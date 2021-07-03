import numpy as np

mode = 1055 # Mode index starting from 1

amp = [0.,100.]

ndisps = 10

emat = np.loadtxt('EMAT')

shape = np.shape(emat)

#print(shape)

natoms = int(shape[0]/3)

eindx = 0

# Read evecs for mode of interest

evecs = []
for n in range(0,natoms):
    evec = []
    for a in range(0,3):
        evec.append(emat[eindx+a][mode-1])
    eindx += 3
    evecs.append(evec)

evecs = np.array(evecs)


#Read LAMMPS data file.
nskip = 0 # Used by np.loadtxt soon.
fh = open("DATA", 'r')
# Ignore first line
line = fh.readline()
nskip+=1
# Next line is empty.
line = fh.readline()
nskip+=1
# Get number of atoms.
line = fh.readline()
nskip+=1
natoms = int(line.split()[0])
print("%d atoms." % (natoms))
# Find line with atom types
while ("types" not in line):
  line = fh.readline()
  nskip+=1
ntypes = int(line.split()[0])
print("%d atom types." % (ntypes))
# Find masses.
while ("Masses" not in line):
  line = fh.readline()
  nskip+=1
# Skip next line.
line = fh.readline()
nskip+=1
# Store masses.
masses_types = []
for t in range(0,ntypes):
  line = fh.readline()
  nskip+=1
  masses_types.append( float(line.split()[1]) )
masses_types = np.array(masses_types)
print("Masses of types:")
print(masses_types)
# Find atoms.
while ("Atoms" not in line):
  line = fh.readline()
  nskip+=1
# Skip next line.
line = fh.readline()
nskip+=1
print("Should skip %d lines to read atoms." % (nskip))
fh.close()

data = np.loadtxt('DATA', skiprows=nskip)
print(data)
types = data[:,1]
print(types)
x = data[:,2:]
print(x)

# Make xyz

fh_xyz = open('mode.xyz', 'w')

amps = np.linspace(amp[0],amp[1],num=ndisps)

timestep = 0
tindx = 0
for d in range(0,ndisps):

    amp = amps[d]
    fh_xyz.write('%d\n' % (natoms))
    timestep = d+tindx
    fh_xyz.write('Atoms. Timestep: %d\n' % (timestep))
    for n in range(0,natoms):

        evec = evecs[n]
        norm = np.linalg.norm(evec)
        mag = amp*norm

        #print(evec)
        disp = evec*mag
        
        xnew = x[n][0]+disp[0]
        ynew = x[n][1]+disp[1]
        znew = x[n][2]+disp[2]
        fh_xyz.write('%d %f %f %f\n' % (types[n], xnew,ynew,znew))

# Go back to zero

tindx += ndisps
flipamps = np.flip(amps,0)[1:]
for d in range(0,ndisps-1):
    
    amp = flipamps[d]
    fh_xyz.write('%d\n' % (natoms))
    timestep = d+tindx
    fh_xyz.write('Atoms. Timestep: %d\n' % (timestep))
    for n in range(0,natoms):

        evec = evecs[n]
        norm = np.linalg.norm(evec)
        mag = amp*norm

        #print(evec)
        disp = evec*mag
        
        xnew = x[n][0]+disp[0]
        ynew = x[n][1]+disp[1]
        znew = x[n][2]+disp[2]
        fh_xyz.write('%d %f %f %f\n' % (types[n], xnew,ynew,znew))
    
# Start negative displacement

tindx += ndisps
for d in range(0,ndisps):
    
    amp = -1.*amps[d]
    fh_xyz.write('%d\n' % (natoms))
    timestep = d+tindx
    fh_xyz.write('Atoms. Timestep: %d\n' % (timestep))
    for n in range(0,natoms):

        evec = evecs[n]
        norm = np.linalg.norm(evec)
        mag = amp*norm

        #print(evec)
        disp = evec*mag
        
        xnew = x[n][0]+disp[0]
        ynew = x[n][1]+disp[1]
        znew = x[n][2]+disp[2]
        fh_xyz.write('%d %f %f %f\n' % (types[n], xnew,ynew,znew))

# Go back to zero

tindx += ndisps
flipamps = np.flip(amps,0)[1:]
for d in range(0,ndisps-1):
    
    amp = -1.*flipamps[d]
    fh_xyz.write('%d\n' % (natoms))
    timestep = d+tindx
    fh_xyz.write('Atoms. Timestep: %d\n' % (timestep))
    for n in range(0,natoms):

        evec = evecs[n]
        norm = np.linalg.norm(evec)
        mag = amp*norm

        #print(evec)
        disp = evec*mag
        
        xnew = x[n][0]+disp[0]
        ynew = x[n][1]+disp[1]
        znew = x[n][2]+disp[2]
        fh_xyz.write('%d %f %f %f\n' % (types[n], xnew,ynew,znew))


#fh_config.close()
fh_xyz.close() 


