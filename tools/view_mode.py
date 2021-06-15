import numpy as np

mode = 20

amp = [0.,100.]

ndisps = 10

emat = np.loadtxt('eigenvectors.csv', delimiter=',')

shape = np.shape(emat)

natoms = shape[0]/3

eindx = 0

# Read evecs for mode of interest

evecs = []
for n in range(0,natoms):
    evec = []
    for a in range(0,3):
        evec.append(emat[eindx+a][mode])
    eindx += 3
    evecs.append(evec)

evecs = np.array(evecs)

# Read config

fh_config = open('CONFIG','r')

line = fh_config.readline()
line = fh_config.readline()
box = [float(a) for a in line.split()]
x = []
types = []
for n in range(0,natoms):
    line = fh_config.readline()
    line_split = [float(a) for a in line.split()]
    x.append(line_split[1:])
    types.append(int(line_split[0]))

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


fh_config.close()
fh_xyz.close() 


