"""
Extract potential and kinetic enregy from log.lammps, and add them to get total energy of each side.
"""

fh = open('log.lammps','r')
fh_e = open('energy.dat','w')
fh_f = open('flux.dat','w')

ecoh = -1.8521648257e+03

fh_e.write('\n')
fh_e.write('\n')
fh_f.write('\n')
fh_f.write('\n')

line = fh.readline()
while ("Step c_pe_A c_ke_A c_pe_B c_ke_B c_fluxA[3] c_fluxB[3]" not in line):
    line = fh.readline()

#line = fh.readline()
e_A = []
e_B = []
hf_A = []
hf_B = []
for l in range(0,101):
    line = fh.readline()
    line_split = line.split()
    step = int(line_split[0])
    pe_A = float(line_split[1])
    ke_A = float(line_split[2])
    pe_B = float(line_split[3])
    ke_B = float(line_split[4])
    flux_A = float(line_split[5])
    flux_B = float(line_split[6])
   
    e_A.append(pe_A+ke_A)
    e_B.append(pe_B+ke_B)
    hf_A.append(flux_A)
    hf_B.append(flux_B) 

    etot_A = pe_A+ke_A - ecoh
    etot_B = pe_B+ke_B - ecoh
    if (l==0):
        eref = etot_A
        fref = flux_A
    ratio_A = etot_A/eref
    ratio_B = etot_B/eref
    ratio_A_flux = flux_A/fref
    ratio_B_flux = flux_B/fref
    fh_e.write('%d %e %e\n' % (step,ratio_A,ratio_B))
    fh_f.write('%d %e %e\n' % (step,ratio_A_flux,ratio_B_flux))


    
    

fh.close()
fh_e.close()
fh_f.close()
