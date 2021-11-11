"""
Read a LAMMPS input script and insert a random number.
"""
import sys

print(sys.argv[1])

fh_r = open("in.run_template", 'r')
fh_w = open("in.run", 'w')

for line_r in fh_r:

  if ("~" in line_r):
    line_r = line_r.replace("~",sys.argv[1])
    print(line_r)
  fh_w.write(line_r)

fh_r.close()
fh_w.close()
