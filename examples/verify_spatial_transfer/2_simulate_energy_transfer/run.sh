lmp_mpi < in.run2
python integrate3.py
mv integrated_power_mode.dat integrated_power_harmonic.dat
lmp_mpi < in.run3
python extract_de.py
python integrate3.py
mv integrated_power_mode.dat integrated_power_anharmonic.dat
gnuplot < gnuplot_echange_compare
