# Edit cor_fft2.cpp and integrate_ave.cpp to have the desired settings like Nmodes (number of columns or pairs), then compile.

g++ cor_fft.cpp -lfftw3 -o calc
g++ cor_fft2.cpp -lfftw3 -o calc2
g++ integrate_ave.cpp -o integrate_ave

N=5 # Number of ensembles
for ((i=1; i<=$N; i++)); do
  echo $i
  ./calc $i
  ./integrate_ave $i
done

python average_ensembles.py $N 100

#gnuplot < gnuplot_script
#gnuplot < gnuplot_script_ave
