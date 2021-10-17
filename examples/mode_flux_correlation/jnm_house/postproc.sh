N=5 # Number of ensembles
for ((i=1; i<=$N; i++)); do
  echo $i
  ./calc $i
  ./integrate_ave $i
done

python average_ensembles.py

gnuplot < gnuplot_script
gnuplot < gnuplot_script_ave
