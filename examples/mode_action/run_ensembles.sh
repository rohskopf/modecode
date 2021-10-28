source ~/venv-3.8.10/bin/activate
#mkdir em_house
#mkdir jnm_house
N=20 # Number of ensembles
temperature=300
for ((i=1; i<=$N; i++)); do
    echo $i
    rm -r e$i
    cp -r files e$i
    cd e$i
    rando_var=$((1000*i + 10101))
    python rando_vel.py ${rando_var}
    lmp_mpi < in.run
    #python integrate6.py
    #cp em.dat ../em_house/em${i}.dat
    cp jnm.dat ../jnm_house/jnm${i}.dat
    cd ..
    #rm -r e$i
done

cd jnm_house
for ((i=1; i<=$N; i++)); do
  echo $i
  #./calc2 $i
  ./calc $i
  ./integrate_ave $i
done

python average_ensembles.py $N ${temperature}
gnuplot < gnuplot_script
gnuplot < gnuplot_script_ave
#gnuplot < gnuplot_script_tnm
#mv plot.png plot_${N}ensembles.png

#cd em_house
#python analyze.py $N
#gnuplot < gnuplot_em_compare
