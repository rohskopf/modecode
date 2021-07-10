source ~/venv3.6/bin/activate
mkdir em_house
N=10 # Number of ensembles
for ((i=1; i<=$N; i++)); do
    cp -r files e$i
    cd e$i
    rando_var=$((1000*i + 10101))
    python rando_vel.py ${rando_var}
    lmp_mpi < in.run
    python integrate6.py
    cp em.dat ../em_house/em${i}.dat
    cp integrated_fv_sum.dat ../em_house/summ${i}.dat
    cd ..
    #rm -r e$i
done
cd em_house
python analyze.py $N
gnuplot < gnuplot_em_compare
