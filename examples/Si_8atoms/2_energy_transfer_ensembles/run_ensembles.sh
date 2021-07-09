source ~/venv3.6/bin/activate
mkdir em_house
for ((i=1; i<=10; i++)); do
    cd e$i
    lmp_mpi < in.run
    python integrate6.py
    cp em.dat ../em_house/em${i}.dat
    cp integrated_fv_sum.dat ../em_house/summ${i}.dat
    cd ..
done
cd em_house
python analyze.py
gnuplot < gnuplot_em_compare
