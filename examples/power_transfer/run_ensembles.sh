source ~/venv3.6/bin/activate
mkdir em_house
N=1 # Number of ensembles
for ((i=1; i<=$N; i++)); do
    rm -r e$i
    cp -r files e$i
    cd e$i
    rando_var=$((30100 + i))
    python rando_vel.py ${rando_var}
    lmp_mpi < in.run
    python integrate5.py
    gnuplot < gnuplot_em_compare
    cd ..
    #rm -r e$i
done
