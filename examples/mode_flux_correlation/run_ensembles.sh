source ~/venv-3.8.10/bin/activate
#mkdir em_house
mkdir jnm_house
N=5 # Number of ensembles
for ((i=1; i<=$N; i++)); do
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

#cd em_house
#python analyze.py $N
#gnuplot < gnuplot_em_compare
