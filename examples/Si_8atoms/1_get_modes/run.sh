#!/bin/bash
#SBATCH -N 30
##SBATCH -n 1
#SBATCH --time=5:00:00
#SBATCH --constraint=centos7
#SBATCH --mem=186G
#SBATCH --ntasks-per-node=40
#SBATCH -p sched_mit_ase
#SBATCH -o job-%x.out 
#SBATCH -e job-%x.err
#SBATCH -J fd2  

module purge
module load openmpi/gcc/64/1.8.1
module load gcc/8.3.0

export PATH=$PATH:/home/rohskopf/fcfd/src
source /home/rohskopf/venv3.6.10/bin/activate



# WORK FLOW FOR MAKING KNM_AB GRID
#####################
# 1. Calc FC2
mpirun -np $SLURM_NTASKS modecode fd 0.001 3.0 2 > outfile_ifc2
# 2. ASR & eigenvector
#fcfd asr 2
#python calc_eig3.py
# 3. Calculate SMCC2
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 5 2 > outfile_smcc2
# 4. Get SMCC2, organize, and sum.
#fcfd ifc2mcc 7 400 # 7 - setting, next number is number of files in /smcc2
#python organize_smcc2.py
#python sum_smcc2.py

#####################


#mpirun -np $SLURM_NTASKS fcfd fd 0.01 2 > outfile_ifc2
#fcfd asr 2
#python calc_eig3.py
#mpirun -np $SLURM_NTASKS fcfd fd 0.01 3 > outfile_ifc3
#mpirun -np $SLURM_NTASKS fcfd fd 0.01 4 > outfile_ifc4

#fcfd ifc2mcc 2 0 1600
##fcfd ifc2mcc 3

# Compute PR, last entry is number of atoms
#fcfd compute pr 4320

# Visualize mode
#python view_mode.py

# Gotta make these directories to run!
#mkdir debug_ifc2mcc
#mkdir mcc3
#mkdir smcc2
##mkdir mcc4
#mpirun -np 1600 fcfd ifc2mcc 3 > outfile

# Extract MCC3s associated with modes i in a range (Option 1)
# Need to edit "nfiles" variable in ifc2mcc.cpp extract function.
#fcfd ifc2mcc 1 9834 9834 # Extracts MCC3s associated with these modes i in a range.
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 1 0 2399 > outfile_extract 

#fcfd ifc2mcc 4 3
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 0 3 > outfile_ifc2mcc_3rd

#SPATIAL MCC2 CALCULATION.
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 5 2 > outfile_smcc2

# Extract SMCC2 into SMCC2 file
#fcfd ifc2mcc 7 400 # 7 - setting, next number is number of files in /smcc2

#Organize SMCC2
#python organize_smcc2.py

#Plot SMCC2
#gnuplot < gnuplot_smcc2

# SPATIAL MCC3 CALCULATION.
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 5 3 > outfile_smcc3

# Average MCC2s
#python avg_mcc2.py

# Single mode (n1) MCC3 calculation.
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 6 9834 > outfile_mcc3_9834

# Calculate TIC
#mpirun -np $SLURM_NTASKS fcfd tic

# Calculate GV and store in "GV" file.
#mkdir gv
#mpirun -np $SLURM_NTASKS fcfd ifc2mcc 8 2

# Calculate TC with QHGK.
#mpirun -np $SLURM_NTASKS fcfd qhgk 0 75

