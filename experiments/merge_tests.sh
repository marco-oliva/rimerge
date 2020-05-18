#!/bin/bash
#SBATCH --job-name=merge                # Job name
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marco.oliva@ufl.edu # Where to send mail
#SBATCH --ntasks=8                      # Number of MPI ranks
#SBATCH --cpus-per-task=16              # Number of cores per MPI rank
#SBATCH --nodes=5                       # Number of nodes
#SBATCH --ntasks-per-node=2             # How many tasks on each node
#SBATCH --ntasks-per-socket=2           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=100mb             # Memory per core
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=%j.log                 # Standard output and error log
#SBATCH --constraint='intel&infiniband'

##----------------------------------------------------------
# Modules
module load gcc/8.2.0
module load python

##----------------------------------------------------------
# Vars
base_path="/ufrc/boucher/marco.oliva/data/r-index/Chr19"
PATH=$PATH:"/home/marco.oliva/bin"

#-----------------------------------------------------------


${base_path}/main --version

/usr/bin/time --verbose ${base_path}/main -a ${base_path}/chr19.1000_s.00.seqs -b ${base_path}/chr19.1000_s.01.seqs -o ${base_path}/chr19.1000_s.out.seqs
for i in {2..9}
do
    /usr/bin/time --verbose ${base_path}/main -a ${base_path}/chr19.1000_s.out.seqs -b ${base_path}/chr19.1000_s.0${i}.seqs -o ${base_path}/chr19.1000_s.out.seqs
done

date