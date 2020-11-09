#!/bin/bash
#SBATCH --job-name=filter_chrs          # Job name
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marco.oliva@ufl.edu # Where to send mail
#SBATCH --ntasks=1                   	  # Run a single task
#SBATCH --cpus-per-task=1           	  # Number of CPU cores per task
#SBATCH --mem=4gb                     	# Job memory request
#SBATCH --time=36:00:00              	  # Time limit hrs:min:sec
#SBATCH --output=%j_performances.log    # Standard output and error log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Modules
module load gcc/8.2.0
module load python3

##----------------------------------------------------------
# Vars
profiler="/usr/bin/time"
base_path="/blue/boucher/marco.oliva/data/r-index/1000GP/fasta"
PATH="/home/marco.oliva/bin":$PATH

CHR=15

echo "Run sed --"
sed -i 's/^>[A-Z]/N/g' "${base_path}/Chr${CHR}/Chr${CHR}_2k.fa"

echo "Run filter --"
filter_file.py -i "${base_path}/Chr${CHR}/Chr${CHR}_2k.fa"


##----------------------------------------------------------
# Print Some statistics, again
pwd; hostname; date