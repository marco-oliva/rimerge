#!/bin/bash
#SBATCH --job-name=performance          # Job name
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marco.oliva@ufl.edu # Where to send mail
#SBATCH --ntasks=1                   	  # Run a single task
#SBATCH --cpus-per-task=32           	  # Number of CPU cores per task
#SBATCH --mem=64gb                     	# Job memory request
#SBATCH --time=168:00:00              	# Time limit hrs:min:sec
#SBATCH --output=%j_performances.log    # Standard output and error log
#SBATCH --exclusive
#SBATCH --nodelist=c27b-s7

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
BWT_THREADS=32
MERGE_JOBS=16
EXEC_BWT=false
EXEC_MERGE=true
CHRS=( 15 16 17 18 19 20 )

##----------------------------------------------------------
# Create index of all files
if [ "${EXEC_BWT}" = true ]
then
  echo "================================================================================"
  echo "Create index all files"
  ${profiler} bigbwt -N 12024 -t ${BWT_THREADS} -s -e "${base_path}/Chr_${CHRS[0]}_${CHRS[-1]}.filtered.seqs"
fi

##----------------------------------------------------------
# Create index by merging
if [ "${EXEC_MERGE}" = true ]
then
  echo "================================================================================"
  echo "Create index by merging"
  # Create index of first chromosome and then merge the rest
  rimerge.py -a "${base_path}/Chr${CHRS[0]}/Chr${CHRS[0]}_2k.filtered.fa" -n 40 -o "${base_path}/index_out"
  for I in $(seq 1 $(( ${#CHRS[@]} - 1)))
  do
    rimerge.py -a "${base_path}/index_out" -b "${base_path}/Chr${CHRS[$I]}/Chr${CHRS[$I]}_2k.filtered.fa" -n 40 -o "${base_path}/index_out"
  done
fi


##----------------------------------------------------------
# Print Some statistics, again
pwd; hostname; date