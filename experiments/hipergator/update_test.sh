#!/bin/bash
#SBATCH --job-name=performance          # Job name
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marco.oliva@ufl.edu # Where to send mail
#SBATCH --ntasks=1                   	  # Run a single task
#SBATCH --cpus-per-task=2           	  # Number of CPU cores per task
#SBATCH --mem=20gb                     	# Job memory request
#SBATCH --time=168:00:00              	# Time limit hrs:min:sec
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
base_path="/blue/boucher/marco.oliva/data/r-index/1000GP/experiments/merge_performance"
PATH="/home/marco.oliva/bin":$PATH
BWT_THREADS=32
MERGE_JOBS=16
SKIP=1000
SEQS_BWT=( 1000 1001 1002 1004 1008 1016 1032 1064 1128 1256 1512 )
SEQS_MRG=( 1 2 4 8 16 32 64 128 256 512 )
ITERATIONS=1

#-----------------------------------------------------------
# Create files and test bigbwt
echo "================================================================================"
echo "Create files and test bigbwt"
echo ""
for S in "${SEQS_BWT[@]}"
do
  cfa -n ${S} -i "${base_path}/Chr_19_2000_bwt.filtered.fa" -o "${base_path}/Chr_19_2000_bwt.filtered.${S}.seqs"
  for i in $(seq 1 $ITERATIONS)
  do
    echo "========================================"
    echo "Iteration ${i}"
    ${profiler} bigbwt -N $(($S + 1)) -t ${BWT_THREADS} -s -e "${base_path}/Chr_19_2000_bwt.filtered.${S}.seqs"
  done
  estw -i "${base_path}/Chr_19_2000_bwt.filtered.${S}.seqs"
  rle  -i "${base_path}/Chr_19_2000_bwt.filtered.${S}.seqs"
done

#-----------------------------------------------------------
# Create files to be merged and test bigbwt
echo "================================================================================"
echo "Create files to be merged and test bigbwt"
echo ""
for S in "${SEQS_MRG[@]}"
do
  cfa -n ${S} -I ${SKIP} -i "${base_path}/Chr_19_2000_mrg.filtered.fa" -o "${base_path}/Chr_19_2000_mrg.filtered.${S}.seqs"
  for i in $(seq 1 $ITERATIONS)
  do
    echo "========================================"
    echo "Iteration ${i}"
    ${profiler} bigbwt -N $(($S + 1)) -t ${BWT_THREADS} -s -e "${base_path}/Chr_19_2000_mrg.filtered.${S}.seqs"
  done
  estw -i "${base_path}/Chr_19_2000_mrg.filtered.${S}.seqs"
  rle  -i "${base_path}/Chr_19_2000_mrg.filtered.${S}.seqs"
done

#-----------------------------------------------------------
# Merge
echo "================================================================================"
echo "Merge"
echo ""
for S in "${SEQS_MRG[@]}"
do
  for i in $(seq 1 $ITERATIONS)
  do
    echo "========================================"
    echo "Iteration ${i}"
    ${profiler} main -a "${base_path}/Chr_19_2000_bwt.filtered.${SKIP}.seqs" -b "${base_path}/Chr_19_2000_mrg.filtered.${S}.seqs" -o "${base_path}/mrg_out_${S}" -j ${MERGE_JOBS}
  done
done

#-----------------------------------------------------------
# Online r-index
echo "================================================================================"
echo "Online r-index"
echo ""
#cfa -i "${base_path}/Chr_19_2000_bwt.filtered.fa" -o "${base_path}/Chr_19_2000_bwt.filtered.seqs"
for i in $(seq 1 $ITERATIONS)
do
  ${profiler} OnlineRindex -i "${base_path}/Chr_19_2000_bwt.filtered.seqs"
done

##----------------------------------------------------------
# Print Some statistics, again
pwd; hostname; date