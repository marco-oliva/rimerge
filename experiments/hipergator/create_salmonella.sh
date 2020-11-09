#!/bin/bash
#SBATCH --job-name=mrg_sal_d            # Job name
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=marco.oliva@ufl.edu # Where to send mail
#SBATCH --ntasks=1                   	# Run a single task
#SBATCH --cpus-per-task=2           	# Number of CPU cores per task
#SBATCH --mem=16gb                     	# Job memory request
#SBATCH --time=168:00:00              	# Time limit hrs:min:sec
#SBATCH --output=%j_mrg_sal_d.log       # Standard output and error log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Modules
module load gcc/8.2.0
module load python3

##----------------------------------------------------------
# Vars
#SIZES=( 1 2 4 8 16 32 64 128 256 512 1024 2048 4096)
#BASE_SIZE=5000
#BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/kuhnle/genomes"
#OUT_BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/experiments/merge/data"
#PATH="/home/marco.oliva/bin":$PATH

SIZES=( 1 2 4 8 16 32 )
BASE_SIZE=100
BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/kuhnle/genomes"
OUT_BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/experiments/merge/data/test_run"
PATH="/home/marco.oliva/bin":$PATH

#-----------------------------------------------------------
# Create dataset

set -o xtrace

mkdir -p ${OUT_BASE_PATH}

# base file
cd "${BASE_PATH}"
find . -not -empty | head -n $(( ${BASE_SIZE} + 2 )) | xargs -I '{}' cat '{}' >> ${OUT_BASE_PATH}/Salmonella_${BASE_SIZE}.fa

# incremental files
for CURR_SIZE in "${SIZES[@]}"
do
    echo "Generating file ${CURR_SIZE}"
    cd "${BASE_PATH}"
    find . -not -empty | head -n $(( ${BASE_SIZE} + ${CURR_SIZE} + 1 )) | tail -n "${CURR_SIZE}" | xargs -I '{}' cat '{}' >> "${OUT_BASE_PATH}/Salmonella_${CURR_SIZE}_from_${BASE_SIZE}.fa"
    echo "Donefile ${CURR_SIZE}"
done

# full size files
for CURR_SIZE in "${SIZES[@]}"
do
    OUT_SIZE=$(( ${BASE_SIZE} + ${CURR_SIZE} ))
    echo "Generating file ${OUT_SIZE}"
    cd "${BASE_PATH}"
    find . -not -empty | head -n $(( ${OUT_SIZE} + 1 )) | xargs -I '{}' cat '{}' >> "${OUT_BASE_PATH}/Salmonella_${OUT_SIZE}.fa"
    echo "Donefile ${CURR_SIZE}"
done


set +o xtrace

