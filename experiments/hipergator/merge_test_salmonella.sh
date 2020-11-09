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
#SBATCH --constraint='(c27|c21)&hpg2'

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Modules
module load gcc/8.2.0
module load python3

##----------------------------------------------------------
# Vars
RIMERGE="rimerge"
BIGBWT="bigbwt"
CFA="cfa.x"
PROFILER="/usr/bin/time --verbose"
ONLIRRINDEX="OnlineRindex"

BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/experiments/merge/data"
PATH="/home/marco.oliva/bin:/home/marco.oliva/clion/r-merge/bin":$PATH
BWT_THREADS=32
MERGE_JOBS=16
EXEC_BWT=true
EXEC_MERGE=true
EXEC_ONLINE=true
SIZES=( 0 1 2 4 8 16 32 64 128 256 512 1024 2048 4096)
BASE_SIZE=5000

#BASE_PATH="/blue/boucher/marco.oliva/data/r-index/Salmonella/experiments/merge/data/test_run"
#PATH="/home/marco.oliva/bin:/home/marco.oliva/clion/r-merge/bin":$PATH
#BWT_THREADS=1
#MERGE_JOBS=16
#EXEC_BWT=true
#EXEC_MERGE=true
#EXEC_ONLINE=true
#SIZES=( 0 1 2 4 8 16 32 )
#BASE_SIZE=100

##----------------------------------------------------------
# Create index of whole files
if [ "${EXEC_BWT}" = true ]
then
  for CURR_SIZE in "${SIZES[@]}"
  do
    echo "================================================================================"
    IN_FILE="${BASE_PATH}/Salmonella_$(( ${BASE_SIZE} + ${CURR_SIZE} )).fa"
    SEQS_FILE="${BASE_PATH}/Salmonella_$(( ${BASE_SIZE} + ${CURR_SIZE} )).seqs"
    ${PROFILER} ${CFA} -i ${IN_FILE} -o ${SEQS_FILE}
    echo "Read file metadata"
    SEQS=$(od -An -t u8 "${SEQS_FILE}.meta" | awk '{print $1}')
    ${PROFILER} ${BIGBWT} -N ${SEQS} -t ${BWT_THREADS} -s -e ${SEQS_FILE}
  done
fi

##----------------------------------------------------------
# Create index by merging
if [ "${EXEC_MERGE}" = true ]
then
  echo "================================================================================"
  echo "Creating Base File (FASTA)"
  BASE_FILE="${BASE_PATH}/Salmonella_${BASE_SIZE}.fa"
  ${RIMERGE} -a ${BASE_FILE} -o "${BASE_PATH}/${BASE_SIZE}"
  for CURR_SIZE in "${SIZES[@]:1}"
  do
    echo "================================================================================"
    IN_FILE_A="${BASE_PATH}/${BASE_SIZE}"
    IN_FILE_B="${BASE_PATH}/Salmonella_${CURR_SIZE}_from_${BASE_SIZE}.fa"
    OUT_FILE="${BASE_PATH}/merged_$(( ${CURR_SIZE} + ${BASE_SIZE} ))"
    echo "Create index by merging"
    ${RIMERGE} -a ${IN_FILE_A} -b ${IN_FILE_B} -j ${MERGE_JOBS} -o ${OUT_FILE}
  done
fi

#-----------------------------------------------------------
# Online r-index
echo "================================================================================"
echo "Online r-index"
echo ""
if [ "${EXEC_ONLINE}" = true ]
then
  MAX_SIZE=${SIZES[-1]}
  SEQS_FILE="${BASE_PATH}/Salmonella_$(( ${BASE_SIZE} + ${MAX_SIZE} )).seqs"
  echo "On: ${SEQS_FILE}"
  ${PROFILER} ${ONLIRRINDEX} -i ${SEQS_FILE}
fi

##----------------------------------------------------------
# Print Some statistics, again
pwd; hostname; date