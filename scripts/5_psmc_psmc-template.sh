#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=72:00:00
#SBATCH --mem=16GB
#SBATCH -o /fast/users/a1645424/2005_psmc/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

export PATH="${PATH}:/home/a1645424/fastdir/tools_fd/psmc"
export PATH="${PATH}:/home/a1645424/fastdir/tools_fd/psmc/utils"

## Variables
PSMCFA=${PROJDIR}/02_consensus/psmcfa/${REF}
CLOCK_CLEAN=$(echo ${CLOCK} | sed 's/\*/x/g;s/\+/_/g')

## Directories
BASE=${PROJDIR}/03_analysis/psmc/${REF}
OUT=${BASE}/${SAMP}
OUTBOOT=${OUT}/bootstrap_${SAMP}_${CLOCK_CLEAN}

mkdir -p ${OUTBOOT}

## Check if file exists
if [[ -f ${BASE}/${SAMP}-${CLOCK_CLEAN}-combined.psmc ]]; then
    echo "File already exists: ${BASE}/${SAMP}-${CLOCK_CLEAN}-combined.psmc"
    exit 0
fi

## Check that the psmcfa file is not 0b
if [[ ! -s ${PSMCFA}/${SAMP}.psmcfa ]]; then
    echo "Empty psmcfa file: ${PSMCFA}/${SAMP}.psmcfa"
    exit 1
fi

## PSMC
echo "PSMC diploid"
psmc -p ${CLOCK} -o ${OUT}/${SAMP}-${CLOCK_CLEAN}.psmc ${PSMCFA}/${SAMP}.psmcfa

echo -e "PSMC bootstrap"
splitfa ${PSMCFA}/${SAMP}.psmcfa > ${PSMCFA}/${SAMP}-split.psmcfa

seq 100 | \
xargs \
-n 1 \
-P ${SLURM_CPUS_PER_TASK} \
-I {} psmc -p ${CLOCK} -b -o ${OUTBOOT}/${SAMP}-${CLOCK_CLEAN}-{}.psmc ${PSMCFA}/${SAMP}-split.psmcfa | \
sh

## Combine diploid and bootstrap output to single file
cat \
${OUT}/${SAMP}-${CLOCK_CLEAN}.psmc ${OUTBOOT}/${SAMP}-${CLOCK_CLEAN}-*.psmc > \
${BASE}/${SAMP}-${CLOCK_CLEAN}-combined.psmc

