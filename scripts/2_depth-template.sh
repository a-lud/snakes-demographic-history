#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=2GB
#SBATCH -o /fast/users/a1645424/2005_psmc/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

## Loading modules
module load SAMtools/1.9-foss-2016b

## Paths
BAMDIR=${PROJDIR}/01_bwa/${REF}

## Data
BAM=${BAMDIR}/${SAMP}.sorted.bam

## Check if file exists
if [[ -f ${BAMDIR}/${SAMP}.avgDepth ]]; then
    echo "Coverage file exists: ${BAMDIR}/${SAMP}.avgDepth"
    exit 0
fi

## Coverage calculations
AVG=$(samtools depth -a ${BAM} | \
awk '{sum+=$3} END { print sum/'${SIZE}' }')

BTH=$(samtools depth -a ${BAM} | \
awk '{if($3>0) total+=1} END { print (total/'${SIZE}') * 100 }')

## Mapping statistics
samtools flagstat ${BAM} > ${BAMDIR}/${SAMP}.flagstat

## To file
echo ${AVG} > ${BAMDIR}/${SAMP}.avgDepth

echo -e "reference\tsample\tgenome_size\tavg_depth\tbreadth_depth
${REF}\t${SAMP}\t${SIZE}\t${AVG}\t${BTH}" > ${BAMDIR}/${SAMP}.coverageStats.tsv

