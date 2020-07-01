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

## Loading modules
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.9-foss-2016b

## Paths
OUTDIR=${PROJDIR}/01_bwa/${REF}
mkdir -p ${OUTDIR}

## Data
R1=$(find ${SAMPDIR} -type f -name "${SAMP}*_1*")
R2=$(find ${SAMPDIR} -type f -name "${SAMP}*_2*")
GENOME=$(find ${REFDIR} -type f -name "${REF}*.fna")

printf "%0.s#" {1..80}
echo -e "Read 1: ${R1}\nRead 2: ${R2}\nGenome: ${GENOME}\n"

## Exit if file exists already
if [[ -f ${OUTDIR}/${SAMP}.sorted.bam ]]; then
    echo -e "BAM file already exists: ${OUTDIR}/${SAMP}.sorted.bam"
    exit 0
fi

## Align the data
bwa mem \
-t ${SLURM_CPUS_PER_TASK} \
-R "@RG\tID:${SAMP}\tLB:${SAMP}_WGS\tPL:ILLUMINA\tSM:${SAMP}" \
-M \
${GENOME} \
${R1} \
${R2} | \
samtools sort | \
samtools view -O BAM -o ${OUTDIR}/${SAMP}.sorted.bam

## Generate index
samtools index -@ ${SLURM_CPUS_PER_TASK} ${OUTDIR}/${SAMP}.sorted.bam

