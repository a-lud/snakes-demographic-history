#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -a 1-7
#SBATCH --time=72:00:00
#SBATCH --mem=16GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

## Loading modules
module load BWA/0.7.15-foss-2017a
module load SAMtools/1.9-foss-2016b

## Directories + data
BASE="${FASTDIR}/2005_psmc"
READS="${FASTDIR}/2005_psmc/data/reads"
REF="${FASTDIR}/databases/ncbi_genomes/GCA_009733165.1_Nana_v5_genomic-subset-Z-removed.fna"
OUT="${BASE}/tree/bams"

## Current sample information
R1=$(find ${READS} -type f -name "*_1.fastq.gz" | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
R2=${R1/_1.fastq.gz/_2.fastq.gz}
SAMP=$(basename ${R1} _1.fastq.gz)

echo -e "Sample: ${SAMP}\nRead 1: ${R1}\nRead 2: ${R2}\n"

## Alignment
bwa mem \
-t ${SLURM_CPUS_PER_TASK} \
-R "@RG\tID:${SAMP}\tLB:${SAMP}_WGS\tPL:ILLUMINA\tSM:${SAMP}" \
-M \
${REF} \
${R1} \
${R2} | \
samtools sort | \
samtools view -O BAM -o ${OUT}/${SAMP}.sorted.bam

## Generate index
samtools index ${OUT}/${SAMP}.sorted.bam
