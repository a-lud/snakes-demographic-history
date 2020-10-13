#!/usr/bin/env bash
#SBATCH --job-name=bcftools_variants
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -a 1-7
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load BCFtools/1.9-foss-2016b

## Directories + data
BASE="${FASTDIR}/2005_psmc"
REF="${FASTDIR}/databases/ncbi_genomes/GCA_009733165.1_Nana_v5_genomic-subset-Z-removed.fna"
BAMS="${BASE}/tree/bams_markdup"
OUT="${BASE}/tree/01_bcftools_wholeGenome"

mkdir -p ${OUT}

BAM=$(find ${BAMS} -type f -name "*.markdup.bam" | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
BN=$(basename ${BAM} .markdup.bam)

mkdir -p ${OUT}/temp_${BN}

bcftools mpileup -O v --gvcf 5,15 -a 'AD,DP' -C 50 -q 20 -Q 20 -f ${REF} ${BAM} | \
bcftools call -c | \
bcftools norm -m- | \
bcftools sort --temp-dir ${OUT}/temp_${BN} -Oz -o ${OUT}/${BN}.vcf.gz

# bcftools convert --gvcf2vcf -f ${REF} | \

tabix -p vcf ${OUT}/${BN}.vcf.gz

