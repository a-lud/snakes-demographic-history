#!/usr/bin/env bash
#SBATCH --job-name=bcftools_merge_filter
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=72:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load BCFtools/1.9-foss-2016b

## Directories + data
BASE="${FASTDIR}/2005_psmc"
REF="${FASTDIR}/databases/ncbi_genomes/GCA_009733165.1_Nana_v5_genomic-subset-Z-removed.fna"
BAMS="${BASE}/tree/01_bcftools_wholeGenome"
OUT="${BASE}/tree/02_merged_filtered_wholeGenome"

mkdir -p ${OUT}

bcftools merge --gvcf ${REF} ${BAMS}/DRR144984_DRR144985.vcf.gz ${BAMS}/ERR2714265.vcf.gz ${BAMS}/SRR10861675_SRR10861676.vcf.gz ${BAMS}/DRR147394.vcf.gz ${BAMS}/KLS0691.vcf.gz ${BAMS}/ERR2714264_ts.vcf.gz ${BAMS}/SRR10428161.vcf.gz | \
bcftools view -m 2 -M 2 -O z -o ${OUT}/merged.tree.vcf.gz

tabix -p vcf ${OUT}/merged.tree.vcf.gz

