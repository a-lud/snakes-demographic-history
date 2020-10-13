#!/usr/bin/env bash
#SBATCH --job-name=filter_whole_genome
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load BCFtools/1.9-foss-2016b
module load VCFtools/0.1.14-GCC-5.3.0-binutils-2.25-Perl-5.22.0

## Directories + data
BASE="${FASTDIR}/2005_psmc"
REF="${FASTDIR}/databases/ncbi_genomes/GCA_009733165.1_Nana_v5_genomic-subset-Z-removed.fna"
VCF="${BASE}/tree/02_merged_filtered_wholeGenome/merged.tree.vcf.gz"
OUT="${BASE}/tree/03_merged_RAXML_ready"

mkdir -p ${OUT}

## Steps below
# 1. Convert GVCF to VCF
# 2. Convert GT of samples with low depth to ./.
# 3. Mark low quality sites (QUAL<20) as 'lowQual'
# 4. Keep only SNPs and remove any sites that are not 'PASS'
# 5. Keep only sites where at most 1 sample is missing data

bcftools convert --gvcf2vcf -f ${REF} ${VCF} | \
bcftools filter -S . -e 'FORMAT/DP<10' | \
bcftools filter -s lowQual -e 'QUAL<20' | \
bcftools view -v snps -f PASS -O v | \
vcftools --vcf - --max-missing-count 1 --recode --recode-INFO-all -c | \
bgzip > ${OUT}/merged.wholeGenome.RAXML_ready.vcf.gz

