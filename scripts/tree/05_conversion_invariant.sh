#!/usr/bin/env bash
#SBATCH --job-name=conversion_invariant
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:10:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

## Directories + data
BASE="${FASTDIR}/2005_psmc"
VCF="${BASE}/tree/03_merged_RAXML_ready"

vcf2phylip.py \
-i ${VCF}/merged.wholeGenome.RAXML_ready.vcf.gz \
-m 7 \
-o SRR10428161 \
-f

./ascbias.py \
-p ${VCF}/merged.wholeGenome.RAXML_ready.min6.phy \
-o ${VCF}/merged.wholeGenome.RAXML_ready.min6.invariantRemoved.phy 
