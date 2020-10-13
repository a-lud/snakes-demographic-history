#!/bin/bash
#SBATCH --job-name=RAXML
#SBATCH -p test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=2:00:00
#SBATCH --mem=8GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load RAxML/8.2.10-foss-2016b-pthreads-avx

BASE="/home/a1645424/fastdir/2005_psmc"
PHY="${BASE}/tree/03_merged_RAXML_ready"
OUT="${BASE}/tree/04_RAXML"

raxmlHPC-PTHREADS-AVX \
-T ${SLURM_CPUS_PER_TASK} \
-f a \
-x 12345 \
-p 56422 \
-# 1000 \
-m ASC_GTRGAMMA \
-n elapidSNP \
-o SRR10428161 \
--asc-corr=lewis \
-s ${PHY}/merged.wholeGenome.RAXML_ready.min7.invariantRemoved.phy
