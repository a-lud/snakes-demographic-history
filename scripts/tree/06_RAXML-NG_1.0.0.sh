#!/bin/bash
#SBATCH --job-name=RAXML_test
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/tree/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load RAxML-NG/1.0.0-foss-2016b

BASE="/home/a1645424/fastdir/2005_psmc"
PHY="${BASE}/tree/03_merged_RAXML_ready"
OUT="${BASE}/tree/04_RAXML"

raxml-ng \
--threads ${SLURM_CPUS_PER_TASK} \
--all \
--seed 12345 \
--msa ${PHY}/merged.wholeGenome.RAXML_ready.min7.invariantRemoved.phy \
--msa-format PHYLIP \
--prefix elapid_raxmlNG \
--outgroup SRR10428161 \
--model GTGTR4+G+ASC_LEWIS \
--bs-trees 1000
