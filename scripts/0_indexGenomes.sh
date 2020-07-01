#!/usr/bin/env bash
#SBATCH --job-name=indexing
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/slurm/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

module load BWA/0.7.17-foss-2016b

REFS='/home/a1645424/fastdir/databases/ncbi_genomes'

for f in ${REFS}/*.fna; do

    bwa index $f

done