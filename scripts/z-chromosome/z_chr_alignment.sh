#!/usr/bin/env bash
#SBATCH --job-name=z_chr_alignment
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH -o /home/a1645424/fastdir/2005_psmc/slurm/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

ZCHR="${FASTDIR}/databases/ncbi_genomes/GCA_009733165.1_Nana_v5_z-chromosome.fasta"
GEN_DIR="${FASTDIR}/databases/ncbi_genomes"
WD='/home/a1645424/fastdir/2005_psmc'

source activate MUMMER_env

mkdir -p ${WD}/z_chr_alignments

for f in ${GEN_DIR}/GCA_009733165.1_Nana_v5_genomic-subset.fna /home/a1645424/fastdir/2005_psmc/data/genome/Hcur1.v1.1.fna; do
    nucmer -t ${SLURM_CPUS_PER_TASK} -p ${WD}/z_chr_alignments/mummer4_alignments/$(basename ${f} .fna) ${ZCHR} ${f}
    show-coords -r -c -l -T -I 80 ${WD}/z_chr_alignments/mummer4_alignments/$(basename ${f} .fna).delta > ${WD}/z_chr_alignments/mummer4_alignments/$(basename ${f} .fna).coords
done

conda deactivate
