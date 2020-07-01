#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH -o /fast/users/a1645424/2005_psmc/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

## Modules
module load BCFtools/1.9-foss-2016b
PSMC="/home/a1645424/fastdir/tools_fd/psmc"

## Files
VCF=${PROJDIR}/01_bwa/vcf/${REF}/${SAMP}.vcf.gz

## Directories
OUT=${PROJDIR}/02_consensus/psmcfa/${REF}

mkdir -p ${OUT}

## Check output file exists
if [[ -f ${OUT}/${SAMP}.psmcfa ]]; then
    echo "File already exists: ${OUT}/${SAMP}.psmcfa"
    exit 0
fi

## Build consensus
bcftools view ${VCF} | \
vcfutils.pl vcf2fq | \
gzip > ${OUT}/${SAMP}.fastq.gz

## psmcfa
${PSMC}/utils/fq2psmcfa \
-q 20 \
${OUT}/${SAMP}.fastq.gz > ${OUT}/${SAMP}.psmcfa

## Remove fastq file
rm ${OUT}/${SAMP}.fastq.gz

## Remove empty psmcfa
if [[ ! -s ${OUT}/${SAMP}.psmcfa ]];then 
    rm ${OUT}/${SAMP}.psmcfa;
fi

## Remove VCF if output file is not empty
if [[ -s ${OUT}/${SAMP}.psmcfa ]]; then
    rm ${VCF}
fi
