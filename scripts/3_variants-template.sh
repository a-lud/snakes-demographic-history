#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH -o /fast/users/a1645424/2005_psmc/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

## Modules
module load BCFtools/1.9-foss-2016b

## Files
BAM=${PROJDIR}/01_bwa/${REF}/${SAMP}.sorted.bam
GENOME=$(find ${REFDIR} -mindepth 1 -maxdepth 1 -type f -name "${REF}*.fna")
DEPTHAVG=$(cat ${PROJDIR}/01_bwa/${REF}/${SAMP}.avgDepth)

DEPTHMAX=$(echo "${DEPTHAVG}*2" | bc)
DEPTHMIN=$(echo "${DEPTHAVG}/3" | bc)

## Directories
OUT=${PROJDIR}/01_bwa/vcf/${REF}

mkdir -p ${OUT}

echo -e "REF: ${REF}\nSAMP: ${SAMP}\nAVG: ${DEPTHAVG}\nMAX: ${DEPTHMAX}\nMIN: ${DEPTHMIN}"

## Check output file exists
if [[ -f ${OUT}/${SAMP}.vcf.gz ]]; then
    echo "File already exists: ${OUT}/${SAMP}.vcf.gz"
    exit 0
fi

## Creating temporary working directory
mkdir -p ${OUT}/temp_${SAMP}

## Calling variants
bcftools mpileup -C 50 -q 20 -Q 25 -Ou -f ${GENOME} ${BAM} | \
bcftools call -c | \
bcftools filter --SnpGap 10 -i "DP>=${DEPTHMIN} & DP<=${DEPTHMAX}" | \
bcftools view --exclude-types indels | \
bcftools sort --temp-dir ${OUT}/temp_${SAMP} -Oz -o ${OUT}/${SAMP}.vcf.gz

rm -R ${OUT}/temp_${SAMP}

tabix -p ${OUT}/${SAMP}.vcf.gz

## VCF statistics
bcftools stats ${OUT}/${SAMP}.vcf.gz > ${OUT}/${SAMP}.vcf.stats

if [[ -s ${OUT}/${SAMP}.vcf.stats ]]; then
    rm ${BAM}
fi
