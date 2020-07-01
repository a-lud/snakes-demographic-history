#!/usr/bin/env bash

for r in aipysurusLaevis GCA_004320025.1_latLat GCA_004320045.1_latCor GCA_004320005.1_hydMel GCF_900518735.1_EBS GCF_900518725.1_TS GCA_009733165.1_Nana; do
    for s in DRR144984_DRR144985 DRR147552 DRR147394 ERR2714264_ts ERR2714265 KLS0691 SRR10428161; do

        ## Submission to SLURM
        sbatch \
        --job-name="bwa_${r}_${s}" \
        --export=REFDIR="/home/a1645424/fastdir/databases/ncbi_genomes",REF=${r},SAMPDIR="/data/biohub/180704_kateSanders/00_sequenceData/wgs/ncbi-samples",SAMP=${s},PROJDIR="/fast/users/a1645424/2005_psmc" \
        /fast/users/a1645424/2005_psmc/scripts/1_bwa-template.sh

    done
done
