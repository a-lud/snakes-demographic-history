#!/usr/bin/env bash

for r in "aipysurusLaevis 1852689089" "GCA_004320045.1_latCor 2024687924" "GCF_900518735.1_EBS 1590035073" "GCA_004320005.1_hydMel 1402639853" "GCA_009733165.1_Nana 1667692826" "GCA_004320025.1_latLat 1558706106" "GCF_900518725.1_TS 1665525958"; do

    set -- ${r}

    for s in DRR144984_DRR144985 DRR147552 DRR147394 ERR2714264_ts ERR2714265 KLS0691 SRR10428161; do

        ## Submission to SLURM
        sbatch \
        --job-name="depth_${1}_${s}" \
        --export=REF=${1},SIZE=${2},SAMP=${s},PROJDIR="/fast/users/a1645424/2005_psmc" \
        /fast/users/a1645424/2005_psmc/scripts/2_depth-template.sh

    done
done
