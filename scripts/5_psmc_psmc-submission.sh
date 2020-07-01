#!/usr/bin/env bash

for r in aipysurusLaevis GCA_004320025.1_latLat GCA_004320045.1_latCor GCA_009733165.1_Nana GCF_900518735.1_EBS GCA_004320005.1_hydMel GCF_900518725.1_TS; do
    for s in DRR144984_DRR144985 DRR147552 DRR147394 ERR2714264_ts ERR2714265 KLS0691 SRR10428161; do
        for c in 4+5*3+4 4+10*3+6+8 4+25*2+4+6 4+30*2+4+6+10 ; do

            ## Submission to SLURM
            sbatch \
            --job-name="psmc_${r}_${s}" \
            --export=REF=${r},SAMP=${s},CLOCK=${c},PROJDIR="/fast/users/a1645424/2005_psmc" \
            /fast/users/a1645424/2005_psmc/scripts/psmc_7_psmc-template.sh

        done
    done
done
