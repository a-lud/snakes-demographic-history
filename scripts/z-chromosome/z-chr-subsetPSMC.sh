#!/usr/bin/env bash

module load seqtk/1.3-foss-2016b

PROJ="/home/a1645424/fastdir/2005_psmc"
PSMCIDS="${PROJ}/z_chr_alignments/psmc_filter_ids"
PSMCFADIR="${PROJ}/02_consensus/psmcfa"
OUTDIR="${PROJ}/02_consensus/psmcfa_50"
REFS="aipysurusLaevis  GCA_004320005.1_hydMel  GCA_004320025.1_latLat  GCA_009733165.1_Nana  GCF_900518725.1_TS  GCF_900518735.1_EBS  Hcur1.v1.1"

for r in ${REFS}; do

    mkdir -pv ${OUTDIR}/${r}

    case ${r} in
        aipysurusLaevis)
            IDS="${PSMCIDS}/aipysurusLaevis_50.txt"
            ;;
        GCA_009733165.1_Nana)
            IDS="${PSMCIDS}/GCA_009733165.1_Nana_v5_genomic-subset_50.txt"
            ;;
        Hcur1.v1.1)
            IDS="${PSMCIDS}/Hcur1.v1.1_50.txt"
            ;;
        GCA_004320005.1_hydMel)
            IDS="${PSMCIDS}/GCA_004320005.1_hydMel_1.0_genomic_50.txt"
            ;;
        GCF_900518725.1_TS)
            IDS="${PSMCIDS}/GCF_900518725.1_TS10Xv2-PRI_genomic_50.txt"
            ;;
        GCA_004320025.1_latLat)
            IDS="${PSMCIDS}/GCA_004320025.1_latLat_1.0_genomic_50.txt"
            ;;
        GCF_900518735.1_EBS)
            IDS="${PSMCIDS}/GCF_900518735.1_EBS10Xv2-PRI_genomic_50.txt"
            ;;
    esac

    for f in ${PSMCFADIR}/${r}/*.psmcfa; do

        BN=$(basename ${f})

        if [[ -s ${OUTDIR}/${r}/${BN} ]]; then
            echo -e "File already exists: ${OUTDIR}/${r}/${BN}"
            :
        elif [[ ${f} == *"split"* ]]; then
            echo -e "File is a split file: ${f}\n"
            :
        else
            echo -e "Reference: ${r}\nFile: ${f}"
            seqtk subseq ${f} ${IDS} > ${OUTDIR}/${r}/${BN}
        fi

    done

done
