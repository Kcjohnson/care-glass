#!/bin/bash

### This script runs GATK/Picard's liftover on VCF files in ssm2filter results ###

# Activate the conda environment that contains GATK4.
conda activate /projects/verhaak-lab/USERS/johnsk/glass4/.snakemake/conda/e8563520

### Input arguments ###
ARRAYID="`expr $1`"

# List of barcodes for which we have VCF files
VCF_FILE_LIST="/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/ssm2filter_liftover/pair_barcodes_to_liftover.txt"

# Target File.
VCF_FILE=$(sed "${ARRAYID}q;d" ${VCF_FILE_LIST})

# Extract the pair barcoe.
PAIR_BARCODE=$(basename ${VCF_FILE} ".filtered.vcf.gz")

### Print out relevant information.
STARTTIME=`date`
echo $STARTTIME
echo "Analyzing $PAIR_BARCODE"
echo ""

gatk --java-options -Xmx12g LiftoverVcf \
    --I=/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/ssm2filter/${PAIR_BARCODE}.filtered.vcf.gz \
    --O=/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/ssm2filter_liftover/${PAIR_BARCODE}.filtered.hg38.vcf \
    --CHAIN=/projects/verhaak-lab/USERS/johnsk/glass4/bin/b37ToHg38.over.chain \
    --REJECT=/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/ssm2filter_liftover/${PAIR_BARCODE}.rejected_variants.vcf \
    --R=/projects/verhaak-lab/USERS/johnsk/glass4/ref/hg38.fa


### END ####
