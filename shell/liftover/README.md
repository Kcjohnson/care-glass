# Running LiftoverVcf for the CARE samples

Background: The CARE project is integrating longitudinal bulk DNAseq and single nucelus RNAseq data on the same samples. As part of the snRNAseq analysis, our collaborator Avishay Spitzer is calling CNVs in the snRNA data, which is a noisy process in some samples. To provide an extra layer of validation on these samples Avishay is running a workflow to call mutations in the 10X data. To do so Avishay needs the bulk DNAseq VCF files (Mutect2 single sample mode) to be reported in hg38 coordinates.

Since the GLASS pipeline is built around hg19 and there is no immediate plan to run hg38, we decided to run "LiftoverVcf" tool from Picard/GATK. The following details what I did to run the workflow.

## Step 1: Download the appropriate reference genome

Avishay initially wanted me to liftover to his specific reference genome that he had sent to K2 to align our snRNAseq data (GRCh38-3.0.0_premrna, located: `/projects/verhaak-lab/scgp2/reference/genomes/hsapiens/Suva_group/GRCh38-3.0.0_premrna`). However, trying to lift over to this version did not immediately work because some differences in the annotation of the genome. Instead, I simply went with the primary hg38 assembly in a few test samples, which Avishay reported as working.

```sh
# Using an rsync command to download the entire directory from UCSC:
mkdir -p /fastscratch/johnsk/hg38
cd /fastscratch/johnsk/hg38
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/ .

# unzip the reference genome and test files of interest.
gunzip hg38.fa.gz
```

I encountered an error when initially running the LiftoverVcf tool. It reported that the reference file must have an associated Dictionary .dict file in the same directory.

```sh
# Activate the conda environment that contains GATK4.
conda activate /projects/verhaak-lab/USERS/johnsk/glass4/.snakemake/conda/e8563520

gatk --java-options -Xmx12g CreateSequenceDictionary \
      --R=/fastscratch/johnsk/hg38/hg38.fa \
      --O=/fastscratch/johnsk/hg38/hg38.dict
```

Since the GLASS data was aligned to b37, I needed a specific chain file to be able to map the coordinates. I could not find the exact file I was looking for in the GATK resource bundle. However, I was able to find it on the Broad's github.

## Step 2: Retrieve the appropriate chain file.


```sh
wget https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/data_sources/gnomAD/b37ToHg38.over.chain .
```

Here's the template command for the liftover:
```sh
# Here is the template command.
 java -jar picard.jar LiftoverVcf \
     I=input.vcf \
     O=lifted_over.vcf \
     CHAIN=b37tohg38.chain \
     REJECT=rejected_variants.vcf \
     R=reference_sequence.fasta
```

I then prepared a test sample for Avishay to see whether this approach would work

```sh
# Prepare the two test samples for Avishay
gatk --java-options -Xmx12g LiftoverVcf \
    --I=/fastscratch/johnsk/vcf/CARE-MD-MD01-R1-01-NB-01D-WXS.filtered.vcf \
    --O=/fastscratch/johnsk/vcf/CARE-MD-MD01-R1-01-NB-01D-WXS.filtered.hg38.vcf \
    --CHAIN=/fastscratch/johnsk/hg38/b37ToHg38.over.chain \
    --REJECT=/fastscratch/johnsk/vcf/CARE-MD-MD01-R1-01-NB-01D-WXS.rejected_variants.vcf \
    --R=/fastscratch/johnsk/hg38/hg38.fa

grep -E "PASS" CARE-MD-MD01-R1-01-NB-01D-WXS.filtered.vcf | wc -l
grep -E "PASS" CARE-MD-MD01-R1-01-NB-01D-WXS.filtered.hg38.vcf | wc -l
```

It appeared that most of the variants lifted over and Avishay had success with implementing his method.


