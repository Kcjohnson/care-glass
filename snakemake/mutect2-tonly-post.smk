## Changes 2022.07.10
## - Adding tumor-only mutation calling for CARE project.
## - To enable comparisons with matched normal samples, we will place all results for each cohort in results/mutect2/tonly/ followed by the same structure as the matched normal blood versions.
## - Introduced a `getSelectedAliquotsTonly()` and getSelectedCasesTonly() functions in the ManifestHandler.py to only grab tumor-only cases and prevent data from being uploaded to the database twice. 
## - Note that these samples were hardcoded in the ManifestHandler.py.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Select variants
## remove selectvariants and only drop GTs
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule selectvariants_tonly:
    input:
        ancient("results/mutect2/tonly/m2filter/{case_barcode}.filtered.vcf.gz")
    output:
        normalized = temp("results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.vcf.gz"),
        sortd = "results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz",
        sortdind = "results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz.tbi",
        vcf = temp("results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz"),
        tbi = temp("results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["selectvariants_tonly"]["mem"]
    threads:
        CLUSTER_META["selectvariants_tonly"]["cpus-per-task"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/tonly/selectvariants/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/selectvariants/{case_barcode}.txt"
    message:
        "Drop genotypes\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref we \
            -m-both \
            {input} | \
         bcftools view \
            -Oz \
            -o {output.normalized} \
            2>> {log};"
       
        "bcftools sort \
            -Oz \
            -o {output.sortd} \
            {output.normalized} \
            2>> {log};"
        
        "bcftools index \
            -t {output.sortd} \
            2>> {log};"

        "bcftools view -Oz -G -o {output.vcf} {output.sortd} >> {log} 2>&1;"
        "bcftools index -t {output.vcf} >> {log} 2>&1;"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single sample Select variants
## Changed directory to tonly and replaced pair_barcode with aliquot_barcode
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssselectvariants_tonly:
    input:
        expand("results/mutect2/tonly/ssm2filter/{aliquot_barcode}.filtered.vcf.gz", aliquot_barcode = manifest.getSelectedAliquotsTonly())
    output:
        normalized = temp("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.vcf.gz"),
        sortd = temp("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.sorted.vcf.gz"),
        sortdind = temp("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.sorted.vcf.gz.tbi"),
        vcf = temp("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.sorted.dropgt.vcf.gz"),
        tbi = temp("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["ssselectvariants_tonly"]["mem"]
    threads:
        CLUSTER_META["ssselectvariants_tonly"]["cpus-per-task"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/tonly/ssselectvariants/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/ssselectvariants/{aliquot_barcode}.txt"
    message:
        "Single Sample Drop genotypes\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref we \
            -m-both \
            {input} | \
         bcftools view \
            -Oz \
            -o {output.normalized} \
            2>> {log};"
       
        "bcftools sort \
            -Oz \
            -o {output.sortd} \
            {output.normalized} \
            2>> {log};"
        
        "bcftools index \
            -t {output.sortd} \
            2>> {log};"

        "bcftools view -Oz -G -o {output.vcf} {output.sortd} >> {log} 2>&1;"
        "bcftools index -t {output.vcf} >> {log} 2>&1;"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule consensusvcf_tonly:
    input:
        vcfs = expand("results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz", case_barcode = manifest.getSelectedCasesTonly()),
        tbis = expand("results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.dropgt.vcf.gz.tbi", case_barcode = manifest.getSelectedCasesTonly())
    output:
        merged = temp("results/mutect2/tonly/consensusvcf/consensus.vcf.gz"),
        normalized = temp("results/mutect2/tonly/consensusvcf/consensus.normalized.vcf.gz"),
        final = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vcf.gz",
        tbi = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi",
    params:
        mem = CLUSTER_META["consensusvcf_tonly"]["mem"]
    threads:
        CLUSTER_META["consensusvcf_tonly"]["cpus-per-task"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/mutect2/tonly/consensusvcf/consensusvcf.log"
    benchmark:
        "benchmarks/mutect2/tonly/consensusvcf/consensusvcf.txt"
    message:
        "Merge consensus variants"
    shell:
        "bcftools merge \
            -m none {input.vcfs} \
            -Oz \
            -o {output.merged} \
            2>> {log};"
        
        "bcftools norm \
            -f {config[reference_fasta]} \
            --check-ref we \
            -m-both \
            {output.merged} | \
         bcftools view \
            -Oz \
            -o {output.normalized} \
            2>> {log};"
       
        "bcftools sort \
            -Oz \
            -o {output.final} \
            {output.normalized} \
            2>> {log};"
        
        "bcftools index \
            -t {output.final} \
            2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annoconsensusvcf_tonly:
    input:
        vcf = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vcf.gz",
        tbi = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vcf.gz.tbi"
    output:
        vcf = temp("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.vcf"),
        idx = temp("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.vcf.idx"),
        gz = protected("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz"),
        tbi = protected("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["annoconsensusvcf_tonly"]["mem"]
    threads:
        CLUSTER_META["annoconsensusvcf_tonly"]["cpus-per-task"]
    conda:
        "../envs/funcotate.yaml"
    log:
        "logs/mutect2/tonly/annoconsensusvcf/annoconsensusvcf.log"
    benchmark:
        "benchmarks/mutect2/tonly/annoconsensusvcf/annoconsensusvcf.txt"
    message:
        "Annotate consensus variants using Funcotator"
    shell:
        "gatk Funcotator \
            --variant {input.vcf} \
            --reference {config[reference_fasta]} \
            --ref-version hg19 \
            --data-sources-path {config[funcotator_dir]} \
            --output {output.vcf} \
            --output-file-format VCF \
            --remove-filtered-variants false \
            2>> {log};"
        
        "bcftools view \
            -Oz \
            -o {output.gz} \
            {output.vcf} \
            2>> {log};"

        "bcftools index \
            -t {output.gz} \
            2>> {log};"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule maf2db_tonly:
    input:
        vcf = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.vcf.gz"
    output:
        tsv = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.funcotated.tsv"
    params:
        mem = CLUSTER_META["maf2db_tonly"]["mem"]
    threads:
        CLUSTER_META["maf2db_tonly"]["cpus-per-task"]
    conda:
        "../envs/variantannotation_r.yaml"
    log:
        "logs/mutect2/tonly/maf2db/maf2db.log"
    benchmark:
        "benchmarks/mutect2/tonly/maf2db/maf2db.txt"
    message:
        "Copy variants to TSV for uploading to database"
    script:
        "../R/snakemake/snv2db.R"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Annotate consensus variant set (with VEP via vcf2maf)
## Excluding the normal-id field here:
##            --tumor-id TUMOR \
##            --normal-id NORMAL \
## Make sure that removing the normal-id doesn't break the script
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule annotate_vep_tonly:
    input:
        vcf = "results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vcf.gz"
    output:
        vcf_uncompressed = temp("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.uncompressed.vcf"),
        vcf_reformatted = temp("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.reformatted.vcf"),
        maf = protected("results/mutect2/tonly/consensusvcf/consensus.normalized.sorted.vep.maf")
    params:
        mem = CLUSTER_META["annotate_vep_tonly"]["mem"]
    threads:
        CLUSTER_META["annotate_vep_tonly"]["cpus-per-task"]
    conda:
        "../envs/vcf2maf.yaml"
    log:
        "logs/mutect2/tonly/annoconsensusvcf/annoconsensusvcf_vep.log"
    benchmark:
        "benchmarks/mutect2/tonly/annoconsensusvcf/annoconsensusvcf_vep.txt"
    message:
        "Annotate consensus variants using VEP"
    shell:
        "bcftools view \
            -Ov \
            -o {output.vcf_uncompressed} {input.vcf} \
            > {log} 2>&1;"

        "vcf2vcf.pl \
            --input-vcf {output.vcf_uncompressed} \
            --output-vcf {output.vcf_reformatted} \
            --ref-fasta {config[reference_fasta]} \
            >> {log} 2>&1;"
        
        "vcf2maf.pl \
            --input-vcf {output.vcf_reformatted} \
            --output-maf {output.maf} \
            --vep-path {config[veppath]} \
            --vep-data {config[vepdata]} \
            --vep-forks 8 \
            --ref-fasta {config[reference_fasta]} \
            --filter-vcf {config[gnomad_vcf]} \
            --tumor-id TUMOR \
            --normal-id NORMAL \
            --species homo_sapiens \
            --ncbi-build GRCh37 \
            >> {log} 2>&1"
            
#Manually performing this step using vep_upload R script
#
#rule vep2db:
#    input:
#        maf = "results/mutect2/annoconsensusvcf/consensus.normalized.sorted.vep.maf",
#        vcf = "results/mutect2/consensusvcf/consensus.normalized.sorted.vcf.gz"
#    output:
#        tsv = "results/mutect2/maf2db/consensus.normalized.sorted.vep.tsv"
#    params:
#        mem = CLUSTER_META["vep2db"]["mem"]
#    threads:
#        CLUSTER_META["vep2db"]["cpus-per-task"]
#    conda:
#        "../envs/vep2db.yaml"
#    log:
#        "logs/mutect2/vep2db/vep2db.log"
#    benchmark:
#        "benchmarks/mutect2/vep2db/vep2db.txt"
#    message:
#        "Copy variants to remote"
#    script:
#        "../R/snakemake/vep_upload.R"
            
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Upload genotype to database
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule geno2db_tonly:
    input:
        vcf = "results/mutect2/tonly/dropgt/{case_barcode}.filtered.normalized.sorted.vcf.gz",
        ssvcf = lambda wildcards: ancient(expand("results/mutect2/tonly/ssdropgt/{aliquot_barcode}.filtered.normalized.sorted.vcf.gz", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode)))
    output:
        geno = "results/mutect2/tonly/geno2db/{case_barcode}.geno.tsv",
        info = "results/mutect2/tonly/geno2db/{case_barcode}.info.tsv",
        mfcov = "results/mutect2/tonly/geno2db/{case_barcode}.mfcov.tsv",
    params:
        mem = CLUSTER_META["geno2db_tonly"]["mem"]
    threads:
        CLUSTER_META["geno2db_tonly"]["cpus-per-task"]
    conda:
        "../envs/variantannotation_r.yaml"
    log:
        "logs/mutect2/tonly/geno2db/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/geno2db/{case_barcode}.txt"
    message:
        "Copy M2 calls to TSV for uploading to database\n"
        "Case: {wildcards.case_barcode}"
    script:
        "../R/snakemake/geno2db.R"
     
# ## END ##