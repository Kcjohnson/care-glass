## Changes 2022.07.07
## - Adding tumor-only mutation calling for CARE project.
## - Removed the tumor-segmentation argument in CalculateContamination since we are only using tumor-only.
## - To enable comparisons with matched normal samples, we will place all results for each cohort in results/mutect2/tonly/ followed by the same structure as the matched normal blood versions.
## - This utilizes existing panel of normals from the same aliquot_batch as the tumor-only samples. May need to adjust if there are ever tumor-only batches.

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Call SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This step uses Mutect2 to call variants on a tumor-only sample, using a panel-of-normals
## constructed from germline samples that were part of the same aliquot_batch (i.e., derived from same population and shared same sequencing runs)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## GATK parameters taken from "GATK4_SomaticSNVindel_worksheet.pdf"
## Modified from GLASS rule callsnv
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule callsnv_tonly:
    input:
        tumor = lambda wildcards: ancient(expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode))),
        pon = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf".format(aliquot_batch = manifest.getBatchByCase(wildcards.case_barcode))),
        ponidx = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf.idx".format(aliquot_batch = manifest.getBatchByCase(wildcards.case_barcode))),
        intervallist = lambda wildcards: ancient("{dir}/{interval}/scattered.interval_list".format(dir = config["mutect2"]["wgs_scatterdir"], interval = wildcards.interval)),
        orientation_priors = lambda wildcards: ancient(expand("results/mutect2/tonly/filterorientation/{aliquot_barcode}.priors.tsv", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode))),
    output:
        vcf = temp("results/mutect2/tonly/m2vcf-scatter/{case_barcode}.{interval}.vcf"),
        idx = temp("results/mutect2/tonly/m2vcf-scatter/{case_barcode}.{interval}.vcf.idx"),
        bam = temp("results/mutect2/tonly/m2bam-scatter/{case_barcode}.{interval}.bam")
    params:
        mem = CLUSTER_META["callsnv_tonly"]["mem"],
        sample_paths = lambda _, input: " ".join(["-I " + s for s in input["tumor"]]),
        priors_paths = lambda _, input: " ".join(["--orientation-bias-artifact-priors " + s for s in input["orientation_priors"]])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["callsnv_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/callsnv/{case_barcode}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/tonly/callsnv/{case_barcode}.{interval}.txt"
    message:
        "Calling tumor-only SNVs (Mutect2)\n"
        "Case: {wildcards.case_barcode}\n"
        "Interval: {wildcards.interval}"
    shell:
        "gatk --java-options -Xmx{params.mem} Mutect2 \
            -R {config[reference_fasta]} \
            {params.sample_paths} \
            {params.priors_paths} \
            -L {input.intervallist} \
            --panel-of-normals {input.pon} \
            --germline-resource {config[mutect2][gnomad_vcf]} \
            --genotyping-mode GENOTYPE_GIVEN_ALLELES \
            --genotype-filtered-alleles true \
            --alleles {config[mutect2][given_alleles]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            -bamout {output.bam} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single-sample SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calls SNVs in a single aliquot rather than a patient in tumor-only mode
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule sscallsnv_tonly:
    input:
        tumor = ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam"),
        pon = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf".format(aliquot_batch = manifest.getBatch(wildcards.aliquot_barcode))),
        ponidx = lambda wildcards: ancient("results/mutect2/pon/{aliquot_batch}.vcf.idx".format(aliquot_batch = manifest.getBatch(wildcards.aliquot_barcode))),
        intervallist = lambda wildcards: ancient("{dir}/{interval}/scattered.interval_list".format(dir = config["mutect2"]["wgs_scatterdir"], interval = wildcards.interval)),
        orientation_priors = ancient("results/mutect2/tonly/filterorientation/{aliquot_barcode}.priors.tsv")
    output:
        vcf = temp("results/mutect2/tonly/ssm2vcf-scatter/{aliquot_barcode}.{interval}.vcf"),
        idx = temp("results/mutect2/tonly/ssm2vcf-scatter/{aliquot_barcode}.{interval}.vcf.idx")
    params:
        mem = CLUSTER_META["sscallsnv_tonly"]["mem"]
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["sscallsnv_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/sscallsnv/{aliquot_barcode}.{interval}.log"
    benchmark:
        "benchmarks/mutect2/tonly/sscallsnv/{aliquot_barcode}.{interval}.txt"
    message:
        "Single Sample (tumor-only) Calling SNVs (Mutect2)\n"
        "Aliquot: {wildcards.aliquot_barcode}\n"
        "Interval: {wildcards.interval}"
    shell:
        "gatk --java-options -Xmx{params.mem} Mutect2 \
            -R {config[reference_fasta]} \
            -I {input.tumor} \
            --orientation-bias-artifact-priors {input.orientation_priors} \
            -L {input.intervallist} \
            --panel-of-normals {input.pon} \
            --germline-resource {config[mutect2][gnomad_vcf]} \
            --genotyping-mode GENOTYPE_GIVEN_ALLELES \
            --genotype-filtered-alleles true \
            --alleles {config[mutect2][given_alleles]} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge SNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## This is an intermediate "gather" step. Like `callpon`, `callsnv` is run in scatter mode, 
## using an interval list and spread across 50 jobs for each normal sample, this step uses
## MergeVCF to merge all intermediate VCF files into one VCF per normal sample
## See: 
## https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule mergesnv_tonly:
    input:
        lambda wildcards: expand("results/mutect2/tonly/m2vcf-scatter/{case_barcode}.{interval}.vcf", case_barcode = wildcards.case_barcode, interval = WGS_SCATTERLIST)
    output:
        vcf = protected("results/mutect2/tonly/m2vcf/{case_barcode}.vcf"),
        idx = protected("results/mutect2/tonly/m2vcf/{case_barcode}.vcf.idx")
    params:
        mem = CLUSTER_META["mergesnv_tonly"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["mergesnv_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/mergesnv/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/mergesnv/{case_barcode}.txt"
    message:
        "Merging VCF files (M2)\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} MergeVcfs \
            {params.input_files} \
            -O {output.vcf} \
            > {log} 2>&1"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge SNV (single sample)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single sample merge
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssmergesnv_tonly:
    input:
        lambda wildcards: expand("results/mutect2/tonly/ssm2vcf-scatter/{aliquot_barcode}.{interval}.vcf", aliquot_barcode = wildcards.aliquot_barcode, interval = WGS_SCATTERLIST)
    output:
        protected("results/mutect2/tonly/ssm2vcf/{aliquot_barcode}.vcf")
    params:
        mem = CLUSTER_META["ssmergesnv_tonly"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["ssmergesnv_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/ssmergesnv/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/ssmergesnv/{aliquot_barcode}.txt"
    message:
        "Single Sample Merging VCF files (M2)\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} MergeVcfs \
            {params.input_files} \
            -O {output} \
            > {log} 2>&1"

rule mergem2bam_tonly:
    input:
        lambda wildcards: expand("results/mutect2/tonly/m2bam-scatter/{case_barcode}.{interval}.bam", case_barcode = wildcards.case_barcode, interval = WGS_SCATTERLIST)
    output:
        protected("results/mutect2/tonly/m2bam/{case_barcode}.bam")
    params:
        mem = CLUSTER_META["mergem2bam_tonly"]["mem"],
        input_files = lambda _, input: " ".join(["-I " + s for s in input])
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["mergem2bam_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/mergem2bam/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/mergem2bam/{case_barcode}.txt"
    message:
        "Merging BAM files (M2)\n"
        "Case: {wildcards.case_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} MergeSamFiles \
            {params.input_files} \
            -O {output} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Summarize read support for known variant sites
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule pileupsummaries_tonly:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        temp("results/mutect2/tonly/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt")
    params:
        mem = CLUSTER_META["pileupsummaries_tonly"]["mem"]
    conda:
        "../envs/gatk4.yaml"
    threads:
        CLUSTER_META["pileupsummaries_tonly"]["cpus-per-task"]
    log:
        "logs/mutect2/tonly/pileupsummaries/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/pileupsummaries/{aliquot_barcode}.txt"
    message:
        "Generating pileupsummaries\n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} GetPileupSummaries \
            -I {input} \
            -V {config[mutect2][tiny_vcf]} \
            -L {config[mutect2][tiny_vcf]} \
            -O {output} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Estimate contamination
## Input: pileup summaries table
## Output: contamination table
## This tool estimates contamination based on the signal from ref reads at hom alt sites
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule calculatecontamination_tonly:
    input:
        tumortable = ancient("results/mutect2/tonly/pileupsummaries/{aliquot_barcode}.pileupsummaries.txt")
    output:
        cont = "results/mutect2/tonly/contamination/{aliquot_barcode}.contamination.txt"
    params:
        mem = CLUSTER_META["calculatecontamination_tonly"]["mem"]
    threads:
        CLUSTER_META["calculatecontamination_tonly"]["cpus-per-task"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/tonly/calculatecontamination/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/calculatecontamination/{aliquot_barcode}.txt"
    message:
        "Computing contamination\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} CalculateContamination \
            -I {input.tumortable} \
            --output {output.cont} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## If variants have not been filtered, filter, else done
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filtermutect_tonly:
    input:
        vcf = ancient("results/mutect2/tonly/m2vcf/{case_barcode}.vcf"),
        tab = lambda wildcards: expand("results/mutect2/tonly/contamination/{aliquot_barcode}.contamination.txt", aliquot_barcode = manifest.getTumorByCase(wildcards.case_barcode))
    output:
        stats = protected("results/mutect2/tonly/m2filter/{case_barcode}.filterstats.tsv"),
        vcf = protected("results/mutect2/tonly/m2filter/{case_barcode}.filtered.vcf.gz"),
        tbi = protected("results/mutect2/tonly/m2filter/{case_barcode}.filtered.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["filtermutect_tonly"]["mem"],
        ttab = lambda _, input: " ".join(["--contamination-table " + s for s in input["tab"]]),
    threads:
        CLUSTER_META["filtermutect_tonly"]["cpus-per-task"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/tonly/filtermutect/{case_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/filtermutect/{case_barcode}.txt"
    message:
        "Filtering Mutect2 calls\n"
        "Case: {wildcards.case_barcode}"
    shell:    
        "gatk --java-options -Xmx{params.mem} FilterMutectCalls \
            -V {input.vcf} \
            {params.ttab} \
            --stats {output.stats} \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Single sample filter
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule ssfiltermutect_tonly:
    input:
        vcf = ancient("results/mutect2/tonly/ssm2vcf/{aliquot_barcode}.vcf"),
        tab = "results/mutect2/tonly/contamination/{aliquot_barcode}.contamination.txt"
    output:
        stats = protected("results/mutect2/tonly/ssm2filter/{aliquot_barcode}.filterstats.tsv"),
        vcf = protected("results/mutect2/tonly/ssm2filter/{aliquot_barcode}.filtered.vcf.gz"),
        tbi = protected("results/mutect2/tonly/ssm2filter/{aliquot_barcode}.filtered.vcf.gz.tbi")
    params:
        mem = CLUSTER_META["ssfiltermutect_tonly"]["mem"]
    threads:
        CLUSTER_META["ssfiltermutect_tonly"]["cpus-per-task"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/tonly/ssfiltermutect/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/ssfiltermutect/{aliquot_barcode}.txt"
    message:
        "Single-Sample (tumor-only) Filtering Mutect2 calls\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:    
        "gatk --java-options -Xmx{params.mem} FilterMutectCalls \
            -V {input.vcf} \
            --contamination-table {input.tab} \
            --stats {output.stats} \
            -O {output.vcf} \
            --seconds-between-progress-updates {config[seconds_between_progress_updates]} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Collect metrics on sequencing context artifacts
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule collectartifacts_tonly:
    input:
        ancient("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam")
    output:
        tab = temp("results/mutect2/tonly/artifacts/{aliquot_barcode}.alt.tsv"),
        ref = temp("results/mutect2/tonly/artifacts/{aliquot_barcode}.ref.metrics"),
        alt = temp("results/mutect2/tonly/artifacts/{aliquot_barcode}.alt.metrics")
    params:
        prefix = "results/mutect2/tonly/artifacts/{aliquot_barcode}",
        mem = CLUSTER_META["collectartifacts_tonly"]["mem"]
    threads:
        CLUSTER_META["collectartifacts_tonly"]["cpus-per-task"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/tonly/collectartifacts/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/collectartifacts/{aliquot_barcode}.txt"
    message:
        "Collecting sequencing artifact metrics\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} CollectF1R2Counts \
            -R {config[reference_fasta]} \
            -I {input} \
            -alt-table {output.tab} \
            -ref-hist {output.ref} \
            -alt-hist {output.alt} \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Filter by orientation bias
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule filterorientation_tonly:
    input:
        tab = "results/mutect2/tonly/artifacts/{aliquot_barcode}.alt.tsv",
        ref = "results/mutect2/tonly/artifacts/{aliquot_barcode}.ref.metrics",
        alt = "results/mutect2/tonly/artifacts/{aliquot_barcode}.alt.metrics"
    output:
        "results/mutect2/tonly/filterorientation/{aliquot_barcode}.priors.tsv"
    params:
        mem = CLUSTER_META["filterorientation_tonly"]["mem"]
    threads:
        CLUSTER_META["filterorientation_tonly"]["cpus-per-task"]
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/mutect2/tonly/filterorientation/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/mutect2/tonly/filterorientation/{aliquot_barcode}.txt"
    message:
        "Calculating orientation bias\n"
        "Aliquot: {wildcards.aliquot_barcode}"
    shell:
        "gatk --java-options -Xmx{params.mem} LearnReadOrientationModel \
            -alt-table {input.tab} \
            -ref-hist {input.ref} \
            -alt-hist {input.alt} \
            -O {output} \
            > {log} 2>&1"

## END ##