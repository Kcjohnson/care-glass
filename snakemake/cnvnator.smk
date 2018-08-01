## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNVnator does not have a good conda installation and needs to be installed manually
## This scripts expects a working install of "cnvnator" command and cnvnator2VCF.pl
## script
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extracting read mappings from BAM file
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_tree:
    input:
        "results/align/bqsr/{aliquot_id}.realn.mdup.bqsr.bam"
    output:
        protected("results/cnvnator/tree/{aliquot_id}.tree.root")
    params:
        mem = CLUSTER_META["cnvnator_tree"]["mem"]
    threads:
        CLUSTER_META["cnvnator_tree"]["ppn"]
    log:
        "logs/cnvnator/tree/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/tree/{aliquot_id}.txt"
    message:
        "Extracting read mappings from BAM file\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "cnvnator -root {output} \
            -tree {input} \
            -chrom {config[cnvnator_chrom]} \
            -unique \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Generating a histogram
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_his:
    input:
        "results/cnvnator/tree/{aliquot_id}.tree.root"
    output:
        temp("results/cnvnator/his/{aliquot_id}.tree.his.root")
    params:
        binsize = 1000,
        contigs = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y",
        mem = CLUSTER_META["cnvnator_his"]["mem"]
    threads:
        CLUSTER_META["cnvnator_his"]["ppn"]
    log:
        "logs/cnvnator/his/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/his/{aliquot_id}.txt"
    message:
        "Generating a histogram\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "(cnvnator -root {input} \
            -his {config[cnvnator_binsize]} \
            -d {config[cnvnator_refdir]} && \
            rsync {input} {output}) \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calculating statistics
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_stat:
    input:
        "results/cnvnator/his/{aliquot_id}.tree.his.root"
    output:
        temp("results/cnvnator/stat/{aliquot_id}.tree.his.stat.root")
    params:
        mem = CLUSTER_META["cnvnator_stat"]["mem"]
    threads:
        CLUSTER_META["cnvnator_stat"]["ppn"]
    log:
        "logs/cnvnator/stat/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/stat/{aliquot_id}.txt"
    message:
        "Calculating statistics\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "(cnvnator -root {input} \
            -stat {config[cnvnator_binsize]} && \
            rsync {input} {output}) \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## RD Signal Partitioning
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_partition:
    input:
        "results/cnvnator/stat/{aliquot_id}.tree.his.stat.root"
    output:
        protected("results/cnvnator/partition/{aliquot_id}.tree.his.stat.partition.root")
    params:
        mem = CLUSTER_META["cnvnator_partition"]["mem"]
    threads:
        CLUSTER_META["cnvnator_partition"]["ppn"]
    log:
        "logs/cnvnator/partition/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/partition/{aliquot_id}.txt"
    message:
        "RD Signal Partitioning\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "(cnvnator -root {input} \
            -partition {config[cnvnator_binsize]} && \
            rsync {input} {output}) \
            > {log} 2>&1"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CNV calling
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator_call:
    input:
        "results/cnvnator/partition/{aliquot_id}.tree.his.stat.partition.root"
    output:
        protected("results/cnvnator/call/{aliquot_id}.call.tsv")
    params:
        mem = CLUSTER_META["cnvnator_call"]["mem"]
    threads:
        CLUSTER_META["cnvnator_call"]["ppn"]
    log:
        "logs/cnvnator/call/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/call/{aliquot_id}.txt"
    message:
        "CNV calling\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "cnvnator -root {input} \
            -call {config[cnvnator_binsize]} \
            >2 {log} 1> {output}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Convert CNVnator calls TSV to VCF
## See: https://github.com/abyzovlab/CNVnator
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

rule cnvnator2vcf:
    input:
        "results/cnvnator/call/{aliquot_id}.call.tsv"
    output:
        "results/cnvnator/vcf/{aliquot_id}.call.vcf"
    params:
        mem = CLUSTER_META["cnvnator2vcf"]["mem"]
    threads:
        CLUSTER_META["cnvnator2vcf"]["ppn"]
    log:
        "logs/cnvnator/vcf/{aliquot_id}.log"
    benchmark:
        "benchmarks/cnvnator/vcf/{aliquot_id}.txt"
    message:
        "Convert CNVnator calls TSV to VCF\n"
        "Sample: {wildcards.aliquot_id}"
    shell:
        "cnvnator2VCF.pl {input} \
            {config[cnvnator_refdir]} \
            >2 {log} 1> {output}"

## END ##