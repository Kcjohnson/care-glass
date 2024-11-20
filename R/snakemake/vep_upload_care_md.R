# This script takes the output of the annotate_vep rule in the Snakemake mutect2-post.smk module reformats it for uploading to the db (variants.vep table)
# Additionally generates a .tsv file for backup
#-----------------------------------------------------

library(VariantAnnotation) # I believe this is only needed for the readVcf function.
library(ensemblVEP)
library(tidyverse)
library(DBI)
library(odbc)


## Parse snakemake. Change folder depending on cohort (e.g, CARE-MD/CARE-PT/CARE-PS).
maff = "results/mutect2/consensusvcf/CARE-MD/consensus.normalized.sorted.vep.maf"
vcff = "results/mutect2/consensusvcf/CARE-MD/consensus.normalized.sorted.vcf.gz"
tsvf = "results/mutect2/maf2db/CARE-MD/consensus.normalized.sorted.vep.tsv"

vcf = readVcf(vcff, "hg19")
maf = read.delim(maff, as.is = T, comment.char = '#')

message("Read file ", basename(vcff))
message("Read file ", basename(maff))

df = data.frame(chrom = as.character(seqnames(vcf)),
				pos = sprintf("[%s,%s)", start(vcf), end(vcf)+1),
                ref = ref(vcf),
                alt = unstrsplit(CharacterList(alt(vcf)), sep=","),
                gene_id = maf$Gene,
                gene_symbol = maf$Hugo_Symbol,
                variant_classification = maf$Variant_Classification,
                variant_type = maf$Variant_Type,
                cdna_position = maf$cDNA_position,
                cds_position = maf$CDS_position,
                protein_position = maf$Protein_position,
                amino_acids = maf$Amino_acids, 
                codons = maf$Codons,
                hgvs_c = maf$HGVSc,
                hgvs_p = maf$HGVSp_Short,
                polyphen = maf$PolyPhen,
                sift = maf$SIFT,
                stringsAsFactors = F)

#Change chromosome X to chromosome 23
df[which(df[,"chrom"]=='X'),"chrom"] <- 23
df[,"chrom"] = as.numeric(df[,"chrom"])

# Manual edit to match GLASS variant_classifications table; this is now done in SQL
df[which(df[,"variant_classification"]=="Splice_Region"),"variant_classification"]  <- "Splice_Site"


# Remove redundant rows
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

old <- dbReadTable(con, Id(schema = "variants",table="vep"))
vid <- dbReadTable(con, Id(schema = "variants",table="anno"))

new <- df %>%
	   inner_join(select(vid, variant_id, chrom, pos, ref, alt), by = c("chrom", "pos", "ref", "alt")) %>%
	   anti_join(old, "variant_id")
new <- new[,c(ncol(new), 1:(ncol(new)-1))]	   

full <- rbind(old, new)
write.table(full, file = tsvf, quote = F, sep = "\t", row.names = F, col.names = T)
message("Wrote output ", basename(tsvf))

dbWriteTable(con, Id(schema="variants",table="vep"), new, append=TRUE)

### END ###