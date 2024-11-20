library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(qusage) 

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
myinf2 <- "/projects/verhaak-lab/USERS/johnsk/glassx/varnf/data/msigdb/h.all.v7.1.symbols.gmt"
myinf3 <- "/projects/verhaak-lab/monitor/reference/genelists/MP_malignant_named.csv"

expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
genesets <- read.gmt(myinf2)
# Load the CARE IDHwt project metadata:
metaprograms <- read.delim(myinf3, sep=",", header = TRUE)

# Convert metaprogram df to list.
mp_lst <- lapply(names(metaprograms), function(x) metaprograms[,x])
names(mp_lst) <- names(metaprograms)


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")
#Read in subtype/rnaseq silver set data
q <- "
SELECT * 
FROM analysis.rna_silver_set rs
JOIN clinical.subtypes cs ON cs.case_barcode = rs.case_barcode
WHERE idh_codel_subtype = 'IDHwt'
"
info <- dbGetQuery(con, q)

myaliquots <- c(info[,"tumor_barcode_a"], info[,"tumor_barcode_b"])
expr <- expr[,myaliquots]
expr <- data.matrix(expr)

# cat("Running...")
# es <- gsva(expr, genesets, method="ssgsea")

# myoutf1 <- "/fastscratch/johnsk/msigdb_hallmark_enrichment_idhwt.txt"
# write.table(es, myoutf1, sep="\t", quote=FALSE, row.names=TRUE)

# Now run for the metaprograms
cat("Running...")
es_mp <- gsva(expr, mp_lst, method="ssgsea")

myoutf2 <- "/fastscratch/johnsk/care_metaprograms_enrichment_idhwt.txt"
write.table(es_mp, myoutf2, sep="\t", quote=FALSE, row.names=TRUE)


