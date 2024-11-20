library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(qusage) 

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
myinf2 <- "/projects/verhaak-lab/USERS/johnsk/glassx/varnf/data/msigdb/h.all.v7.1.symbols.gmt"

expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
genesets <- read.gmt(myinf2)

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

cat("Running...")
es <- gsva(expr, genesets, method="ssgsea")

myoutf <- "/fastscratch/johnsk/msigdb_hallmark_enrichment_idhwt.txt"
write.table(es, myoutf, sep="\t", quote=FALSE, row.names=TRUE)