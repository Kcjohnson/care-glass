#######################################################
# Necessary packages:
library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

## Current project.
projdir <- "/projects/verhaak-lab/USERS/johnsk/glass4"
setwd(projdir)

## Connect to the correct* internal database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

# For the GLSS-LX cohort remove the aliquots already in the db.
ext_gencode_coverage = dbReadTable(con,  Id(schema="analysis",table="gencode_coverage"))
ext_gencode_aliquots <- unique(ext_gencode_coverage$aliquot_barcode)[grep("GLSS-LX",unique(ext_gencode_coverage$aliquot_barcode))]
to_remove <- paste0(ext_gencode_aliquots, ".gencode-coverage.tsv")

myDir1 <- "results/align/gencode-coverage/"
myinf_tmp <- dir(myDir1)

#Specify cohort of interest if adding one specific cohort
#myinf1 <- myinf1[grep("GLSS-SN-",myinf1)]
#myinf1 <- myinf1[grep("CARE-TK|CARE-TO",myinf1)]
myinf_tmp <- myinf_tmp[grep("GLSS-LX",myinf_tmp)]
# Subset and remove elements
myinf1 <- myinf_tmp[!(myinf_tmp %in% to_remove)]


aliquot_barcode <- gsub(".gencode-coverage.tsv","",myinf1)
myinf1 <- paste(myDir1,myinf1,sep="")

dat <- list()
for(i in 1:length(myinf1))
{
	aliquot_gencode <- read.delim(myinf1[i],header=FALSE)
	aliquot_gencode <- cbind(aliquot_barcode[i], aliquot_gencode)
	dat[[i]] <- aliquot_gencode
}

new_gencode <- do.call(rbind,dat)

#Remove sum of exon sizes column:
final_gencode <- new_gencode[,c(1,2,4)]
colnames(final_gencode) <- c("aliquot_barcode","ensembl_gene_id","gene_coverage")


message("starting upload of analysis.gencode_coverage to the database.")

## Write the genecode_coverage tables to the database.
dbWriteTable(con, Id(schema="analysis",table="gencode_coverage"), final_gencode, append=T)

message("Completed upload of analysis.gencode_coverage to the database.")

### END ###
