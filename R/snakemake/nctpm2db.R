#!/usr/bin/env Rscript

#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

transcript_tpm <- args[1]

cat("Reading in inputs...")
ttpm <- read.delim(transcript_tpm,stringsAsFactors=FALSE)

output <- "/projects/verhaak-lab/USERS/johnsk/glass4/results/kallisto/noncoding/analysis_noncoding_tpm.txt"
write.table(ttpm,output,sep="\t",quote=FALSE,row.names=FALSE)

#dbWriteTable(con, Id(schema="analysis", table="noncoding_tpm"), ttpm, overwrite=TRUE)

### END ###