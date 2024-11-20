#######################################################
# Add ssm2 coverage counts to db
# Date: 2023.03.22
# Author: Kevin Johnson
#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(DBI)
library(odbc)

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

files = list.files("results/mutect2/geno2db", recursive = T, pattern = ".mfcov.tsv", full.names = T)

## Double check to make sure the correct number of samples are returned.
files <- files[grep("GLSS-LX",files)]


for(i in 1:length(files))
{
	mfcov <- read.delim(files[i],header=FALSE,stringsAsFactor=FALSE)
	colnames(mfcov) <- c("aliquot_barcode","ad_depth","ssm2_call_count")
	dbWriteTable(con, Id(schema="variants", table="ssm2_count"), mfcov, append=TRUE)
}


### END ###