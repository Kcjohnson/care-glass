##################################
### Upload new files data for CARE WXS files
### Author: Kevin Johnson
### Date updated: 2022.02.18
##################################
library(tidyverse)
library(odbc)
library(DBI)
library(stringi)

## Establish connection with GLASS 4 database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

## Load in tables that will need to be generated for the CARE data.
files = dbReadTable(con,  Id(schema="analysis",table="files"))
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))

## CARE WXS file size and md5 files were determined on sumner.
realn_size <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/care_batch2_realigned_filesize.txt", header=FALSE)
colnames(realn_size) <- c("file_size", "file_path")
# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
realn_size$file_name <- sapply(strsplit(as.character(realn_size$file_path), "/"), last)
realn_size$aliquot_barcode <- sapply(strsplit(realn_size$file_name, "\\."), "[", 1)
realn_size$file_format <- "aligned BAM"

realn_md5 <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/care_batch2_realigned_md5.txt", sep=" ", header=FALSE)
realn_md5 <- realn_md5[ ,c(1,3)]
colnames(realn_md5) <- c("file_md5sum", "file_path")
realn_md5$file_name <- sapply(strsplit(as.character(realn_md5$file_path), "/"), last)


## Create re-aligned file list to be uploaded to the database so that a pairs table can be created.
care_realn_files <- realn_size %>% 
  inner_join(realn_md5, by=c("file_name", "file_path")) %>% 
  dplyr::select(aliquot_barcode, file_name, file_size, file_md5sum, file_format, file_path)

## Upload these data to the glass4 database.
dbWriteTable(con, Id(schema="analysis", table="files"), care_realn_files, append=T)

#### END ###### 