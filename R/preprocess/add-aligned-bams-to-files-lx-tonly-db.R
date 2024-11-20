##################################
### Upload tumor-only files data for GLSS-LX WGS files
### Author: Kevin Johnson
### Date updated: 2022.07.29
##################################
library(tidyverse)
library(odbc)
library(DBI)
library(stringi)

## Establish connection with GLASS 4 database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

## Load in tables that will need to be generated for the GLSS-LX data.
files = dbReadTable(con,  Id(schema="analysis",table="files"))
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

## Derive a list of tumor-only aliquots
# Identify all WGS aliquots
lx_tumor_aliquots <- aliquots %>% 
  filter(grepl("GLSS-LX", aliquot_barcode), grepl("-WGS-", aliquot_barcode))


# Examine the files to make sure that it is indeed correct (i.e., missing 'aligned bam' designation)
lx_tonly_files <- files %>% 
  filter(aliquot_barcode%in%lx_tumor_aliquots$aliquot_barcode) %>% 
  filter(file_format=="aligned BAM")

## GLSS-LX WGS file size and md5 files were determined on sumner.
realn_size <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/glss_lx_realigned_filesize.txt", header=FALSE)
colnames(realn_size) <- c("file_size", "file_path")
realn_size <- realn_size %>% 
  filter(grepl("\\/projects\\/verhaak-lab\\/USERS\\/johnsk\\/glass4\\/results\\/align\\/bqsr\\/", file_path)) %>% 
  distinct() %>% 
  arrange(file_path)
# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
realn_size$file_name <- sapply(strsplit(as.character(realn_size$file_path), "/"), last)
realn_size$aliquot_barcode <- sapply(strsplit(realn_size$file_name, "\\."), "[", 1)
realn_size$file_format <- "aligned BAM"

# The GLSS-LX cohort
realn_md5_lx <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/glss_lx_realigned_md5.txt", sep=" ", header=FALSE)
realn_md5_lx <- realn_md5_lx[ ,c(1,3)]
colnames(realn_md5_lx) <- c("file_md5sum", "file_path")
realn_md5_lx$file_name <- sapply(strsplit(as.character(realn_md5_lx$file_path), "/"), last)



## Create re-aligned file list to be uploaded to the database so that a pairs table can be created.
lx_realn_files <- realn_size %>% 
  inner_join(realn_md5_lx, by=c("file_name", "file_path")) %>% 
  dplyr::select(aliquot_barcode, file_name, file_size, file_md5sum, file_format, file_path)

## Upload these data to the glass4 database.
dbWriteTable(con, Id(schema="analysis", table="files"), lx_realn_files, append=T)

#### END ###### 