##################################
### Upload tumor-only files data for CARE WXS files
### Author: Kevin Johnson
### Date updated: 2022.07.08
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
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

## Derive a list of tumor-only aliquots
# Identify all tumor aliquots
care_tumor_aliquots <- aliquots %>% 
  filter(grepl("CARE", aliquot_barcode), !grepl("-NB", aliquot_barcode))
# Identify all tumor aliquots with a matched normal blood sample
care_tumor_aliquots_matched <- pairs %>% 
  filter(grepl("CARE", tumor_barcode))

# Restrict the analysis to those without a normal blood sample.
care_tumor_aliquots_tonly <- care_tumor_aliquots %>% 
  filter(!aliquot_barcode%in%care_tumor_aliquots_matched$tumor_barcode) %>% 
  # excluding CARE-TK-TK06 which was already uploaded to the database
  filter(!grepl("CARE-TK-TK06", aliquot_barcode))

# Examine the files to make sure that it is indeed correct (i.e., missing 'aligned bam' designation)
care_tonly_files <- files %>% 
  filter(aliquot_barcode%in%care_tumor_aliquots_tonly$aliquot_barcode) %>% 
  filter(file_format=="aligned BAM")

## CARE WXS file size and md5 files were determined on sumner.
realn_size <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/care_wxs_tonly_realigned_filesize.txt", header=FALSE)
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

# The France cohort
realn_md5_ps <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/care_ps_tonly_realigned_md5.txt", sep=" ", header=FALSE)
realn_md5_ps <- realn_md5_ps[ ,c(1,3)]
colnames(realn_md5_ps) <- c("file_md5sum", "file_path")
realn_md5_ps$file_name <- sapply(strsplit(as.character(realn_md5_ps$file_path), "/"), last)

# The MD Anderson cohort
realn_md5_md <- read.delim("/projects/verhaak-lab/USERS/johnsk/glass4/results/db_uploads/care_md_tonly_realigned_md5.txt", sep=" ", header=FALSE)
realn_md5_md <- realn_md5_md[ ,c(1,3)]
colnames(realn_md5_md) <- c("file_md5sum", "file_path")
realn_md5_md$file_name <- sapply(strsplit(as.character(realn_md5_md$file_path), "/"), last)

# Combine together
realn_md5 <- bind_rows(realn_md5_md, realn_md5_ps)

## Create re-aligned file list to be uploaded to the database so that a pairs table can be created.
care_realn_files <- realn_size %>% 
  inner_join(realn_md5, by=c("file_name", "file_path")) %>% 
  dplyr::select(aliquot_barcode, file_name, file_size, file_md5sum, file_format, file_path)

# Filter to those tumor-only aliquots (excluding CARE-TK-TK06 which was already included)
care_realn_files_filt <- care_realn_files %>% 
  filter(aliquot_barcode%in%care_tumor_aliquots_tonly$aliquot_barcode)

## Upload these data to the glass4 database.
dbWriteTable(con, Id(schema="analysis", table="files"), care_realn_files_filt, append=T)

#### END ###### 