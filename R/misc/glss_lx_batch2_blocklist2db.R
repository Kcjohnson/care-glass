#################################
## Create blocklist for GLSS-LX WGS batch 2
## Updated: 2023.03.20
## Author: Kevin Johnson
#################################

library(DBI)
library(tidyverse)

## Connect to the latest version of the database that contains the GLSS-LX samples.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

## Load that tables that will be needed, including the blocklist that will be edited.
blocklist = dbReadTable(con,  Id(schema="analysis",table="blocklist"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

## Some GLSS-LX samples are already in the blocklist
ext_blocklist_lx <- blocklist %>% 
  filter(grepl("GLSS-LX", aliquot_barcode))

## Take the data introduced in batch2 GLSS-LX
aliquots_lx <- aliquots %>% 
  filter(grepl("GLSS-LX", aliquot_barcode), aliquot_analysis_type=="WGS") %>% 
  filter(!aliquot_barcode%in%ext_blocklist_lx$aliquot_barcode)

## Create an empty blocklist data.frame.
blocklist_lx <- data.frame(matrix(ncol = dim(blocklist)[2], nrow = dim(aliquots_lx)[1]))
colnames(blocklist_lx) <- colnames(blocklist)
blocklist_lx$aliquot_barcode <- aliquots_lx$aliquot_barcode

## Start with clinical exclusion
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))
surgeries_lx <- surgeries %>% 
  filter(grepl("GLSS-LX", case_barcode), !grepl("GLSS-LX-0654", case_barcode))

## Check to see whether any surgical intervals were less than 2 months, which may indicate not a true recurrence.
table(surgeries_lx$surgical_interval_mo) 
## There are no clinical exclusions
blocklist_lx$clinical_exclusion <- "allow"

## There were no reported fingerprint exclusions for this batch.
blocklist_lx$fingerprint_exclusion <- "allow"

## Manual review
## Enter default of "allow" for all fields. 
blocklist_lx$cnv_exclusion <- "allow"


## Change from "allow" default and specify reasoning. Available options pulled from current blocklist.
block_aliquots <- c("GLSS-LX-0409-TP-01D-WGS-QO2ESR", "GLSS-LX-0677-R2-01D-WGS-KXWO0E", "GLSS-LX-0833-R1-01D-WGS-N4818R", "GLSS-LX-0959-R1-01D-WGS-LSTK3U")

## All the samples in this list are likely of low tumor purity and it's difficult to observe CNV events.
blocklist_lx$cnv_exclusion[blocklist_lx$aliquot_barcode%in%block_aliquots] <- "block"

## Provide the rationale. In some cases there appears like there is likely "insufficient_signal" and in other cases, there is not even a faint CNV signal.
review_rationale <- c("unexpected_genome_stability")
blocklist_lx$cnv_exclusion_reason[blocklist_lx$aliquot_barcode%in%block_aliquots] <- review_rationale

## There are no samples with low sequencing coverage in these batches.
blocklist_lx$coverage_exclusion <- "allow"

blocklist_lx %>% 
  filter(cnv_exclusion=="block")

## Upload these data to the glass4 database. Successfully done on 2023.03.20.
dbWriteTable(con, Id(schema="analysis", table="blocklist"), blocklist_lx, append=T)


### END ###
