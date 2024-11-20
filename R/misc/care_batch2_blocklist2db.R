#################################
## Create blocklist for batch2 of CARE
## Updated: 2022.07.05
## Author: Kevin Johnson
#################################

library(DBI)
library(tidyverse)

## Connect to the latest version of the database that contains the CARE samples.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

## Load that tables that will be needed, including the blocklist that will be edited.
blocklist = dbReadTable(con,  Id(schema="analysis",table="blocklist"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))

## Take the data introduced in batch2 CARE - cohorts CARE-TO and CARE-TU.
aliquots_batch2 <- aliquots %>% 
  filter(grepl("CARE-TO|CARE-TK", aliquot_barcode))

## Create an empty blocklist data.frame.
blocklist_batch2 <- data.frame(matrix(ncol = dim(blocklist)[2], nrow = dim(aliquots_batch2)[1]))
colnames(blocklist_batch2) <- colnames(blocklist)
blocklist_batch2$aliquot_barcode <- aliquots_batch2$aliquot_barcode

## Start with clinical exclusion
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))
surgeries_batch2 <- surgeries %>% 
  filter(grepl("CARE-TO|CARE-TK", case_barcode))

## Check to see whether any surgical intervals were less than 2 months, which may indicate not a true recurrence.
table(surgeries_batch2$surgical_interval_mo) 
## There are no clinical exclusions
blocklist_batch2$clinical_exclusion <- "allow"

## There were no reported fingerprint exclusions for these batches,
blocklist_batch2$fingerprint_exclusion <- "allow"

## Manual review
## Enter default of "allow" for all fields. 
blocklist_batch2$cnv_exclusion <- "allow"


## Change from "allow" default and specify reasoning. Available options pulled from current blocklist.
review_aliquots <- c("CARE-TK-TK01-R2-01D-WXS-J60B3L", "CARE-TK-TK04-R1-01D-WXS-TG719Y",
                     "CARE-TK-TK05-R1-01D-WXS-6LPEAD",
                     "CARE-TO-SM02-R1-01D-WXS-IS4Z60", "CARE-TO-SM02-TP-01D-WXS-YDNDU3",
                     "CARE-TO-SM06-TP-01D-WXS-QIS1HH")
## All the samples in this list are likely of low tumor purity and it's difficult to observe CNV events.
blocklist_batch2$cnv_exclusion[blocklist_batch2$aliquot_barcode%in%review_aliquots] <- "review"

## Provide the rationale. In some cases there appears like there is likely "insufficient_signal" and in other cases, there is not even a faint CNV signal.
review_rationale <- c("insufficient_signal", "unexpected_genome_stability",
                      "insufficient_signal",
                      "insufficient_signal", "unexpected_genome_stability",
                      "unexpected_genome_stability")
blocklist_batch2$cnv_exclusion_reason[blocklist_batch2$aliquot_barcode%in%review_aliquots] <- review_rationale

## There are no samples with low sequencing coverage in these batches.
blocklist_batch2$coverage_exclusion <- "allow"


## Upload these data to the glass4 database.
dbWriteTable(con, Id(schema="analysis", table="blocklist"), blocklist_batch2, append=T)


### END ###
