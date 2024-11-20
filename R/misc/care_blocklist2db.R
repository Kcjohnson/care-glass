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

## Retrieve the blocklist table and CARE aliquots.
blocklist = dbReadTable(con,  Id(schema="analysis",table="blocklist"))
aliquots = dbReadTable(con,  Id(schema="biospecimen",table="aliquots"))
surgeries = dbReadTable(con,  Id(schema="clinical",table="surgeries"))

## What's the minimum surgical interval? If longer than 2 months, clinical_exclusion = "allow".
surgeries %>% 
  filter(grepl("CARE-", case_barcode), surgical_interval_mo!=0) %>% 
  summarise(min_interval = min(surgical_interval_mo, na.rm = T))

## Create CARE blocklist.
care_aliquot_blocklist <- aliquots %>% 
  filter(grepl("CARE-", aliquot_barcode)) %>% 
  dplyr::select(aliquot_barcode) %>% 
  left_join(blocklist, by="aliquot_barcode") %>% 
  mutate(fingerprint_exclusion = "allow",
         coverage_exclusion = ifelse(grepl("NB", aliquot_barcode), NA, "allow"),
         cnv_exclusion = "allow",
         clinical_exclusion = "allow")

## Perform manual review of GATK CNV plots. 
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD03-R2-01D-WXS-K3WXX0"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD03-R2-01D-WXS-K3WXX0"] <- "insufficient_signal"
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD05-TP-01D-WXS-K9CF1W"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD05-TP-01D-WXS-K9CF1W"] <- "insufficient_signal"
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD07-R1-01D-WXS-S9K78J" ] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-MD-MD07-R1-01D-WXS-S9K78J" ] <- "insufficient_signal"

## CARE-PS had lower quality on average. Take caution when interpreting some of these CNV results.
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR03-R1-01D-WXS-YIGKHB"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR03-R1-01D-WXS-YIGKHB"] <- "noisy_signal"
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR04-TP-01D-WXS-ELWNZ7"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR04-TP-01D-WXS-ELWNZ7"] <- "noisy_signal"
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR04-R1-01D-WXS-95SIV0" ] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR04-R1-01D-WXS-95SIV0"] <- "noisy_signal"
care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR07-R1-01D-WXS-HWNFNP"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-PS-FR07-R1-01D-WXS-HWNFNP"] <- "noisy_signal"

care_aliquot_blocklist$cnv_exclusion[care_aliquot_blocklist$aliquot_barcode=="CARE-PT-DU05-R1-01D-WXS-I00DOZ"] <- "review"
care_aliquot_blocklist$cnv_exclusion_reason[care_aliquot_blocklist$aliquot_barcode=="CARE-PT-DU05-R1-01D-WXS-I00DOZ"] <- "insufficient_signal"

## Write the blocklist table to the database.
dbWriteTable(con, Id(schema="analysis",table="blocklist"), care_aliquot_blocklist, append = TRUE)

### END ###