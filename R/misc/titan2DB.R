#####################################
## Write TITAN segments to database
## Author: Kevin Johnson
## Updated: 2022.02.25
#####################################

library(tidyverse)
library(DBI)
library(odbc)

## Current project.
projdir <- "/projects/verhaak-lab/USERS/johnsk/glass4"
setwd(projdir)

## Connect to the correct* internal database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")
pairs = dbReadTable(con,  Id(schema="analysis",table="pairs"))

## The segment files outputted by TITAN
segfiles <- list.files('results/cnv/titanfinal/seg', full.names = TRUE)

## Query the database to get a table of pairs.
q <- "SELECT * FROM analysis.pairs"
pair <- dbGetQuery(con, q)

## Write each new TITAN segment file to the database. Note that the append argument is used here.
lapply(segfiles, function(f){
  message(f)
  dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
  df <- dat %>%
    transmute(pair_barcode = Sample,
              chrom = Chromosome,
              pos = sprintf("[%s,%s]",Start_Position.bp.,End_Position.bp.),
              num_snp = Length.snp.,
              median_ratio = Median_Ratio,
              median_logr = Median_logR,
              titan_state = TITAN_state,
              titan_call = TITAN_call,
              copy_number = Copy_Number,
              major_cn = MajorCN,
              minor_cn = MinorCN,
              clonal_cluster = Clonal_Cluster,
              cellular_prevalence = Cellular_Prevalence,
              logr_copy_number = logR_Copy_Number,
              corrected_copy_number = Corrected_Copy_Number,
              corrected_call = Corrected_Call) %>% 
    inner_join(pairs, by="pair_barcode") %>% 
    dplyr::select(pair_barcode:corrected_call, aliquot_barcode = tumor_barcode)
	
	df[which(df[,"chrom"]=="X"),"chrom"] <- 23
	df[which(df[,"chrom"]=="Y"),"chrom"] <- 24

    dbWriteTable(con, Id(schema="variants",table="titan_seg"), df, append=T)
    Sys.sleep(1)
})

### END ###