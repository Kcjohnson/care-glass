#####################################
## Write GATK segments to database
## Author: Kevin Johnson
## Updated: 2022.02.25
#####################################

library(DBI)
library(tidyverse)

## Current project
projdir <- "/projects/verhaak-lab/USERS/johnsk/glass4"
setwd(projdir)

## Connect to the correct* internal database.
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

## Find segments for aliquots. If directory already has segments on the database be sure to select new files.
segfiles <- list.files("results/cnv/callsegments", full.names = TRUE)
# Ignore the igv seg files.
segfiles <- segfiles[grep(".called.seg",segfiles)]


segs <- parallel::mclapply(segfiles, function(f) {
  dat <- read.delim(f, comment.char = "@", as.is= TRUE)
  dat <- dat %>%
    mutate(aliquot_barcode = substr(basename(f),1,30), pos = sprintf("[%s,%s]", START, END)) %>%
    select(aliquot_barcode, chrom = CONTIG, pos, num_points = NUM_POINTS_COPY_RATIO, log2_copy_ratio = MEAN_LOG2_COPY_RATIO, cnv_call = CALL)
  return(dat)
}, mc.cores = 8)
segs <- data.table::rbindlist(segs) %>% as.data.frame()

# If the data in the directory is mixed new and old, add the new cohort (e.g., CARE).
# Change this grep to be whatever new string you would like to select.
new_segs <- segs[grep("CARE-",segs[,"aliquot_barcode"]),]

# Change sex chromsomes to numbers.
new_segs[which(new_segs[,"chrom"] == "X"),"chrom"] <- 23
new_segs[which(new_segs[,"chrom"] == "Y"),"chrom"] <- 24

# Does this number match the expected number of aliquots?
n_distinct(new_segs$aliquot_barcode)

# Check for conflicts with the current aliquots table
q <-"SELECT * FROM biospecimen.aliquots WHERE aliquot_batch LIKE 'CARE-%'"
aliquots <- dbGetQuery(con, q)

# If aliquots are already present and you don't want to override, you can do something along these lines to select for the batch of interest.
new_segs_final <- new_segs[which(new_segs[,"aliquot_barcode"] %in% aliquots[,"aliquot_barcode"]),]

## Write out to the database.
dbWriteTable(con, Id(schema="variants",table="gatk_seg"), new_segs_final, append=TRUE)


### END ###