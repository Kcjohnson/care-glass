#################################
# Push sequenza seg into db
# Updated: 2022.03.01
# Author: Kevin Johnson
#################################
library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB4")

seg <- read.delim("results/sequenza/glass_seqz_segments.tsv", as.is = TRUE)
pp <- read.delim("results/sequenza/glass_seqz_purity_ploidy.tsv", as.is = TRUE)

pp <- pp %>% select(pair_barcode, cellularity, ploidy, slpp = SLPP)
seg <- seg %>% transmute(pair_barcode,
                         chrom = ifelse(chromosome=='X',23,as.integer(chromosome)),
                         pos = sprintf("[%s,%s]",start.pos,end.pos),
                         baf = Bf,
                         baf_n = N.BAF,
                         baf_sd = sd.BAF,
                         ratio = depth.ratio,
                         ratio_n = N.ratio,
                         ratio_sd = sd.ratio,
                         copy_number = CNt,
                         major_cn = A,
                         minor_cn = B,
                         log_posterior_proba = LPP)


dbWriteTable(con, Id(schema="variants",table="seqz_seg"), seg, append=T)
dbWriteTable(con, Id(schema="variants",table="seqz_params"), pp, append=T)

### END ###