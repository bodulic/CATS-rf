#!/usr/bin/env Rscript
#Script to add uncovered bases to pysamstats output
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TR_LEN_DT_PATH <- ext_args[1]

#Importing transcript lengths from file
tr_lengths <- fread(TR_LEN_DT_PATH, header = F)
setnames(tr_lengths, c("transcript", "tr_len"))

#Importing pysamstats output from STDIN
per_base_cov_acc <- fread('file:///dev/stdin', header = T)
setnames(per_base_cov_acc, c("transcript", "position", "coverage", "matches", "insertions"))

#Excluding transcripts not present in the imported pysamstats file
tr_lengths <- tr_lengths[transcript %in% unique(per_base_cov_acc[, transcript])]

#Creating a table with all bases
all_bases_dt <- list()
for (i in 1 : tr_lengths[, .N]) {
 all_bases_dt[[i]] <- data.table(transcript = tr_lengths[i, transcript], position = 1 : tr_lengths[i, tr_len])
}
all_bases_dt <- rbindlist(all_bases_dt)

#Adding uncovered bases to the original pysamstats output
all_bases_dt <- merge(all_bases_dt, per_base_cov_acc, all.x = T, by = c("transcript", "position"))
rm(per_base_cov_acc)
all_bases_dt[is.na(coverage), ':=' ("coverage" = 0, "matches" = 0, "insertions" = 0)]
setorder(all_bases_dt, transcript, position)

#Writing pysamstats output with uncovered bases to STDOUT
write.table(all_bases_dt, file = "", sep = "\t", row.names = F, col.names = F, quote = F)