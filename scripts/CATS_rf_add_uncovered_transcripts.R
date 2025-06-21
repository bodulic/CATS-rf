#!/usr/bin/env Rscript
#Script to convert uncovered transcript coordinates to pysamstats output
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
UNCOV_TR_LEN_DT_PATH <- ext_args[2]
OUT_PREF <- ext_args[3]

#Importing uncovered transcript lengths from file
setDTthreads(THREAD_N)
uncov_tr_lengths <- fread(UNCOV_TR_LEN_DT_PATH, header = F)
setnames(uncov_tr_lengths, c("transcript", "tr_len"))

#Creating a pysamstats table with uncovered transcripts
uncov_bases_dt <- list()
for (i in 1 : uncov_tr_lengths[, .N]) {
 uncov_bases_dt[[i]] <- data.table(transcript = uncov_tr_lengths[i, transcript], position = 1 : uncov_tr_lengths[i, tr_len], coverage = 0, matches = 0, insertions = 0)
}
uncov_bases_dt <- rbindlist(uncov_bases_dt)

#Writing the pysamstats table with uncovered transcripts to file
write.table(uncov_bases_dt, file = paste(OUT_PREF, "uncov_trans", sep = "_"), sep = "\t", row.names = F, col.names = F, quote = F)