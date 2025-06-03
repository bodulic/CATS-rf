#!/usr/bin/env Rscript
#Read assignment script (single-end mode)
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
EDIT_DIST_FACT <- as.numeric(ext_args[1])
SEED <- as.numeric(ext_args[2])
OUT_PREF <- ext_args[3]

#Importing read mappings from STDIN
read_mappings <- fread('file:///dev/stdin', header = F)
setnames(read_mappings, c("read_id", "transcript", "row_N", "ed_dist"))

#Filtering read mappings based on edit distance
read_mappings[, "keep_mapping" := fifelse(ed_dist <= min(ed_dist) + EDIT_DIST_FACT, T, F), by = "read_id"]
read_mappings <- read_mappings[keep_mapping == T]
read_mappings[, c("ed_dist", "keep_mapping") := NULL]

#Calculating the number of transcripts to which each read maps
read_mappings <- unique(read_mappings, by = c("read_id", "transcript"))
read_mappings[, "read_map_N" := .N, by = "read_id"]

#Saving reads mapped to a single transcript and continuing with multimapped reads
single_mapped_read_row_N <- read_mappings[read_map_N == 1, row_N]
read_mappings <- read_mappings[read_map_N != 1]
read_mappings[, "read_map_N" := NULL]

if (read_mappings[, .N] != 0) {
#Importing transcript counts from file
 tr_counts <- fread(paste(OUT_PREF, "abundance_tpm.tsv", sep = "_"), header = T)
 tr_counts[, "tpm" := sqrt(tpm)]
 tr_counts[tpm == 0, "tpm" := 1e-100]
  
#Merging read mappings with transcript counts
 read_mappings <- merge(read_mappings, setnames(tr_counts, c("transcript", "tpm")), by = "transcript")
 rm(tr_counts)
 invisible(gc(full = T))
 
#Probabilistic assignment of reads with probabilities proportional to TPM
 read_mappings[, "sample_prob" := tpm / sum(tpm), by = "read_id"]
 set.seed(SEED)
 read_mappings[, "keep_tr" := fifelse(sample(transcript, size = 1, prob = sample_prob) == transcript, T, F), by = "read_id"]
 samp_mapping_row_N <- read_mappings[keep_tr == T, row_N]
 
#Writing row numbers (ids) of the sampled read-transcript combinations to STDOUT
#1 Reads mapped to a single transcript
#2 Sampled read-transcript combinations 
 write.table(data.table(row_N = c(single_mapped_read_row_N, samp_mapping_row_N)), file = "", sep = "\n", row.names = F, col.names = F, quote = F)
} else {
 write.table(data.table(row_N = single_mapped_read_row_N), file = "", sep = "\n", row.names = F, col.names = F, quote = F)
}