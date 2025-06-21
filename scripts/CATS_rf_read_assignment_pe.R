#!/usr/bin/env Rscript
#Read assignment script (paired-end mode)
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
EDIT_DIST_FACT <- as.numeric(ext_args[1])
SEED <- as.numeric(ext_args[2])
OUT_PREF <- ext_args[3]

#Importing read mappings from STDIN
read_mappings <- fread('file:///dev/stdin', header = F)
setnames(read_mappings, c("read_id", "transcript", "pair_N", "frag_id", "row_N", "ed_dist"))

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
 read_mappings[, "read_id" := NULL]
 
#Saving sampled reads with an unmapped pair and continuing with mapped pairs
 read_mappings[, "frag_total_N" := length(unique(pair_N)), by = "frag_id"]
 unmapped_pair_samp_map_row_N <- read_mappings[keep_tr == T & frag_total_N == 1, row_N]
 read_mappings <- read_mappings[frag_total_N == 2]
 read_mappings[, "frag_total_N" := NULL]
 
#Merging pairs mapped to the same transcript
 if(read_mappings[, .N] != 0) {
  setkey(read_mappings, "transcript", "frag_id")
  read_mappings <- merge(read_mappings[pair_N == 1], read_mappings[pair_N == 2], all = T, by = c("transcript", "frag_id"))
  read_mappings[, c("transcript", "tpm.y") := NULL]
  
#Turning transcript counts into ranks to avoid identical maximum TPM values
  set.seed(SEED)
  read_mappings[, "tpm.x" := rank(tpm.x, ties.method = "first")]
    
#Determining pairs mapped to the same transcript with one read sampled on this transcript
  read_mappings[, "diff_tr_sampled_pair" := fifelse((keep_tr.x == T & keep_tr.y == F) | (keep_tr.x == F & keep_tr.y == T), T, F)]
    
#Choosing the most abundant of these transcripts
  read_mappings[diff_tr_sampled_pair == T, "diff_tr_sampled_pair" := fifelse(tpm.x == max(tpm.x), T, F), by = "frag_id"]
  
#Saving these pairs
  diff_tr_sampled_frag <- read_mappings[diff_tr_sampled_pair == T, frag_id]

#Saving pairs sampled on the same transcript
  sane_tr_sampled_frag <- read_mappings[keep_tr.x == T & keep_tr.y == T, frag_id]
    
#Determining pairs mapped to the same transcript with both reads sampled on transcripts to which their pair is not mapped
  read_mappings[frag_id %in% c(diff_tr_sampled_frag, sane_tr_sampled_frag) == F & keep_tr.x == F & keep_tr.y == F, "same_tr_unsampled_pair" := T]
    
#Choosing the most abundant of these transcripts
  read_mappings[same_tr_unsampled_pair == T, "same_tr_unsampled_pair" := fifelse(tpm.x == max(tpm.x), T, F), by = "frag_id"]
  
#Determining pairs mapped to different transcripts, without common transcripts
  sane_tr_unsampled_frag <- read_mappings[same_tr_unsampled_pair == T, frag_id]
  read_mappings[frag_id %in% c(diff_tr_sampled_frag, sane_tr_sampled_frag, sane_tr_unsampled_frag)  == F, "wo_same_tr_sampled_pair" := T]
    
#Writing row numbers (ids) of the sampled read-transcript combinations to STDOUT
#1 Reads mapped to a single transcript
#2 Sampled reads with an unmapped pair
#3 Pairs sampled on the same transcript
#4 Pairs mapped to the same transcript with one read sampled on this transcript - the transcript with the highest abundance is chosen for both reads if more such transcripts exist
#5 Pairs mapped to the same transcript with both reads sampled on transcripts on which their pair is not sampled - the transcript with the highest abundance is chosen for both reads if more such transcripts exist
#6 Sampled pairs mapped to different transcripts, without common transcripts
  write.table(data.table(row_N = c(single_mapped_read_row_N, unmapped_pair_samp_map_row_N, read_mappings[keep_tr.x == T & keep_tr.y == T, c(row_N.x, row_N.y)], read_mappings[diff_tr_sampled_pair == T, c(row_N.x, row_N.y)], read_mappings[same_tr_unsampled_pair == T, c(row_N.x, row_N.y)], read_mappings[wo_same_tr_sampled_pair == T & keep_tr.x == T, row_N.x], read_mappings[wo_same_tr_sampled_pair == T & keep_tr.y == T, row_N.y])), file = "", sep = "\n", row.names = F, col.names = F, quote = F)
 } else {
  write.table(data.table(row_N = c(single_mapped_read_row_N, unmapped_pair_samp_map_row_N)), file = "", sep = "\n", row.names = F, col.names = F, quote = F)
 }
} else {
 write.table(data.table(row_N = single_mapped_read_row_N), file = "", sep = "\n", row.names = F, col.names = F, quote = F)
}