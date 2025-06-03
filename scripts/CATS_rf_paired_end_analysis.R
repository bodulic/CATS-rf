#!/usr/bin/env Rscript
#Paired-end analysis script
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
TR_LEN_DT_PATH <- ext_args[2]
INCOMP_BRIDGE_TR_END_DIST <- as.numeric(ext_args[3])
STRANDNESS <- ext_args[4]
DIST_THR_CORR_FACT <- as.numeric(ext_args[5])
DIST_THR_FACT_LOWER <- as.numeric(ext_args[6])
DIST_THR_FACT_HIGHER <- as.numeric(ext_args[7])
OUT_PREF <- ext_args[8]
IMP_WITHIN_READ_PROP_BREAKS <- unlist(strsplit(ext_args[9], ","))
DIFF_TR_PAIR_READ_PROP_BREAKS <- unlist(strsplit(ext_args[10], ","))
ALPHA_COMP_FACT <- as.numeric(ext_args[11])
BETA_COMP_FACT <- as.numeric(ext_args[12])
FRAG_TR_BRIDGE_N <- as.numeric(ext_args[13])

#Cutting a categorical variable into intervals (function)
get_category_intervals <- function(input_dt, variable_to_cut, variable_breaks, cut_variable_name) {
 input_dt[get(variable_to_cut) < variable_breaks[1], (cut_variable_name) := paste0("<", variable_breaks[1])]
 input_dt[get(variable_to_cut) >= variable_breaks[1], (cut_variable_name) := cut(get(variable_to_cut), breaks = variable_breaks, include.lowest = T)]
 input_dt[is.na(get(cut_variable_name)), (cut_variable_name) := paste0(">", variable_breaks[length(variable_breaks)])]
}

#Calculating bridging event penalty components (function)
get_bridge_penalty_comp <- function(tr_len, tr_start, tr_end, bridge_penalty_comp) {
 imp_map_pair_diff_tr[, "tr_len_med" := (get(tr_len) + 1) / 2]
 imp_map_pair_diff_tr[, "read_len_med" := (get(tr_start) + get(tr_end)) / 2]
 imp_map_pair_diff_tr[read_len_med <= tr_len_med, (bridge_penalty_comp) := 1 - (get(tr_start) - 1) / (floor(0.5 * get(tr_len)) - ceiling((get(tr_end) - get(tr_start) + 1) / 2))]
 imp_map_pair_diff_tr[read_len_med > tr_len_med, (bridge_penalty_comp) := 1 - (get(tr_len) - get(tr_end)) / (floor(0.5 * get(tr_len)) - ceiling((get(tr_end) - get(tr_start) + 1) / 2))]
 imp_map_pair_diff_tr[read_len_med <= tr_len_med & get(tr_len) %% 2 == 1 & (get(tr_end) - get(tr_start) + 1) %% 2 == 1, (bridge_penalty_comp) := 1 - (get(tr_start) - 1) / (floor(0.5 * get(tr_len)) - floor((get(tr_end) - get(tr_start) + 1) / 2))]
 imp_map_pair_diff_tr[read_len_med > tr_len_med & get(tr_len) %% 2 == 1 & (get(tr_end) - get(tr_start) + 1) %% 2 == 1, (bridge_penalty_comp) := 1 - (get(tr_len) - get(tr_end)) / (floor(0.5 * get(tr_len)) - floor((get(tr_end) - get(tr_start) + 1) / 2))]
 imp_map_pair_diff_tr[, c("tr_len_med", "read_len_med") := NULL]
}

#Importing transcript lengths from file
setDTthreads(THREAD_N)
tr_lengths <- fread(TR_LEN_DT_PATH, header = F)
setnames(tr_lengths, c("transcript", "tr_len"))

#Defining empty variables
mapped_read_N_vec <- c()
end_mapped_read_N_vec <- c()
imp_map_pair_strand_ori_dt <- list()
pair_dist_q3 <- NA
imp_map_pair_dist_dt <- list()
local_fidelity_stats_dt <- list()
imp_map_pair_diff_tr_dt <- list()
integrity_stats_dt <- list()

#Defining read mapping table import pattern
files_in_directory <- grep("x[0-9]+$", list.files(path = "."), value = T)

#Looping over split files
for (i in 1 : length(files_in_directory)) {
#Importing split read mapping tables from files
 read_mappings <- fread(files_in_directory[i], header = F)
 setnames(read_mappings, c("transcript", "tr_start", "tr_end", "tr_strand", "pair_N", "frag_id"))
  
#Adding transcript lengths to read mapping table
 read_mappings <- merge(read_mappings, tr_lengths, by = "transcript")
  
#Calculating the number of mapped reads per transcript
 tr_mapped_read_N <- read_mappings[, .(read_N = .N), by = "transcript"]
  
#Calculating the number of mapped reads from the same fragment
 read_mappings[, "frag_total_N" := .N, by = "frag_id"]
  
#Calculating the number of reads with pair not mapped to the assembly per transcript
 unmapped_pair_read_N <- read_mappings[frag_total_N == 1, .(unmap_pair_read_N = .N), by = "transcript"]
 if (unmapped_pair_read_N[, .N] != 0) {
  local_fidelity_stats <- merge(tr_mapped_read_N, unmapped_pair_read_N, by = "transcript", all.x = T)
 } else {
  local_fidelity_stats <- tr_mapped_read_N
  local_fidelity_stats[, "unmap_pair_read_N" := 0]
 }
 rm(unmapped_pair_read_N)
  
#Calculating the number of reads with pair not mapped to the assembly on transcript ends (reads suggesting transcript end incompleteness or fragmentation)
 unmapped_pair_tr_end_read_N <- read_mappings[frag_total_N == 1 & (tr_start <= INCOMP_BRIDGE_TR_END_DIST | tr_end >= (tr_len - INCOMP_BRIDGE_TR_END_DIST +  1)), .(unmap_pair_end_read_N = .N), by = "transcript"]
 if (unmapped_pair_tr_end_read_N[, .N] != 0) {
  local_fidelity_stats <- merge(local_fidelity_stats, unmapped_pair_tr_end_read_N, by = "transcript", all.x = T)
 } else {
  local_fidelity_stats[, "unmap_pair_end_read_N" := 0]
 }
 rm(unmapped_pair_tr_end_read_N)
  
#Calculating the number of reads mapped to transcript ends per transcript
 tr_first_end_read_N <-  read_mappings[tr_start <= INCOMP_BRIDGE_TR_END_DIST, .N, by = "transcript"]
 tr_second_end_read_N <- read_mappings[tr_end >= tr_len - INCOMP_BRIDGE_TR_END_DIST + 1, .N, by = "transcript"]
 if (tr_first_end_read_N[, .N] != 0 & tr_second_end_read_N[, .N] != 0) {
  tr_end_read_N <- merge(tr_first_end_read_N, tr_second_end_read_N, by = "transcript", all = T)
  tr_end_read_N[is.na(N.x), "N.x" := 0]
  tr_end_read_N[is.na(N.y), "N.y" := 0]                      
 } else if (tr_first_end_read_N[, .N] == 0 & tr_second_end_read_N[, .N] != 0 ) {
  tr_end_read_N <- setnames(tr_second_end_read_N, c("transcript", "N.y"))
  tr_end_read_N[, "N.x" := 0]
 } else if (tr_first_end_read_N[, .N] != 0 & tr_second_end_read_N[, .N] == 0 ) {
  tr_end_read_N <- setnames(tr_first_end_read_N, c("transcript", "N.x"))
  tr_end_read_N[, "N.y" := 0]
 } else {
  tr_end_read_N <- data.table(transcript = local_fidelity_stats[, transcript], N.x = 0, N.y = 0)
 }
 rm(tr_first_end_read_N, tr_second_end_read_N)
 tr_end_read_N[, "end_mapped_read_N" := N.x + N.y]
 tr_end_read_N[, c("N.x", "N.y") := NULL]
  
#Calculating the number of mapped reads 
 mapped_read_N_vec[i] <- read_mappings[, .N]
  
#Calculating the number of reads mapped to transcript ends 
 end_mapped_read_N_vec[i] <- tr_end_read_N[, sum(end_mapped_read_N)]
  
#Continuing with mapped pairs
 read_mappings <- read_mappings[frag_total_N == 2]
  
#Preparing for the analysis of mapped pairs
 if (read_mappings[, .N] != 0) {
  read_mappings[, "frag_total_N" := NULL]
  setkey(read_mappings, "frag_id")
  read_mappings <- merge(read_mappings[pair_N == 1], read_mappings[pair_N == 2], by = "frag_id")
  read_mappings[, c("pair_N.x", "pair_N.y", "frag_id") := NULL]
    
#Identifying pairs which map to different transcripts
  imp_map_pair_diff_tr <- read_mappings[transcript.x != transcript.y]
  read_mappings <- read_mappings[transcript.x == transcript.y]
    
#Identifying pairs which map to an inappropriate strand or in an unexpected orientation
  if (read_mappings[, .N] != 0) {
   if (STRANDNESS == "fr") {
    read_mappings[, "strand_ori" := fifelse(tr_strand.x == "+" & tr_strand.y == "-" & tr_start.x <= tr_start.y, "correct", "incorrect")]
   } else if (STRANDNESS == "rf") {
    read_mappings[, "strand_ori" := fifelse(tr_strand.x == "-" & tr_strand.y == "+" & tr_start.y <= tr_start.x, "correct", "incorrect")]
   } else {
    read_mappings[, "strand_ori" := fifelse(tr_strand.x != tr_strand.y, "correct", "incorrect")]
   }
   imp_map_pair_strand_ori <- read_mappings[strand_ori == "incorrect"]
   imp_map_pair_strand_ori[, "strand_ori" := NULL]
   read_mappings <- read_mappings[strand_ori == "correct"]
   read_mappings[, "strand_ori" := NULL]
      
   if (imp_map_pair_strand_ori[, .N] != 0) {
    local_fidelity_stats <- merge(local_fidelity_stats, setnames(imp_map_pair_strand_ori[, 2 * .N, by = "transcript.x"], c("transcript", "imp_strand_ori_read_N")), by = "transcript", all.x = T)
    imp_map_pair_strand_ori_dt[[i]] <- imp_map_pair_strand_ori
   } else {
    local_fidelity_stats[, "imp_strand_ori_read_N" := 0]
   }
   rm(imp_map_pair_strand_ori)
      
#Calculating distance between pairs
   if (read_mappings[, .N] != 0) {
    if (STRANDNESS == "fr") {
     read_mappings[, "pair_distance" := tr_start.y - tr_end.x]
    } else if (STRANDNESS == "rf") {
     read_mappings[, "pair_distance" := tr_start.x - tr_end.y]
    } else {
     read_mappings[, "pair_distance" := fifelse(tr_start.x < tr_start.y, tr_start.y - tr_end.x, tr_start.x - tr_end.y)]
    }
    read_mappings[pair_distance < 0, "pair_distance" := 0]
    read_mappings[, "pair_distance" := abs(pair_distance)]   
        
#Analysing pair distance distribution from a subset of fragments
    if (is.na(pair_dist_q3))  {
     pair_dist_q3 <- quantile(read_mappings[, pair_distance], 0.75) 
     pair_dist_iqr_corr <- pair_dist_q3 - quantile(read_mappings[, pair_distance], 0.25) + DIST_THR_CORR_FACT
    }
        
#Identifying pairs which map too far apart and calculating distance penalty
    imp_map_pair_dist <- read_mappings[pair_distance >= pair_dist_q3 + DIST_THR_FACT_LOWER * pair_dist_iqr_corr]
    rm(read_mappings)
    if (imp_map_pair_dist[, .N] != 0) {
     imp_map_pair_dist[, "distance_penalty" := pair_distance / (pair_dist_q3 + DIST_THR_FACT_HIGHER * pair_dist_iqr_corr)]
     imp_map_pair_dist[distance_penalty > 1, "distance_penalty" := 1]
     local_fidelity_stats <- merge(local_fidelity_stats, setnames(imp_map_pair_dist[, 2 * .N, by = "transcript.x"], c("transcript", "imp_dist_read_N")), by = "transcript", all.x = T)
     local_fidelity_stats <- merge(local_fidelity_stats, setnames(imp_map_pair_dist[, 2 * sum(distance_penalty), by = "transcript.x"], c("transcript", "distance_pen_sum")), by = "transcript", all.x = T)
     imp_map_pair_dist_dt[[i]] <- imp_map_pair_dist
    } else {
     local_fidelity_stats[, ':=' ("imp_dist_read_N" = 0, "distance_pen_sum" = 0)]
    }
    rm(imp_map_pair_dist)
   } else {
    rm(read_mappings)
    local_fidelity_stats[, ':=' ("imp_dist_read_N" = 0, "distance_pen_sum" = 0)]
   }
  } else {
   rm(read_mappings)
   local_fidelity_stats[, ':=' ("imp_strand_ori_read_N" = 0, "imp_dist_read_N" = 0, "distance_pen_sum" = 0)]
  }
    
#Adding the number of reads mapped to transcript ends per transcript
  local_fidelity_stats <- merge(local_fidelity_stats, tr_end_read_N, by = "transcript", all.x = T)
    
#Saving local fidelity statistics
  local_fidelity_stats[is.na(local_fidelity_stats)] <- 0
  local_fidelity_stats_dt[[i]] <- local_fidelity_stats
  rm(local_fidelity_stats)
    
#Saving the number of mapped reads and the number of reads mapped to transcript ends per transcript
  integrity_stats <- merge(tr_mapped_read_N, tr_end_read_N, by = "transcript", all.x = T)
  rm(tr_mapped_read_N, tr_end_read_N)
    
#Calculating the per-transcript number of reads with pairs mapped to another transcript
  if (imp_map_pair_diff_tr[, .N] != 0) {  
   imp_map_pair_diff_tr_N <- rbindlist(list(imp_map_pair_diff_tr[, .(transcript.x)], setnames(imp_map_pair_diff_tr[, .(transcript.y)], "transcript.x")))[, .(diff_tr_pair_read_N = .N), by = "transcript.x"]
   integrity_stats <- merge(integrity_stats, setnames(imp_map_pair_diff_tr_N, c("transcript", "diff_tr_pair_read_N")), by ="transcript", all.x = T)
   rm(imp_map_pair_diff_tr_N)
      
#Identifying bridging events
   imp_map_pair_diff_tr[, "bridging_event" := fifelse((tr_start.x <= INCOMP_BRIDGE_TR_END_DIST & tr_start.y <= INCOMP_BRIDGE_TR_END_DIST) | (tr_start.x <= INCOMP_BRIDGE_TR_END_DIST & tr_end.y >= tr_len.y - INCOMP_BRIDGE_TR_END_DIST + 1) | (tr_start.y <= INCOMP_BRIDGE_TR_END_DIST & tr_end.x >= tr_len.x - INCOMP_BRIDGE_TR_END_DIST + 1) | (tr_end.x >= tr_len.x - INCOMP_BRIDGE_TR_END_DIST + 1 & tr_end.y >= tr_len.y - INCOMP_BRIDGE_TR_END_DIST + 1), T, F)]
   tr_bridge_N <- rbindlist(list(imp_map_pair_diff_tr[bridging_event == T, .(transcript.x)], setnames(imp_map_pair_diff_tr[bridging_event == T, .(transcript.y)], "transcript.x")))[, .(bridge_N = .N), by = "transcript.x"]
   if (tr_bridge_N[, .N] != 0) {
    integrity_stats <- merge(integrity_stats, setnames(tr_bridge_N, c("transcript", "bridge_N")), by = "transcript", all.x = T)
   } else {
    integrity_stats[, "bridge_N" := 0]
   }
   rm(tr_bridge_N)
      
#Calculating bridge penalty
   get_bridge_penalty_comp("tr_len.x", "tr_start.x", "tr_end.x", "bridge_penalty_comp.x")
   get_bridge_penalty_comp("tr_len.y", "tr_start.y", "tr_end.y", "bridge_penalty_comp.y")
   imp_map_pair_diff_tr[, "bridge_penalty" := (bridge_penalty_comp.x + bridge_penalty_comp.y) / 2]
   imp_map_pair_diff_tr[, c("bridge_penalty_comp.x", "bridge_penalty_comp.y") := NULL]
   bridge_penalty_sum <- rbindlist(list(imp_map_pair_diff_tr[, .(transcript.x, bridge_penalty)], setnames(imp_map_pair_diff_tr[, .(transcript.y, bridge_penalty)], c("transcript.x", "bridge_penalty"))))[, .(bridge_pen_sum = sum(bridge_penalty)), by = "transcript.x"]
   integrity_stats <- merge(integrity_stats, setnames(bridge_penalty_sum, c("transcript", "bridge_pen_sum")), by = "transcript", all.x = T)
   rm(bridge_penalty_sum)
   imp_map_pair_diff_tr_dt[[i]] <- imp_map_pair_diff_tr
   rm(imp_map_pair_diff_tr)
  } else {
   integrity_stats[, ':=' ("diff_tr_pair_read_N" = 0, "bridge_N" = 0, "bridge_pen_sum" = 0)]
  }
    
#Saving integrity statistics
  integrity_stats[is.na(integrity_stats)] <- 0
  integrity_stats_dt[[i]] <- integrity_stats
  rm(integrity_stats)
    
 } else {
  local_fidelity_stats[, ':=' ("imp_strand_ori_read_N" = 0, "imp_dist_read_N" = 0, "distance_pen_sum" = 0)]
  local_fidelity_stats <- merge(local_fidelity_stats, tr_end_read_N, by = "transcript", all.x = T)
  local_fidelity_stats[is.na(local_fidelity_stats)] <- 0
  local_fidelity_stats_dt[[i]] <- local_fidelity_stats
  rm(local_fidelity_stats)
  integrity_stats <- merge(tr_mapped_read_N, tr_end_read_N, by = "transcript", all.x = T)
  rm(tr_mapped_read_N, tr_end_read_N)
  integrity_stats[is.na(end_mapped_read_N), "end_mapped_read_N" := 0]
  integrity_stats[, ':=' ("diff_tr_pair_read_N" = 0, "bridge_N" = 0, "bridge_pen_sum" = 0)]
  integrity_stats_dt[[i]] <- integrity_stats
  rm(integrity_stats)
 }
}

#Saving all transcripts
tr_vec <- tr_lengths[, transcript]
rm(tr_lengths)

#Summing the number of mapped reads and the number of reads mapped to transcript ends
mapped_read_N_sum <- sum(mapped_read_N_vec)
end_mapped_read_N_sum <- sum(end_mapped_read_N_vec)

#Writing pairs which map to an inappropriate strand or in an unexpected orientation to file
imp_map_pair_strand_ori_dt <- rbindlist(imp_map_pair_strand_ori_dt)
if (imp_map_pair_strand_ori_dt[, .N] != 0) {
 setorder(imp_map_pair_strand_ori_dt, transcript.x)
 write.table(imp_map_pair_strand_ori_dt, file = paste(OUT_PREF, "read_pairs_mapping_to_inconsistent_strand_or_orientation.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
}
rm(imp_map_pair_strand_ori_dt)

#Writing pairs mapping too far apart to file
imp_map_pair_dist_dt <- rbindlist(imp_map_pair_dist_dt)
if (imp_map_pair_dist_dt[, .N] != 0) {
 setorder(imp_map_pair_dist_dt, transcript.x)
 write.table(imp_map_pair_dist_dt, file = paste(OUT_PREF, "read_pairs_mapping_too_far_apart.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
}
rm(imp_map_pair_dist_dt)

#Writing pairs mapping to different transcripts to file
imp_map_pair_diff_tr_dt <- rbindlist(imp_map_pair_diff_tr_dt)
if (imp_map_pair_diff_tr_dt[, .N] != 0) {
 setorder(imp_map_pair_diff_tr_dt, transcript.x)
 write.table(imp_map_pair_diff_tr_dt, file = paste(OUT_PREF, "read_pairs_mapping_to_different_transcripts.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
}
rm(imp_map_pair_diff_tr_dt)

#Processing local fidelity statistics
local_fidelity_stats_dt <- rbindlist(local_fidelity_stats_dt)
local_fidelity_stats_dt <- local_fidelity_stats_dt[, lapply(.SD, sum), by = "transcript", .SDcols = -"transcript"]

#Saving fully uncovered transcripts
uncov_tr_vec <- tr_vec[tr_vec %in% local_fidelity_stats_dt[, transcript] == F]

#Calculating the number of transcripts
tr_N <- length(tr_vec)
rm(tr_vec)

#Calculating the proportion of reads with unmapped pairs per transcript
local_fidelity_stats_dt[, "unmap_pair_read_prop" := unmap_pair_read_N / read_N]

#Calculating the per-transcript proportion of reads with unmapped pairs on transcript ends
local_fidelity_stats_dt[end_mapped_read_N > 0, "unmap_pair_end_read_prop" := unmap_pair_end_read_N / end_mapped_read_N]
local_fidelity_stats_dt[, "end_mapped_read_N" := NULL]

#Calculating the proportion of reads with pairs which map to an inappropriate strand or in an unexpected orientation per transcript
local_fidelity_stats_dt[, "imp_strand_ori_read_prop" := imp_strand_ori_read_N / read_N]

#Normalizing the total distance penalty per transcript
local_fidelity_stats_dt[, "norm_distance_pen" := distance_pen_sum / read_N]
local_fidelity_stats_dt[, "distance_pen_sum" := NULL]

#Calculating the number and proportion of improperly paired reads within a transcript
local_fidelity_stats_dt[, "imp_within_read_N" := unmap_pair_read_N + imp_strand_ori_read_N + imp_dist_read_N]
local_fidelity_stats_dt[, "imp_within_read_prop" := imp_within_read_N / read_N]
local_fidelity_stats_dt[, "read_N" := NULL]
get_category_intervals(local_fidelity_stats_dt, "imp_within_read_prop", IMP_WITHIN_READ_PROP_BREAKS, "imp_within_read_prop_cat")

#Calculating local fidelity score component (Sl)
local_fidelity_stats_dt[, "loc_fid_score_comp" := 1 - sqrt(unmap_pair_read_prop + imp_strand_ori_read_prop + norm_distance_pen)]

#Reordering the local fidelity statistics table
setcolorder(local_fidelity_stats_dt, c("transcript", "unmap_pair_read_N", "unmap_pair_read_prop", "unmap_pair_end_read_N", "unmap_pair_end_read_prop", "imp_strand_ori_read_N", "imp_strand_ori_read_prop", "imp_dist_read_N", "norm_distance_pen", "imp_within_read_N", "imp_within_read_prop", "imp_within_read_prop_cat", "loc_fid_score_comp"))

#Calculating the number of reads with unmapped pairs
unmap_pair_read_N_sum <- local_fidelity_stats_dt[, sum(unmap_pair_read_N)]

#Calculating the number of reads with unmapped pairs on transcript ends
unmap_pair_end_read_N_sum <- local_fidelity_stats_dt[, sum(unmap_pair_end_read_N)]

#Calculating the number of reads whose pairs map to an inappropriate strand or in an unexpected orientation
imp_strand_ori_read_N_sum <- local_fidelity_stats_dt[, sum(imp_strand_ori_read_N)] 

#Calculating the number of reads whose pairs map too far apart
imp_dist_read_N_sum <- local_fidelity_stats_dt[, sum(imp_dist_read_N)]

#Calculating the number of improperly paired reads within a transcript
imp_within_read_N_sum <- unmap_pair_read_N_sum + imp_strand_ori_read_N_sum + imp_dist_read_N_sum

#Calculating summaries of local fidelity metrics
unmap_pair_end_read_prop_summary <- summary(local_fidelity_stats_dt[, unmap_pair_end_read_prop])
imp_within_read_prop_summary <- summary(local_fidelity_stats_dt[, imp_within_read_prop])
loc_fid_score_comp_summary <- summary(local_fidelity_stats_dt[, loc_fid_score_comp])

#Adding uncovered transcripts
if(length(uncov_tr_vec) > 0) {
 local_fidelity_stats_dt <- rbindlist(list(local_fidelity_stats_dt, data.table(transcript = uncov_tr_vec, unmap_pair_read_N = 0, unmap_pair_read_prop = NA, unmap_pair_end_read_N = 0, unmap_pair_end_read_prop = NA, imp_strand_ori_read_N = 0, imp_strand_ori_read_prop = NA, imp_dist_read_N = 0, norm_distance_pen = NA, imp_within_read_N = 0, imp_within_read_prop = NA, imp_within_read_prop_cat = NA, loc_fid_score_comp = NA)))
}

#Writing local fidelity statistics to file
setnames(local_fidelity_stats_dt, c("transcript", "unmapped_pair_read_N", "unmapped_pair_read_prop", "unmapped_pair_tr_end_read_N" ,"unmapped_pair_tr_end_read_prop", "improp_pair_strand_orientation_read_N", "improp_pair_strand_orientation_read_prop", "improp_pair_distance_read_N", "norm_distance_penalty", "improp_pair_within_tr_read_N", "improp_pair_within_tr_read_prop", "improp_pair_within_tr_read_prop_category", "local_fidelity_score_component"))
setorder(local_fidelity_stats_dt, transcript)
write.table(local_fidelity_stats_dt, file = paste(OUT_PREF, "local_fidelity_stats.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(local_fidelity_stats_dt)

#Processing integrity statistics
integrity_stats_dt <- rbindlist(integrity_stats_dt)
integrity_stats_dt <- integrity_stats_dt[, lapply(.SD, sum), by = "transcript", .SDcols = -"transcript"]

#Calculating the proportion of reads whose pairs map to other transcripts per transcript
integrity_stats_dt[, "diff_tr_pair_read_prop" := diff_tr_pair_read_N / read_N]
get_category_intervals(integrity_stats_dt, "diff_tr_pair_read_prop", DIFF_TR_PAIR_READ_PROP_BREAKS, "diff_tr_pair_read_prop_cat")

#Calculating the proportion of reads representing bridging events on transcript ends
integrity_stats_dt[end_mapped_read_N > 0, "bridge_prop" := bridge_N / end_mapped_read_N]
integrity_stats_dt[, "end_mapped_read_N" := NULL]

#Calculating integrity score component (Si)
integrity_stats_dt[, "bridge_idx" := 1 - (bridge_pen_sum / read_N)]
integrity_stats_dt[, "integ_score_comp" := bridge_idx^ALPHA_COMP_FACT / (bridge_idx^ALPHA_COMP_FACT + (1 - bridge_idx)^BETA_COMP_FACT)]
integrity_stats_dt[, c("read_N", "bridge_pen_sum", "bridge_idx") := NULL]

#Reordering the integrity statistics table
setcolorder(integrity_stats_dt, c("transcript", "diff_tr_pair_read_N", "diff_tr_pair_read_prop", "diff_tr_pair_read_prop_cat", "bridge_N", "bridge_prop", "integ_score_comp"))

#Calculating the number of reads whose pairs map to other transcripts
diff_tr_pair_read_N_sum <- integrity_stats_dt[, sum(diff_tr_pair_read_N)]

#Calculating the number of fragmented transcripts
frag_tr_N <- integrity_stats_dt[bridge_N >= FRAG_TR_BRIDGE_N, .N]

#Calculating the number of bridging events
bridge_N_sum <- integrity_stats_dt[, sum(bridge_N)]

#Calculating summaries of integrity metrics
diff_tr_pair_read_prop_summary <- summary(integrity_stats_dt[, diff_tr_pair_read_prop])
bridge_prop_summary <- summary(integrity_stats_dt[, bridge_prop])
integ_score_comp_summary <- summary(integrity_stats_dt[, integ_score_comp])

#Adding uncovered transcripts
if(length(uncov_tr_vec) > 0) {
 integrity_stats_dt <- rbindlist(list(integrity_stats_dt, data.table(transcript = uncov_tr_vec, diff_tr_pair_read_N = 0, diff_tr_pair_read_prop = NA, diff_tr_pair_read_prop_cat = NA, bridge_N = 0, bridge_prop = NA, integ_score_comp = NA)))
}

#Writing integrity statistics to file
setnames(integrity_stats_dt, c("transcript", "pair_mapped_to_other_tr_N", "pair_mapped_to_other_tr_prop", "pair_mapped_to_other_tr_prop_category", "bridge_N", "bridge_prop", "integrity_score_component"))
setorder(integrity_stats_dt, transcript)
write.table(integrity_stats_dt, file = paste(OUT_PREF, "integrity_stats.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(integrity_stats_dt)

#Writing paired-end analysis summary to file
write.table(data.table(parameter = c("N, % reads with unmapped pairs", "N, % reads with unmapped pairs on transcript ends", "% of reads with unmapped pairs on transcript ends per transcript (mean, IQR)", "N, % reads with pairs mapping to an inappropriate strand or in an unexpected orientation", "N, % reads mapping too far from their pair", "N, % improperly paired reads within a transcript", "% of improperly paired reads within a transcript per transcript (mean, IQR)", "Local fidelity score component (mean, IQR)", "N, % reads with pairs mapping to different transcripts", "% of reads with pairs mapping to different transcripts per transcript (mean, IQR)", "N, % fragmented transcripts", "N, % reads representing bridging events on transcript ends", "% of reads representing bridging events on transcript ends per transcript (mean, IQR)", "Integrity score component (mean, IQR)"), value = c(paste0(unmap_pair_read_N_sum, ", ", round(100 * unmap_pair_read_N_sum / mapped_read_N_sum, 2), "%"), paste0(unmap_pair_end_read_N_sum, ", ", round(100 * unmap_pair_end_read_N_sum / end_mapped_read_N_sum, 2), "%"), paste0(round(100 * unmap_pair_end_read_prop_summary[4], 2), "%, ", round(100 * unmap_pair_end_read_prop_summary[2], 2), "%-", round(100 * unmap_pair_end_read_prop_summary[5], 2), "%"), paste0(imp_strand_ori_read_N_sum, ", ", round(100 * imp_strand_ori_read_N_sum / mapped_read_N_sum, 2), "%"), paste0(imp_dist_read_N_sum, ", ", round(100 * imp_dist_read_N_sum / mapped_read_N_sum, 2), "%"), paste0(imp_within_read_N_sum, ", ", round(100 * imp_within_read_N_sum / mapped_read_N_sum, 2), "%"), paste0(round(100 * imp_within_read_prop_summary[4], 2), "%, ", round(100 * imp_within_read_prop_summary[2], 2), "%-", round(100 * imp_within_read_prop_summary[5], 2), "%"), paste0(round(loc_fid_score_comp_summary[4], 3), ", ", round(loc_fid_score_comp_summary[2], 3), "-", round(loc_fid_score_comp_summary[5], 3)), paste0(diff_tr_pair_read_N_sum, ", ", round(100 * diff_tr_pair_read_N_sum / mapped_read_N_sum, 2), "%"), paste0(round(100 * diff_tr_pair_read_prop_summary[4], 2), "%, ", round(100 * diff_tr_pair_read_prop_summary[2], 2), "%-", round(100 * diff_tr_pair_read_prop_summary[5], 2), "%"), paste0(frag_tr_N, ", ", round(100 * frag_tr_N / tr_N, 2), "%"), paste0(bridge_N_sum, ", ", round(100 * bridge_N_sum / end_mapped_read_N_sum, 2), "%"), paste0(round(100 * bridge_prop_summary[4], 2), "%, ", round(100 * bridge_prop_summary[2], 2), "%-", round(100 * bridge_prop_summary[5], 2), "%"), paste0(round(integ_score_comp_summary[4], 3), ", ", round(integ_score_comp_summary[2], 3), "-", round(integ_score_comp_summary[5], 3)))), file = paste(OUT_PREF, "paired_end_analysis_summmary.tsv", sep = "_"), sep = "\t", row.names = F, col.names = F, quote = F)