#!/usr/bin/env Rscript
#Coverage and accuracy analysis script
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
OUT_PREF <- ext_args[1]
THREAD_N <- as.numeric(ext_args[2])
COV_PROP_BREAKS <- unlist(strsplit(ext_args[3], ","))
TR_COV_MEAN_BREAKS <- unlist(strsplit(ext_args[4], ","))
BASE_COV_BREAKS <- unlist(strsplit(ext_args[5], ","))
LOCAL_COV_WINDOW_SIZE <- as.numeric(ext_args[6])
LOCAL_COV_THR <- as.numeric(ext_args[7])
LCR_PROP_BREAKS <- unlist(strsplit(ext_args[8], ","))
BASE_COV_WEIGHT <- as.numeric(ext_args[9])
LCR_EX_PEN <- as.numeric(ext_args[10])
POS_REL_COV_PROP <- as.numeric(ext_args[11])
COV_MEAN_END_PROP <- as.numeric(ext_args[12])
POS_ACC_PROP <- as.numeric(ext_args[13])
BASE_ACC_THR <- as.numeric(ext_args[14])
ACC_BASE_PROP_BREAKS <- unlist(strsplit(ext_args[15], ","))
BASE_ACC_BREAKS <- unlist(strsplit(ext_args[16], ","))
LOCAL_ACC_WINDOW_SIZE <- as.numeric(ext_args[17])
LOCAL_ACC_THR <- as.numeric(ext_args[18]) 
LAR_PROP_BREAKS <- unlist(strsplit(ext_args[19], ","))
LAR_EX_PEN <- as.numeric(ext_args[20])

#Cutting a categorical variable into intervals (function)
get_category_intervals <- function(input_dt, variable_to_cut, variable_breaks, cut_variable_name) {
 input_dt[get(variable_to_cut) < variable_breaks[1], (cut_variable_name) := paste0("<", variable_breaks[1])]
 input_dt[get(variable_to_cut) >= variable_breaks[1], (cut_variable_name) := cut(get(variable_to_cut), breaks = variable_breaks, include.lowest = T)]
 input_dt[is.na(get(cut_variable_name)), (cut_variable_name) := paste0(">", variable_breaks[length(variable_breaks)])]
}

#Defining empty variables
fully_cov_tr_N_vec  <- c()
fully_uncov_tr_N_vec <- c()
per_base_cov_category_dt <- list()
per_base_higher_cov_dt <-list()
lcr_dt <- list()
pos_cov_dt <- list()
coverage_stats_dt <- list()
base_N_vec <- c()
cov_base_N_vec <- c()
pos_acc_dt <- list()
per_base_acc_category_dt <- list()
per_base_higher_acc_dt <- list()
lar_dt <- list()
accuracy_stats_dt <- list()

#Defining pysamstats output import pattern
files_in_directory <- grep(paste("x[0-9][0-9]_sorted", OUT_PREF, "pysam_out_with_uncov_bases_and_uncov_tr", sep = "_"), list.files(path = "."), value = T)
setDTthreads(THREAD_N)

#Looping over split files
for (i in 1 : length(files_in_directory)) {
#Importing split pysamsstats outputs from files
 per_base_cov_acc <- fread(files_in_directory[i], header = F)
 setnames(per_base_cov_acc, c("transcript", "position", "coverage", "matches", "insertions"))
 
#Accounting for insertions
 per_base_cov_acc[, "matches" := matches - insertions]
 per_base_cov_acc[, "insertions" := NULL]
 per_base_cov_acc[matches < 0, "matches" := 0]

#Coverage analysis
#Calculating transcript lengths
 tr_lengths <- per_base_cov_acc[, .(tr_len = max(position)), by = "transcript"]
 
#Calculating the number and proportion of covered bases per transcript
 tr_cov_base_N <- per_base_cov_acc[coverage != 0, .(cov_N = .N), by = "transcript"]
 tr_cov_base_N <- merge(tr_cov_base_N, tr_lengths, by = "transcript", all.y = T)
 tr_cov_base_N[is.na(cov_N), "cov_N" := 0]
 tr_cov_base_N[, "cov_prop" := cov_N / tr_len]
 tr_cov_base_N[, "tr_len" := NULL]
 get_category_intervals(tr_cov_base_N, "cov_prop", COV_PROP_BREAKS, "cov_prop_cat")

#Saving fully covered and fully uncovered transcripts
 fully_cov_tr_N_vec[i] <- tr_cov_base_N[cov_prop == 1, .N]
 fully_uncov_tr <- tr_cov_base_N[cov_prop == 0, transcript]
 fully_uncov_tr_N_vec[i] <- length(fully_uncov_tr)

#Calculating mean transcript coverage
 tr_cov_mean <- per_base_cov_acc[, .(cov_mean = mean(coverage)), by = "transcript"]
 get_category_intervals(tr_cov_mean, "cov_mean", TR_COV_MEAN_BREAKS, "cov_mean_cat")
 coverage_stats <- merge(tr_cov_base_N, tr_cov_mean, by = "transcript")
 rm(tr_cov_base_N, tr_cov_mean)

#Calculating per-base coverage distribution
 get_category_intervals(per_base_cov_acc, "coverage", BASE_COV_BREAKS, "coverage_cat")
 per_base_cov_category_dt[[i]] <- per_base_cov_acc[, .(base_N = .N), by = "coverage_cat"]
 per_base_cov_acc[, "coverage_cat" := NULL]

#Calculating the number of bases with >= coverage than coverage breaks
 BASE_COV_BREAKS2 <- BASE_COV_BREAKS[-1]
 per_base_higher_cov <- list()
 for (j in 1 : length(BASE_COV_BREAKS2)) {
  per_base_higher_cov[[j]] <- data.table(cov_threshold = BASE_COV_BREAKS2[j], base_N = per_base_cov_acc[coverage >= BASE_COV_BREAKS2[j], .N])
 }
 per_base_higher_cov_dt[[i]] <- rbindlist(per_base_higher_cov)
 rm(per_base_higher_cov)

#Defining low-coverage regions (LCRs)
 per_base_cov_acc[, "cov_rollmean" := frollmean(coverage, n = LOCAL_COV_WINDOW_SIZE, fill = NA, align = "center", algo = "fast"), by = "transcript"]
 per_base_cov_acc[, "lcr" := fifelse(cov_rollmean <= LOCAL_COV_THR, T, F)]
 per_base_cov_acc[, "cov_rollmean" := NULL]
 per_base_cov_acc[is.na(lcr), "lcr" := fifelse(coverage <= LOCAL_COV_THR, T, F)]

#Calculating the proportion of LCR bases per transcript
 lcr_stats <- per_base_cov_acc[lcr == T, .(lcr_N = .N), by = "transcript"]
 lcr_stats <- merge(lcr_stats, tr_lengths, all.y = T)
 lcr_stats[is.na(lcr_N), "lcr_N" := 0]
 lcr_stats[, "lcr_prop" := lcr_N / tr_len]
 lcr_stats[, "tr_len" := NULL]
 get_category_intervals(lcr_stats, "lcr_prop", LCR_PROP_BREAKS, "lcr_prop_cat")
 coverage_stats <- merge(coverage_stats, lcr_stats, by = "transcript")
 rm(lcr_stats)
 
#Saving the list of LCR coordinates
 if (per_base_cov_acc[lcr == T, .N] != 0) {
  per_base_cov_acc[, "lcr_mark" := rleid(paste0(transcript, lcr))]
  lcr <- per_base_cov_acc[lcr == T, .(lcr_start = min(position), lcr_end = max(position)), by = c("transcript", "lcr_mark")]
  lcr[, "lcr_mark" := NULL]
  lcr[, "lcr_len" := lcr_end - lcr_start + 1]
  lcr_dt[[i]] <- lcr
  rm(lcr)

#Calculating coverage score component (Sc)
  per_base_cov_acc[lcr == T, "cov_fun" := fifelse(coverage == 0, 1, 1 / (BASE_COV_WEIGHT * coverage))]
  cov_score_comp_stats <- per_base_cov_acc[lcr == T, .(lcr_cov_pen = sum(cov_fun) + .N * LCR_EX_PEN), by = c("lcr_mark", "transcript")][, .(cov_pen = sum(lcr_cov_pen)), by = "transcript"]
  per_base_cov_acc[, c("lcr", "lcr_mark", "cov_fun") := NULL]
  cov_score_comp_stats <- merge(cov_score_comp_stats, tr_lengths, by = "transcript", all.y = T)
  cov_score_comp_stats[, "cov_score_comp" := 1 - sqrt(cov_pen / (tr_len * (LCR_EX_PEN + 1)))]
  cov_score_comp_stats[, c("tr_len", "cov_pen") := NULL]
  cov_score_comp_stats[is.na(cov_score_comp), "cov_score_comp" := 1]
  coverage_stats <- merge(coverage_stats, cov_score_comp_stats, by = "transcript")
  rm(cov_score_comp_stats)
 } else {
  per_base_cov_acc[, "lcr" := NULL]
  coverage_stats[, "cov_score_comp" := 1]
 }
 rm(tr_lengths)
 
#Identifying uncovered regions
 if (per_base_cov_acc[coverage == 0, .N] != 0) { 
  per_base_cov_acc[, "ucr_mark" := rleid(paste0(transcript, coverage))]
   
#Calculating the maximum uncovered region length per transcript
  tr_ucr_length_max <- per_base_cov_acc[coverage == 0, .(ucr_len = .N), by = c("ucr_mark", "transcript")][, .(ucr_len_max = max(ucr_len)), by = "transcript"]
  per_base_cov_acc[, "ucr_mark" := NULL]
  coverage_stats <- merge(coverage_stats, tr_ucr_length_max, by = "transcript", all.x = T)
  rm(tr_ucr_length_max)
  coverage_stats[is.na(ucr_len_max), "ucr_len_max" := 0]
 } else {
  coverage_stats[, "ucr_len_max" := 0]
 }
 
#Calculating relative coverage by transcript position
 per_base_cov_acc[, "position_prop" := position / max(position), by = "transcript"]
 per_base_cov_acc[, "position_prop_cat" := cut(position_prop, breaks = seq(0, 1, by = POS_REL_COV_PROP), include.lowest = T)]
 per_base_cov_acc[, "coverage_max" := max(coverage), by = "transcript"]
 pos_cov <- per_base_cov_acc[, .(rel_cov_mean = mean(coverage / coverage_max)), by = c("transcript", "position_prop_cat")]
 per_base_cov_acc[, c("coverage_max", "position_prop_cat") := NULL]
 pos_cov_dt[[i]] <- pos_cov
 rm(pos_cov)
 
#Calculating transcript end coverage
 tr_end_cov_mean <- per_base_cov_acc[position_prop <= COV_MEAN_END_PROP | position_prop >= 1 - COV_MEAN_END_PROP + 1e-10, .(end_cov_mean = mean(coverage)), by = "transcript"]
 coverage_stats <- merge(coverage_stats, tr_end_cov_mean)
 rm(tr_end_cov_mean)
 
#Saving coverage statistics
 coverage_stats_dt[[i]] <- coverage_stats
 rm(coverage_stats)

#Accuracy analysis
#Cutting the transcript length into intervals for positional accuracy analysis
 per_base_cov_acc[, "position_prop_cat" := cut(position_prop, breaks = seq(0, 1, by = POS_ACC_PROP), include.lowest = T)]
 per_base_cov_acc[, "position_prop" := NULL]
 
#Calculating the number of bases
 base_N_vec[i] <- per_base_cov_acc[, .N]
 
#Removing uncovered bases
 per_base_cov_acc <- per_base_cov_acc[coverage != 0]
 
#Calculating the number of covered bases
 cov_base_N_vec[i] <- per_base_cov_acc[, .N]
 
#Calculating accuracy
 per_base_cov_acc[, "accuracy" := matches / coverage]
 per_base_cov_acc[, c("matches", "coverage") := NULL]
 
#Calculating accuracy by transcript position
 pos_acc <- per_base_cov_acc[, .(acc_mean = mean(accuracy)), by = c("transcript", "position_prop_cat")]
 per_base_cov_acc[, "position_prop_cat" := NULL]
 pos_acc_dt[[i]] <- pos_acc
 rm(pos_acc)
 
#Calculating the per-transcript number of bases with accuracy >= than the selected threshold
 tr_accurate_base_N <- per_base_cov_acc[accuracy >= BASE_ACC_THR, .(acc_base_N = .N),  by = "transcript"]
 
#Calculating the number of covered bases per transcript
 tr_cov_len <- per_base_cov_acc[, .(cov_len = .N), by = "transcript"]
 
#Calculating the per-transcript proportion of bases with accuracy >= than the selected threshold
 if (tr_accurate_base_N[, .N] != 0) {
  accuracy_stats <- merge(tr_accurate_base_N, tr_cov_len, by = "transcript", all.y = T)
  accuracy_stats[is.na(acc_base_N), "acc_base_N" := 0]
  accuracy_stats[, "acc_base_prop" := acc_base_N / cov_len]
  accuracy_stats[, "cov_len" := NULL]
 } else {
  accuracy_stats <- data.table(transcript = tr_cov_len[, transcript], acc_base_N = 0, acc_base_prop = 0)
 }
 rm(tr_accurate_base_N)
 get_category_intervals(accuracy_stats, "acc_base_prop", ACC_BASE_PROP_BREAKS, "acc_base_prop_cat")
 
 if(length(fully_uncov_tr) != 0) {
  accuracy_stats <- rbindlist(list(accuracy_stats, data.table(transcript = fully_uncov_tr, acc_base_N = 0, acc_base_prop = NA, acc_base_prop_cat = NA)))
 }

#Calculating per-base accuracy distribution
 get_category_intervals(per_base_cov_acc, "accuracy", BASE_ACC_BREAKS, "accuracy_cat")
 per_base_acc_category_dt[[i]] <- per_base_cov_acc[, .(base_N = .N), by = "accuracy_cat"]
 per_base_cov_acc[, "accuracy_cat" := NULL]

#Calculating the number of bases with >= accuracy than accuracy breaks
 BASE_ACC_BREAKS2 <- BASE_ACC_BREAKS[-1]
 per_base_higher_acc <- list()
 for (j in 1 : length(BASE_ACC_BREAKS2)) {
  per_base_higher_acc[[j]] <- data.table(acc_threshold = BASE_ACC_BREAKS2[j], base_N = per_base_cov_acc[accuracy >= BASE_ACC_BREAKS2[j], .N])
 }
 per_base_higher_acc_dt[[i]] <- rbindlist(per_base_higher_acc)
 rm(per_base_higher_acc)

#Defining low-accuracy regions (LARs)
 per_base_cov_acc[, "acc_rollmean" := frollmean(accuracy, n = LOCAL_ACC_WINDOW_SIZE, fill =  NA, align = "center", algo = "fast"), by = "transcript"]
 per_base_cov_acc[, "lar" := fifelse(acc_rollmean <= LOCAL_ACC_THR, T, F)]
 per_base_cov_acc[, "acc_rollmean" := NULL]
 per_base_cov_acc[is.na(lar), "lar" := fifelse(accuracy <= LOCAL_ACC_THR, T, F)]

#Removing fully accurate bases on LAR ends
 if (per_base_cov_acc[lar == T, .N] != 0) {
  per_base_cov_acc[, "lar_mark" := rleid(paste0(transcript, lar))]
  per_base_cov_acc[lar == T, "acc_mark" := rleid(accuracy), by = "lar_mark"]
  per_base_cov_acc[lar == T, "lar" := fifelse(accuracy == 1 & (acc_mark == min(acc_mark) | acc_mark == max(acc_mark)), F, T), by = "lar_mark"]
  per_base_cov_acc[, "acc_mark" := NULL]
 }

#Calculating the proportion of LAR bases per transcript
 lar_stats <- per_base_cov_acc[lar == T, .(lar_N = .N), by = "transcript"]
 lar_stats <- merge(lar_stats, tr_cov_len, all.y = T)
 lar_stats[is.na(lar_N), "lar_N" := 0]
 lar_stats[, "lar_prop" := lar_N / cov_len]
 lar_stats[, "cov_len" := NULL]
 get_category_intervals(lar_stats, "lar_prop", LAR_PROP_BREAKS, "lar_prop_cat")
 
 if(length(fully_uncov_tr) != 0) {
  lar_stats <- rbindlist(list(lar_stats, data.table(transcript = fully_uncov_tr, lar_N = 0, lar_prop = NA, lar_prop_cat = NA)))
 }
 accuracy_stats <- merge(accuracy_stats, lar_stats, by = "transcript")
 rm(lar_stats)

#Saving the list of LAR coordinates
 if (per_base_cov_acc[lar == T, .N] != 0) {
  lar <- per_base_cov_acc[lar == T, .(lar_start = min(position), lar_end = max(position)), by = c("transcript", "lar_mark")]
  lar[, "lar_mark" := NULL]
  lar[, "lar_len" := lar_end - lar_start + 1]
  lar_dt[[i]] <- lar
  rm(lar)
  
#Calculating accuracy score component (Sa)
  per_base_cov_acc[lar == T, "acc_compl" := 1 - accuracy]
  per_base_cov_acc[, "accuracy" := NULL]
  acc_score_comp_stats <- per_base_cov_acc[lar == T, .(lar_acc_pen = sum(acc_compl) + .N * LAR_EX_PEN), by = c("lar_mark", "transcript")][, .(acc_pen = sum(lar_acc_pen)), by = "transcript"]
  rm(per_base_cov_acc)
  acc_score_comp_stats <- merge(acc_score_comp_stats, tr_cov_len, by = "transcript", all.y = T)
  rm(tr_cov_len)
  acc_score_comp_stats[, "acc_score_comp" := 1 - sqrt(acc_pen / (cov_len * (LAR_EX_PEN + 1)))]
  acc_score_comp_stats[, c("cov_len", "acc_pen") := NULL]
  acc_score_comp_stats[is.na(acc_score_comp), "acc_score_comp" := 1]
 } else {
  acc_score_comp_stats <- data.table(transcript = tr_cov_len[, transcript], acc_score_comp = 1)
  rm(per_base_cov_acc, tr_cov_len)
 }
 if(length(fully_uncov_tr) != 0) {
  acc_score_comp_stats <- rbindlist(list(acc_score_comp_stats, data.table(transcript = fully_uncov_tr, acc_score_comp = NA)))
 }
 rm(fully_uncov_tr)
 accuracy_stats <- merge(accuracy_stats, acc_score_comp_stats, by = "transcript")
 rm(acc_score_comp_stats)
 
#Saving accuracy statistics
 accuracy_stats_dt[[i]] <- accuracy_stats
 rm(accuracy_stats)
}

#Binding tables
per_base_cov_category_dt <- rbindlist(per_base_cov_category_dt)
per_base_higher_cov_dt <- rbindlist(per_base_higher_cov_dt)
lcr_dt <- rbindlist(lcr_dt)
pos_cov_dt <- rbindlist(pos_cov_dt)
coverage_stats_dt <- rbindlist(coverage_stats_dt)
pos_acc_dt <- rbindlist(pos_acc_dt)
per_base_acc_category_dt <- rbindlist(per_base_acc_category_dt)
per_base_higher_acc_dt <- rbindlist(per_base_higher_acc_dt)
lar_dt <- rbindlist(lar_dt)
accuracy_stats_dt <- rbindlist(accuracy_stats_dt)

#Calculating the number of fully covered and fully uncovered transcripts
fully_cov_tr_N_sum <- sum(fully_cov_tr_N_vec)
fully_uncov_tr_N_sum <- sum(fully_uncov_tr_N_vec)

#Calculating the number of bases
base_N_sum <- sum(base_N_vec)

#Calculating the number of covered bases
cov_base_N_sum <- sum(cov_base_N_vec)

#Calculating summaries of coverage and accuracy metrics
cov_prop_summary <- summary(coverage_stats_dt[, cov_prop])
cov_mean_summary <- summary(coverage_stats_dt[, cov_mean])
lcr_prop_summary <- summary(coverage_stats_dt[, lcr_prop])
lcr_len_summary <- summary(lcr_dt[, lcr_len])
cov_score_comp_summary <- summary(coverage_stats_dt[, cov_score_comp])
ucr_len_max_summary <- summary(coverage_stats_dt[, ucr_len_max])
end_cov_mean_summary <- summary(coverage_stats_dt[, end_cov_mean])
acc_base_prop_summary <- summary(accuracy_stats_dt[, acc_base_prop])
lar_prop_summary <- summary(accuracy_stats_dt[, lar_prop])
lar_len_summary <- summary(lar_dt[, lar_len])
acc_score_comp_summary <- summary(accuracy_stats_dt[, acc_score_comp])

#Calculating the number of LCR and LAR bases
lcr_len_sum <- coverage_stats_dt[, sum(lcr_N)]
lar_len_sum <- accuracy_stats_dt[, sum(lar_N)]

#Writing coverage statistics to file
setnames(coverage_stats_dt, c("transcript", "covered_base_N", "covered_base_prop", "covered_base_prop_category", "coverage_mean", "coverage_mean_category", "lcr_base_N", "lcr_base_prop" ,"lcr_base_prop_category", "coverage_score_component", "uncov_region_length_max", "transcript_end_coverage_mean"))
setcolorder(coverage_stats_dt, c("transcript", "covered_base_N", "covered_base_prop", "covered_base_prop_category", "coverage_mean", "coverage_mean_category", "uncov_region_length_max", "transcript_end_coverage_mean", "lcr_base_N", "lcr_base_prop", "lcr_base_prop_category", "coverage_score_component"))
setorder(coverage_stats_dt, transcript)
write.table(coverage_stats_dt, file = paste(OUT_PREF, "coverage_stats.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(coverage_stats_dt)

#Writing per-base coverage to file
per_base_cov_category_dt <- per_base_cov_category_dt[, .(base_N = sum(base_N)), by = "coverage_cat"]
per_base_cov_category_dt[, "base_prop" := base_N / base_N_sum]
per_base_cov_category_dt[, "coverage_cat_ord" := as.numeric(sub(".*?(-?\\d+\\.?\\d*).*", "\\1", coverage_cat))]
setorder(per_base_cov_category_dt, coverage_cat_ord)
per_base_cov_category_dt[, "coverage_cat_ord" := NULL]
setnames(per_base_cov_category_dt, c("coverage_category", "base_N", "base_prop"))
write.table(per_base_cov_category_dt, file = paste(OUT_PREF, "per_base_coverage_distribution.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(per_base_cov_category_dt)

#Forming the table with the number of bases with >= coverage than coverage breaks
per_base_higher_cov_dt <- per_base_higher_cov_dt[, .(base_N = sum(base_N)) , by = "cov_threshold"]
per_base_higher_cov_dt[, "base_prop" := base_N / base_N_sum]

#Writing median relative coverage by transcript position to file
pos_cov_dt <- pos_cov_dt[, .(rel_cov_mean_median = median(rel_cov_mean, na.rm = T)), by = "position_prop_cat"]
pos_cov_dt[, "position_prop_cat_ord" := as.numeric(sub(".*?(-?\\d+\\.?\\d*).*", "\\1", position_prop_cat))]
setorder(pos_cov_dt, position_prop_cat_ord)
pos_cov_dt[, "position_prop_cat_ord" := NULL]
setnames(pos_cov_dt, c("position_prop_category", "relative_coverage_median"))
write.table(pos_cov_dt, file = paste(OUT_PREF, "relative_coverage_median_by_transcript_position.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(pos_cov_dt)

#Writing LCR coordinates to file
if(lcr_dt[, .N] != 0) {
 setnames(lcr_dt, c("transcript", "lcr_start", "lcr_end", "lcr_length"))
 setorder(lcr_dt, transcript)
 write.table(lcr_dt, file = paste(OUT_PREF, "lcr_list.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
} else {
 write.table(data.table(), file = paste(OUT_PREF, "lcr_list.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
}
rm(lcr_dt)

#Writing accuracy statistics to file
setnames(accuracy_stats_dt, c("transcript", "acc_base_N", "acc_base_prop", "acc_base_prop_category", "lar_base_N", "lar_base_prop", "lar_base_prop_category", "accuracy_score_component"))
setorder(accuracy_stats_dt, transcript)
write.table(accuracy_stats_dt, file = paste(OUT_PREF, "accuracy_stats.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)

#Calculating the number of transcripts
tr_N <- accuracy_stats_dt[, .N]
rm(accuracy_stats_dt)

#Writing per-base accuracy to file
per_base_acc_category_dt <- per_base_acc_category_dt[, .(base_N = sum(base_N)), by = "accuracy_cat"]
per_base_acc_category_dt[, "base_prop" := base_N / cov_base_N_sum]
per_base_acc_category_dt[, "accuracy_cat_ord" := as.numeric(sub(".*?(-?\\d+\\.?\\d*).*", "\\1", accuracy_cat))]
setorder(per_base_acc_category_dt, accuracy_cat_ord)
per_base_acc_category_dt[, "accuracy_cat_ord" := NULL]
setnames(per_base_acc_category_dt, c("accuracy_category", "base_N", "base_prop"))
write.table(per_base_acc_category_dt, file = paste(OUT_PREF, "per_base_accuracy_distribution.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(per_base_acc_category_dt)

#Forming the table with the number of bases with >= accuracy than accuracy breaks
per_base_higher_acc_dt <- per_base_higher_acc_dt[, .(base_N = sum(base_N)), by = "acc_threshold"]
per_base_higher_acc_dt[, "base_prop" := base_N / cov_base_N_sum]

#Writing median accuracy by transcript position to file
pos_acc_dt <- pos_acc_dt[, .(acc_mean_median = median(acc_mean, na.rm = T)), by = "position_prop_cat"]
pos_acc_dt[, "position_prop_cat_ord" := as.numeric(sub(".*?(-?\\d+\\.?\\d*).*", "\\1", position_prop_cat))]
setorder(pos_acc_dt, position_prop_cat_ord)
pos_acc_dt[, "position_prop_cat_ord" := NULL]
setnames(pos_acc_dt, c("position_prop_category", "accuracy_median"))
write.table(pos_acc_dt, file = paste(OUT_PREF, "accuracy_median_by_transcript_position.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
rm(pos_acc_dt)

#Writing LAR coordinates to file
if(lar_dt[, .N] != 0) {
 setnames(lar_dt, c("transcript", "lar_start", "lar_end", "lar_length"))
 setorder(lar_dt, transcript)
 write.table(lar_dt, file = paste(OUT_PREF, "lar_list.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
} else {
 write.table(data.table(), file = paste(OUT_PREF, "lar_list.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)
}
rm(lar_dt)

#Forming the vector with the number of bases with >= coverage than coverage breaks
generate_cov_break_labels <- c()
for(i in 1 : length(BASE_COV_BREAKS2)) {
 generate_cov_break_labels[i] <- paste("N, % bases with coverage equal to or higher than", BASE_COV_BREAKS2[i])
}

#Forming the vector with the number of bases with >= accuracy than accuracy breaks
generate_acc_break_labels <- c()
for(i in 1 : length(BASE_ACC_BREAKS2[])) {
 generate_acc_break_labels[i] <- paste("N, % bases with accuracy equal to or higher than", BASE_ACC_BREAKS2[i])
}

#Writing coverage and accuracy summary to file
write.table(data.table(parameter = c("% of covered bases per transcript (mean, IQR)", "N, % fully covered transcripts", "N, % fully uncovered transcripts", "Mean coverage per transcript (mean, IQR)", generate_cov_break_labels, "Maximum uncovered region length per transcript (mean, IQR)", "Mean end coverage per transcript (mean, IQR)", "N, % assembly bases in LCR", "% of bases in LCR per transcript (mean, IQR)", "LCR length (mean, IQR)", "Coverage score component (mean, IQR)", paste0("% of accurate bases (bases with accuracy higher than or equal to ", BASE_ACC_THR, ") per transcript (mean, IQR)"), generate_acc_break_labels, "N, % assembly bases in LAR", "% of bases in LAR per transcript (mean, IQR)", "LAR length (mean, IQR)", "Accuracy score component  (mean, IQR)"), value = c(paste0(round(100 * cov_prop_summary[4], 2), "%, ", round(100 * cov_prop_summary[2], 2), "%-", round(100 * cov_prop_summary[5], 2), "%"), paste0(fully_cov_tr_N_sum, ", ", round(100 * fully_cov_tr_N_sum / tr_N, 2), "%"), paste0(fully_uncov_tr_N_sum, ", ", round(100 * fully_uncov_tr_N_sum / tr_N, 2), "%"), paste0(round(cov_mean_summary[4], 2), ", ", round(cov_mean_summary[2], 2), "-", round(cov_mean_summary[5], 2)), paste0(per_base_higher_cov_dt[, base_N], ", ", round(100 * per_base_higher_cov_dt[, base_prop], 2), "%"), paste0(round(ucr_len_max_summary[4], 2), ", ", ucr_len_max_summary[2], "-", ucr_len_max_summary[5]), paste0(round(end_cov_mean_summary[4], 2), ", ", round(end_cov_mean_summary[2], 2), "-", round(end_cov_mean_summary[5], 2)), paste0(lcr_len_sum, ", ", round(100 * lcr_len_sum / base_N_sum, 2), "%"), paste0(round(100 * lcr_prop_summary[4], 2), "%, ", round(100 * lcr_prop_summary[2], 2), "%-", round(100 * lcr_prop_summary[5], 2), "%"), paste0(round(lcr_len_summary[4], 2), ", ", lcr_len_summary[2], "-", lcr_len_summary[5]), paste0(round(cov_score_comp_summary[4], 3), ", ", round(cov_score_comp_summary[2], 3), "-", round(cov_score_comp_summary[5], 3)), paste0(round(100 * acc_base_prop_summary[4], 2), "%, ", round(100 * acc_base_prop_summary[2], 2), "%-", round(100 * acc_base_prop_summary[5], 2), "%"), paste0(per_base_higher_acc_dt[, base_N], ", ", round(100 * per_base_higher_acc_dt[, base_prop], 2), "%"), paste0(lar_len_sum, ", ", round(100 * lar_len_sum / cov_base_N_sum, 2), "%"), paste0(round(100 * lar_prop_summary[4], 2), "%, ", round(100 * lar_prop_summary[2], 2), "%-", round(100 * lar_prop_summary[5], 2), "%"), paste0(round(lar_len_summary[4], 2), ", ", lar_len_summary[2], "-", lar_len_summary[5]), paste0(round(acc_score_comp_summary[4], 3), ", ", round(acc_score_comp_summary[2], 3), "-", round(acc_score_comp_summary[5], 3)))), file = paste(OUT_PREF, "coverage_and_accuracy_analysis_summary.tsv", sep = "_"), sep = "\t", row.names = F, col.names = F, quote = F)
