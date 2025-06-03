#!/usr/bin/env Rscript
#Transcriptome assembly score generation script
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
COV_STATS_DT_PATH <- ext_args[2]
OUT_PREF <- ext_args[3]

#Importing transcript score components from files
setDTthreads(THREAD_N)
tr_scores <- fread(COV_STATS_DT_PATH, header = T, select = c(1, 12))
for (i in c("accuracy_stats.tsv", "local_fidelity_stats.tsv", "integrity_stats.tsv")) {
 import_filename <- paste(OUT_PREF, i, sep = "_")
 if (i == "accuracy_stats.tsv") {
  selected_cols <- c(1, 8)
 } else if (i == "local_fidelity_stats.tsv") {
  selected_cols <- c(1, 13)
 } else {
  selected_cols <- c(1, 7)
 }
 tr_scores <- merge(tr_scores, fread(import_filename, header = T, select = selected_cols), by = "transcript")
}

#Calculating transcript scores
tr_scores[, "transcript_score" := coverage_score_component * accuracy_score_component * local_fidelity_score_component * integrity_score_component]
tr_scores[is.na(transcript_score), "transcript_score" := 0]

#Calculating score component summaries
coverage_score_component_summary <- summary(tr_scores[, coverage_score_component])
accuracy_score_component_summary <- summary(tr_scores[, accuracy_score_component])
local_fidelity_score_component_summary <- summary(tr_scores[, local_fidelity_score_component])
integrity_score_component_summary <- summary(tr_scores[, integrity_score_component])

#Calculating transcriptome assembly score
assembly_score <- tr_scores[, mean(transcript_score)]

#Writing transcript scores to file
write.table(tr_scores, file = paste(OUT_PREF, "transcript_scores.tsv", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)

#Writing transcriptome assembly score summary to file
write.table(data.table(parameter = c("Coverage score component (mean, IQR)", "Accuracy score component (mean, IQR)", "Local fidelity score component (mean, IQR)", "Integrity score component (mean, IQR)", "Assembly score"), value = c(paste0(round(coverage_score_component_summary[4], 3), ", ", round(coverage_score_component_summary[2], 3), "-", round(coverage_score_component_summary[5], 3)), paste0(round(accuracy_score_component_summary[4], 3), ", ", round(accuracy_score_component_summary[2], 3), "-", round(accuracy_score_component_summary[5], 3)), paste0(round(local_fidelity_score_component_summary[4], 3), ", ", round(local_fidelity_score_component_summary[2], 3), "-", round(local_fidelity_score_component_summary[5], 3)), paste0(round(integrity_score_component_summary[4], 3), ", ", round(integrity_score_component_summary[2], 3), "-", round(integrity_score_component_summary[5], 3)), round(assembly_score, 3))), file = paste(OUT_PREF, "assembly_score_summary.tsv", sep = "_"), sep = "\t", row.names = F, col.names = F, quote = F)