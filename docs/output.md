# Output explanation

## Summary table

Summary files represent the main output of CATS-rf. In paired-end mode, four summary files are produced:

`assembly_score_summary.tsv` contains descriptive statistics of transcript score components and the overall assembly score. The content of this file is printed when CATS-rf finishes running in paired-end mode.

`general_statistics_table.tsv` contains descriptive statistics of transcript length (mean, median, interquartile range, range, N50, L50, N90, L90), GC content, and read mapping rate.

`coverage_and_accuracy_analysis_summary.tsv` contains summarized results of coverage and accuracy analysis. The content of this file is printed when CATS-rf finishes running in single-end mode.

`paired_end_read_analysis_summary.tsv` contains summarized results of paired-end read analysis, including local fidelity and integrity analysis.

CATS-rf also produces several .tsv files containing detailed per-transcript metrics:

## Transcript score components

`transcript_scores.tsv` contains CATS-rf score components and transcript score for each transcript.

## Coverage analysis

`coverage_stats.tsv` contains coverage analysis results for each transcript:

| **Column**                    | **Description**                                      |
|-------------------------------|------------------------------------------------------|
| `transcript`                  | Transcript name                                      |
| `covered_base_N`              | Number of covered bases                              |
| `covered_base_prop`           | Proportion of covered bases                          |
| `covered_base_prop_category`  | Proportion of covered bases category                 |
| `coverage_mean`               | Mean transcript coverage                             |
| `coverage_mean_category`      | Mean transcript coverage category                    |
| `uncov_region_length_max`     | Maximum uncovered region length                      |
| `transcript_end_coverage_mean`| Mean transcript end coverage                         |
| `lcr_base_N`                  | Number of bases in low-coverage regions              |
| `lcr_base_prop`               | Proportion of bases in low-coverage regions          |
| `lcr_base_prop_category`      | Proportion of bases in low-coverage regions category |
| `coverage_score_component`    | Coverage score component                             |

`per_base_coverage_distribution.tsv` contains distribution of assembly-level per-base coverage.

`relative_coverage_median_by_transcript_position.tsv` contains median values of mean relative coverage per transcript fraction.

`lcr_list.tsv` contains low-coverage region coordinates.

## Accuracy analysis

`accuracy_stats.tsv` contains accuracy analysis results for each transcript:

| **Column**                 | **Description**                                      |
|----------------------------|------------------------------------------------------|
| `transcript`               | Transcript name                                      |
| `acc_base_N`               | Number of accurate bases                             |
| `acc_base_prop`            | Proportion of accurate bases                         |
| `acc_base_prop_category`   | Proportion of accurate bases category                |
| `lar_base_N`               | Number of bases in low-accuracy regions              |
| `lar_base_prop`            | Proportion of bases in low-accuracy regions          |
| `lar_base_prop_category`   | Proportion of bases in low-accuracy regions category |
| `accuracy_score_component` | Accuracy score component                             |

`per_base_accuracy_distribution.tsv` contains distribution of assembly-level per-base accuracy.

`accuracy_median_by_transcript_position.tsv` contains median values of mean accuracy per transcript fraction.

`lar_list.tsv` contains low-accuracy region coordinates.

## Local fidelity analysis

`local_fidelity_stats.tsv` contains local fidelity analysis results for each transcript:

| **Column**                                 | **Description**                                                              |
|--------------------------------------------|------------------------------------------------------------------------------|
| `transcript`                               | Transcript name                                                              |
| `unmapped_pair_read_N`                     | Number of reads with pair not mapped to the assembly                         |
| `unmapped_pair_read_prop`                  | Proportion of reads with pair not mapped to the assembly                     |
| `unmapped_pair_tr_end_read_N`              | Number of reads with pair not mapped to the assembly on transcript ends      |
| `unmapped_pair_tr_end_read_prop`           | Proportion of reads with pair not mapped to the assembly on transcript ends  |
| `improp_pair_orientation_read_N`           | Number of reads with pair mapped in an unexpected orientation                |
| `improp_pair_orientation_read_prop`        | Proportion of reads with pair mapped in an unexpected orientation            |
| `improp_pair_distance_read_N`              | Number of reads with pair mapped too far apart                               |
| `transcript_distance_penalty`              | Transcript distance penalty                                                  |
| `improp_pair_within_tr_read_N`             | Number of improperly paired reads within a transcript                        |
| `improp_pair_within_tr_read_prop`          | Proportion of improperly paired reads within a transcript                    |
| `improp_pair_within_tr_read_prop_category` | Proportion of improperly paired reads within a transcript category           |
| `local_fidelity_score_component`           | Local fidelity score component                                               |

`read_pairs_mapping_in_unexpected_orientation.tsv` contains coordinates of read pairs mapping in an unexpected orientation.

`read_pairs_mapping_too_far_apart.tsv` contains coordinates of read pairs mapping too far apart.

## Integrity analysis

`integrity_stats.tsv` contains integrity analysis results for each transcript:

| **Column**                              | **Description**                                                      |
|-----------------------------------------|----------------------------------------------------------------------|
| `transcript`                            | Transcript name                                                      |
| `pair_mapped_to_other_tr_N`             | Number of reads with pair mapped to another transcript               |
| `pair_mapped_to_other_tr_prop`          | Proportion of reads with pair mapped to another transcript           |
| `pair_mapped_to_other_tr_prop_category` | Proportion of reads with pair mapped to another transcript category  |
| `bridge_N`                              | Number of reads representing bridging events                         |
| `bridge_prop`                           | Proportion of reads representing bridging events on transcript ends  |
| `integrity_score_component`             | Integrity score component                                            |

`read_pairs_mapping_to_different_transcripts.tsv` contains coordinates of read pairs mapping to different transcripts.
