# CATS-rf

<img src="cats_rf_logo.png" alt="Logo" width="750" height="160"/>

# Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Test data](#test-data)
- [Example usage](#example-usage)
- [Detailed options](#detailed-options)
- [Output explanation](#output-explanation)
- [Assembly comparison with `CATS_rf_compare`](#assembly-comparison-with-cats_rf_compare)
- [Citation](#citation)
- [Troubleshooting](#troubleshooting)
- [Changelog](#changelog)

# Introduction 

CATS-rf is the reference-free module of the CATS (Comprehensive Assessment of Transcript Sequences) framework. It evaluates the quality of transcriptomes assembled from short reads without using reference data, relying solely on RNA-seq reads used in the assembly construction. The pipeline maps reads back to the assembled transcripts and examines mapping evidence suggesting misassembly. Quality evaluation is performed at the transcript level, integrating four score components each targeting specific assembly errors:

| **Score Component**                                            | **Evidence**                                 | **Targeted Assembly Errors**                                      |
|----------------------------------------------------------------|----------------------------------------------|-------------------------------------------------------------------|
| Coverage&nbsp;component&nbsp;(<i>S<sub>c</sub></i>)            | Low-coverage regions                         | Insertions, redundancy                                            |
| Accuracy&nbsp;component&nbsp;(<i>S<sub>a</sub></i>)            | Low-accuracy regions                         | Sequence inaccuracy                                               |
| Local&nbsp;fidelity&nbsp;component&nbsp;(<i>S<sub>l</sub></i>) | Inconsistent pair mapping within transcripts | Structural errors (e.g. deletions, translocations, inversions...) |
| Integrity&nbsp;component&nbsp;(<i>S<sub>i</sub></i>)           | Pairs mapping to different transcripts       | Transcript fragmentation                                          |

Transcript quality score S<sub><i>t</i></sub> is calculated as the product of the described score components, equally weighting detected assembly errors. Assembly score <i>S</i> is computed as the mean of individual transcript scores.

In addition to transcript scores, CATS-rf provides a comprehensive set of assembly metrics, including transcript length and composition statistics, read mapping rates, positional coverage and accuracy profiles, and pair mapping consistency metrics.

CATS-rf consistently displays stronger performance than currently existing reference-free transcriptome assembly evaluation tools. For detailed benchmarks and methodology, please refer to the CATS [preprint](test)


# Installation 

## Compatibility

### Linux and Windows

For the best compatibility and performance, we recommend running CATS-rf on:
- Any modern Linux distribution (e.g. Ubuntu, Debian, Fedora, etc.)
- WSL (i.e. Ubuntu on Windows)

### MacOS
If you are using MacOS, Bash (version >= 4.0) and GNU versions of core utilities are required. In this case, `PATH` variable should be adjusted so that CATS-rf uses GNU versions of core utilities:

- Install Bash ≥ 4.0 via [Homebrew](https://formulae.brew.sh/formula/bash):

```bash
brew install bash
```

- Install GNU core utilities:

```bash
brew install coreutils findutils gnu-sed gawk grep
```

- Add Bash and GNU utilities to your `PATH` (adjust path depending on your architecture):

For Apple Silicon:
```bash
export PATH="/opt/homebrew/bin:$PATH"
export PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/findutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/grep/libexec/gnubin:$PATH"
```

For Intel-based configurations:
```bash
export PATH="/usr/local/bin:$PATH"
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/findutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/grep/libexec/gnubin:$PATH"
```

- Run CATS-rf using the installed Bash version:
  
```bash
bash CATS_rf
```

## Installing CATS-rf via conda

CATS-rf and its dependencies can be directly installed via [Bioconda](https://bioconda.github.io/):

```bash
conda install -c bioconda cats_rf
```

## Installing CATS-rf from source

CATS-rf consists of Bash and R scripts located in the `scripts` directory of this repository. After cloning the repository, all CATS-rf scripts must be included in the `PATH` environment variable. 

The following dependencies are required:

| **Dependency**      | **Tested Version** | **Homepage**                                   | **Conda Installation**                    |
|---------------------|--------------------|------------------------------------------------|-------------------------------------------|
| R                   | 4.4.2              | https://www.r-project.org/                     | `conda install conda-forge::r-base`       |
| data.table (R)      | 1.16.4             | https://cran.r-project.org/package=data.table  | `conda install conda-forge::r-data.table` |
| Bowtie2             | 2.5.4              | https://github.com/BenLangmead/bowtie2         | `conda install -c bioconda bowtie2`       |
| Samtools            | 1.21               | https://www.htslib.org/                        | `conda install -c bioconda samtools`      |
| kallisto            | 0.50.1             | https://github.com/pachterlab/kallisto         | `conda install -c bioconda kallisto`      |
| GNU Parallel        | 20220922           | https://www.gnu.org/software/parallel/         | `conda install conda-forge::parallel`     |
| bedtools (bamToBed) | 2.31.1             | https://github.com/arq5x/bedtools2             | `conda install -c bioconda bedtools`      |
| pysamstats          | 1.1.2              | https://github.com/alimanfoo/pysamstats        | `conda install -c bioconda pysamstats`    |

R, Bowtie2, Samtools, kallisto, GNU Parallel, bedtools (bamToBed), and pysamstats executables must be included in `PATH`. R package data.table can be installed via conda or directly in R with `install.packages("data.table")`

# Test data

CATS-rf installation can be tested using instructions and files located in `test_data` directory.

# Example usage 

CATS-rf requires a transcriptome assembly in FASTA format, along with short RNA-seq reads used during assembly in either FASTQ or FASTA format. Compressed (.gz) read files are supported.

CATS-rf supports both paired-end and single-end library configurations.

Example paired-end mode usage:

```bash
CATS_rf [OPTIONS] TRANSCRIPTOME READS1 READS2
```

Example single-end mode usage:

```bash
CATS_rf -C se -m MEAN_INS_SIZE -s SD_INS_SIZE [OTHER_OPTIONS] TRANSCRIPTOME READS1
```

Single-end mode requires three options to be specified: `C` for library configuration, `m` for mean fragment size, and `s` for standard deviation of fragment size. Note that single-end runs will output only general assembly statistics, read mapping metrics, and positional coverage and accuracy analysis.

# Detailed options

CATS-rf offers a comprehensive list of options which allow users to control the analysis parameters.

## Library type options

`-C`: Paired- vs. single-end library configuration: pe = paired-end, se = single-end, default: pe

`-S`: Library strandness, fr = forward-reverse, rf = reverse-forward, u = unstranded, a = automatic detection, default: u

CATS-rf can leverage strandness information when quantifying transcripts and calculating local fidelity score component. When the automatic detection option is enabled, strandness is estimated using the first 100 000 read mappings.

While CATS-rf was primarily tested on Illumina data, the analysis can be run on assemblies generated from other short-read platforms. In such scenario, `S` should be adjusted accordingly. If strandness of the data is unknown, it is recommended to use either unstranded mode or automatic detection. Note that in unstranded mode, read pairs are expected to map to opposite strands. This is consistent with the behavior of virtually every short-read sequencing technology.

`-Q`: Phred quality encoding of FASTQ files, 33 = phred33, 64 = phred64, default: 33

## Read mapping, transcript quantification, and read assignment options

`-R`: Random seed for read mapping, transcript quantification, and read assignment, default: 12345

Random seed is defined to ensure reproducible CATS-rf runs.

`-N`: Maximum number of distinct mappings per read, default: 10

The value of `N` should be increased for complex transcriptome assemblies that contain a large number of isoforms, and decreased for simpler assemblies with fewer isoforms to maximize performance and accuracy. Note that Bowtie2 mapping parameters are optimized to detect transcript errors, while minimizing the number of false-positive mappings. Furthermore, secondary mappings of each read are filtered based on edit distance.

`-m`: Estimated mean of fragment length needed for transcript quantification (single-end mode only)

`-s`: Estimated standard deviation of fragment length needed for transcript quantification (single-end mode only)

Fragment length distribution parameters `m` and `n` are required in single-end mode for transcript quantification by kallisto.

## Coverage analysis options

`-i`: Per-base coverage distribution breakpoints (specified with x,y,z...), default: "0,5,10,20,40,60,80,100"

Per-base coverage is split into intervals defined by `i` (e.g. [0-5>, [5-10>...). This category variable is used for plotting by the `CATS_rf_compare` script.

All category variable breaks (`i`, `p`, `r`, `u`, `I`, `P`, `U`, `y`, and `F`) should be supplied as strings separated with commas and enclosed in quotes (e.g. "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1").

`-p`: Per-transcript proportion of covered bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of covered bases is split into intervals defined by `p` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-r`: Mean transcript coverage distribution breakpoints (specified with x,y,z...), default: "0,5,10,20,40,60,80,100"

Mean transcript coverage is split into intervals defined by `r` (e.g. [0-5>, [5-10>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-l`: Proportion of transcript length for positional relative coverage distribution analysis, default: 0.01

Transcripts are split into fractional segments of size `l` for positional relative coverage distribution analysis. Coverage is expressed relative to the base with the highest coverage within the same transcript. Relative coverage for each segment is calculated as mean relative coverage within the segment. Positional analysis output contains assembly-level median relative coverage for each transcript segment.

`-n`: Proportion of transcript length for transcript end definition when calculating mean transcript end coverage, default: 0.02

Relative size of transcript end regions when calculating mean transcript end coverage is controlled by `n`.

`-k`: Rolling window size for local coverage calculation (in bp) when defining low-coverage regions (LCR), default: 10

`-z`: Local coverage threshold for LCR characterization, default: 3

LCRs are defined as rolling windows of size `k` with mean coverage lower than or equal to `z`.

`-u`: Per-transcript proportion of LCR bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of LCR bases is split into intervals defined by `u` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-w`: Base coverage weight, default: 1.5

`-e`: LCR extension penalty, default: 0.5

Coverage penalties assigned to LCRs are controlled by `w` and `e`. Lower values of `w` and higher values of `e` increase the relative impact of LCR length on coverage penalty.

## Accuracy analysis options

`-I`: Per-base accuracy distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Accuracy is defined as the proportion of aligned read bases matching the transcript base. Per-base accuracy is split into intervals defined by `I` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-A`: Minimum accuracy for a base to be considered accurate, default: 0.95

`-P`: Per-transcript proportion of accurate bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Per-transcript proportion of accurate bases (bases with accuracy higher or equal to `A`) is split into intervals defined by `P` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-L`: Proportion of transcript length for positional accuracy distribution analysis, default: 0.01

Transcripts are split into fractional segments of size `L` for positional accuracy distribution analysis. Accuracy for each segment is calculated as mean accuracy within the segment. Positional analysis output contains assembly-level median accuracy for each transcript segment.

`-K`: Rolling window size for local accuracy calculation (in bp) when defining low-accuracy regions (LAR), default: 10

`-Z`: Local accuracy threshold for LAR characterization, default: 0.98

LARs are defined as rolling windows of size `K` with mean accuracy lower than or equal to `Z`.

`-U`: Per-transcript proportion of LAR bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Per-transcript proportion of LAR bases is split into intervals defined by `U` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-E`: LAR extension penalty, default: 0.1

Accuracy penalties assigned to LARs are controlled with `E`. Higher values of `E` increase the relative impact of LAR length on accuracy penalty.

## Paired-end read analysis options

These options should only be supplied in paired-end mode.

`-d`: Maximum distance from transcript ends for reads with unmapped pair to be considered evidence of transcript end incompleteness or fragmentation (in bp), default: 40

Reads with unmapped pair mapping to transcript ends are considered evidence for transcript end incompleteness or fragmentation. Relative size of transcript end regions when identifying such reads is controlled by `d`.

`-x`: Multiplicative factor for lower distance outlier threshold calculation, default: 8

`-X`: Multiplicative factor for higher distance outlier threshold calculation, default: 10

`-c`: Correction factor for distance outlier threshold calculation, default: 5

Read pair distance penalty calculation is controlled by `x`, `X`, and `c`. Read pairs are classified as mapping too far apart if their distance exceeds the lower distance threshold, defined as D<sub>1</sub> = Q<sub>3</sub>(d) + x * (IQR(d) + c). These reads are assigned a distance penalty P<sub>d</sub> = d / D<sub>2</sub>, where D<sub>2</sub> = Q<sub>3</sub>(d) + X * (IQR(d) + c), with the penalty capped at 1. Higher values of `x` increase the threshold for classifying read pairs as too distant, while `X` controls the scaling of the distance penalty. Higher values of `c` increase penalty robustness in libraries with a high proportion of overlapping read pairs.

`-y`: Per-transcript proportion of improperly paired reads within a transcript distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Improperly paired reads include reads with pair not mapped to the assembly, reads with pair mapped in an unexpected orientation, and reads with pair mapped too far apart. Per-transcript proportion of improperly paired reads within a transcript is split into intervals defined by `y` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-f`: Minimum number of bridging events for transcripts to be considered fragmented, default: 3

A transcript is considered fragmented if more than `f` reads representing bridging events map to transcript end regions.

`-F`: Per-transcript proportion of reads with pair mapped to another transcript distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of reads with pair mapped to another transcript is split into intervals defined by `F` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the `CATS_rf_compare` script.

`-a`: Alpha compression factor for sigmoid transformation applied to bridge index during integrity score component calculation, default: 7

`-b`: Beta compression factor for sigmoid transformation applied to bridge index during integrity score component calculation, default: 0.5

Bridge index measures the proportion of reads with pair mapped to a different transcript and considers the mapping distance of such reads from the ends of their respective transcript. This definition gives more weight to bridging events near transcript ends. Integrity score component is calculated using a sigmoid transformation of bridge index. Compression factors `a` and `b` control the shape of the transformation: higher values of `a` increase sensitivity to fragmentation, while higher values of `b` reduce the likelihood of false-positive fragmentation penalties in transcripts with minimal bridging evidence.

## General options

`-t`: Number of CPU threads, default: 10

Several steps of CATS-rf pipeline are parallelized. This includes read mapping, transcript quantification, read assignment, SAM/BAM file processing, positional coverage and accuracy calculation and analysis, as well as positional paired-end analysis. Recommended number of threads: 10-20.

`-G`: Percentage of available RAM used by GNU sort, default: 50

CATS-rf utilizes GNU sort in several steps of the pipeline. Higher values of `G` will ensure faster sorting, but may exhaust available RAM. In such scenarios, CATS-rf will resort to sorting with minimal RAM usage.

`-M`: Memory block size for GNU Parallel, default: 512M

Block size used by GNU Parallel when splitting the mapping table for read assignment is controlled by `M`. If sufficient RAM is available, increasing the value of `M` is recommended to minimize artifacts introduced by file splitting.

`-T`: Number of splits performed on positional and read pair mapping tables, default: 3

Positional and read pair mapping tables are split before analysis to reduce RAM usage. Increase the value of `T` when working with limited memory to further reduce RAM demands.

`-D`: CATS-rf output directory name, default: TRANSCRIPTOME_CATS_rf_dir

`-o`: CATS-rf output file prefix, default: TRANSCRIPTOME

`-O`: Overwrite the CATS-rf output directory, default: off

`-h`: Show usage information

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

# Assembly comparison with `CATS_rf_compare`

CATS-rf also supports direct comparison of multiple analysed assemblies. The `CATS_rf_compare` script generates summary tables and visualizations that compare the most significant CATS-rf results of each assembly. As such, CATS-rf should be run on each individual assembly and the resulting CATS-rf output directories should then act as input to `CATS_rf_compare`.

## `CATS_rf_compare` dependencies

`CATS_rf_compare` requires the following dependencies:

| **Dependency**           | **Tested Version** | **Homepage**                                     | **Conda Installation**                    | **R installation**               |
|--------------------------|--------------------|--------------------------------------------------|-------------------------------------------|----------------------------------|
| R                        | 4.4.2              | https://www.r-project.org/                       | `conda install conda-forge::r-base`       | /                                |
| knitr (R)                | 1.49               | https://cran.r-project.org/web/packages/knitr/   | `conda install conda-forge::r-knitr`      | `install.packages("knitr")`      |
| data.table (R)           | 1.16.4             | https://cran.r-project.org/package=data.table    | `conda install conda-forge::r-data.table` | `install.packages("data.table")` |
| ggplot2 (R)              | 3.5.1              | https://cran.r-project.org/web/packages/ggplot2/ | `conda install conda-forge::r-ggplot2`    | `install.packages("ggplot2")`    |
| ggdist (R)               | 3.3.2              | https://cran.r-project.org/web/packages/ggdist/  | `conda install conda-forge::r-ggdist`     | `install.packages("ggdist")`     | 

R (Rscript) executable must be included in `PATH`. Tools denoted with (R) correspond to R packages and can be installed via conda or directly in R with the supplied commands.

## `CATS_rf_compare` example usage 

`CATS_rf_compare` requires one or more CATS-rf output directories as input.

While `CATS_rf_compare` is primarily designed to compare multiple transcriptome assemblies, it can also be used with a single assembly to visualize its CATS-rf results.

Example `CATS_rf_compare` usage:

```bash
CATS_rf_compare [OPTIONS] CATS_RF_DIR ...
```

## Detailed `CATS_rf_compare` options

`CATS_rf_compare` offers a comprehensive list of options which allow users to control the graphical and general comparison parameters.

### Graphical options

`-x`: Figure extension, default: png

`-d`: Figure DPI, default: 600

Extension (device) and DPI of each plotted figure are controlled with `x` and `d`, respectively.

`-r`: Raincloud plot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

Raincloud plot densities are normalized for each transcriptome assembly. Boxplots within raincloud plots mark the distribution median, Q<sub>1</sub>, and Q<sub>3</sub>, with whiskers extending from Q<sub>1</sub> - 1.5 * IQR to Q<sub>3</sub> + 1.5 * IQR of the distribution.

All color sets (`r`, `l`, `H`, and `b`) should be supplied as R color names or hexadecimal codes separated with commas and enclosed in quotes (e.g. "#FDAF4A,#DC151D"). R color cheatsheet is available [here](https://sites.stat.columbia.edu/tzheng/files/Rcolor.pdf).

`-l`: Lineplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

`-H`: Histogram colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted Set1 palette from RColorBrewer package

`-b`: Barplot colors (quoted hexadecimal codes or R color names, specified with x,y,z...), default: adjusted YlOrRd palette from RColorBrewer package

`-q`: Maximum right-tail distribution quantile for histograms, default: 0.98"

Histograms show relative density per transcriptome assembly and omit right-tail extreme values for visualization purposes. The x-axis in all histograms is square-root scaled.

### General options

`-t`: Number of CPU threads, default: 10

Several steps of `CATS_rf_compare` are parallelized. This mainly includes operations performed by the data.table package. Recommended number of threads: 8-12.

`-D`: Comparison output directory name, default: CATS_rf_comparison

`-O`: Overwrite the comparison output directory, default: off

`-h`: Show usage information

## `CATS_rf_compare` output explanation

The analysis is summarized in the `CATS_rf_comparison.html` HTML file. 
An example of the HTML output is provided [here](CATS_rf_compare_output_example.html).

### Summary tables

`CATS_rf_compare` aggregates individual summary tables into comprehensive joint tables encompassing all analyzed transcriptome assemblies:

`CATS_rf_general_statistics.tsv` contains aggregated CATS-rf general statistics table.

`CATS_rf_assembly_scores.tsv` contains aggregated CATS-rf score component statistics and overall assembly score table.

`CATS_rf_coverage_accuracy_statistics.tsv` contains aggregated CATS-rf coverage and accuracy analysis table.

`CATS_rf_local_fidelity_integrity_statistics.tsv` contains aggregated CATS-rf paired-end read analysis table.

### Figures

`CATS_rf_compare` produces several figures, providing a detailed visualization of CATS-rf quality metrics. 

`transcript_score` visualizes the distribution of transcript scores.

`base_coverage` and `base_accuracy` visualize the distribution of per-base coverage/accuracy.

`proportion_of_covered_bases` visualizes the distribution of the proportion of covered bases per transcript.

`mean_transcript_coverage` visualizes the distribution of mean transcript coverage.

`positional_relative_coverage_median` and `positional_accuracy_median` visualize the positional relative coverage/accuracy distribution.

`maximum_uncovered_region_length` visualizes the distribution of maximum uncovered region length per transcript.

`mean_transcript_end_coverage` visualizes the distribution of mean transcript end coverage.

`proportion_of_bases_in_lcrs` and `proportion_of_bases_in_lars` visualize the distribution of the proportion of bases in LCRs/LARs per transcript.

`lcr_length` and `lar_length` visualize the distribution of LCR/LAR length.

`coverage_score_component` visualizes the distribution of coverage score component per transcript.

`proportion_of_accurate_bases` visualizes the distribution of the proportion of accurate bases per transcript.

`accuracy_score_component` visualizes the distribution of accuracy score component per transcript.

`proportion_of_improperly_paired_reads` visualizes the per-transcript distribution of the proportion of improperly paired reads within a transcript.

`local_fidelity_score_component` visualizes the distribution of local fidelity score component per transcript.

`prop_reads_with_pair_mapped_to_another_tr` visualizes the per-transcript distribution of the proportion of reads with pair mapped to another transcript.

`integrity_score_component` visualizes the distribution of integrity score component per transcript.

# Citation

CATS is an academic software distributed under the MIT license. 

Copyright © 2025 Kristian Bodulić

if you use CATS, please cite the CATS preprint:

(add reference)

# Troubleshooting

Please report all potential bugs in the Issues tracker.

# Changelog

Version 1.0.0. Initial commit, June 3, 2025.
