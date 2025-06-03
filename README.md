# CATS-rf

<img src="cats_rf_logo.png" alt="Logo" width="750" height="160"/>

# Introduction 

CATS-rf is the reference-free module of the CATS (Comprehensive Assessment of Transcript Sequences) framework. It evaluates transcriptome assembly quality without using reference data, relying solely on the short-read RNA-seq reads used in the assembly construction. The pipeline maps reads back to the assembled transcripts and examines mapping evidence suggesting misassembly. Quality evaluation is performed at the transcript level, integrating four score components each targeting specific assembly errors:

| **Score Component**                                            | **Evidence**                                 | **Targeted Assembly Errors**                                      |
|----------------------------------------------------------------|----------------------------------------------|-------------------------------------------------------------------|
| Coverage&nbsp;component&nbsp;(<i>S<sub>c</sub></i>)            | Low-coverage regions                         | Insertions, redundancy                                            |
| Accuracy&nbsp;component&nbsp;(<i>S<sub>a</sub></i>)            | Low-accuracy regions                         | Sequence inaccuracy                                               |
| Local&nbsp;fidelity&nbsp;component&nbsp;(<i>S<sub>l</sub></i>) | Inconsistent pair mapping within transcripts | Structural errors (e.g. deletions, translocations, inversions...) |
| Integrity&nbsp;component&nbsp;(<i>S<sub>i</sub></i>)           | Pairs mapping to different transcripts       | Transcript fragmentation                                          |

Transcript quality scores S<sub><i>t</i></sub> are calculated as products of the described score components, equally weighting detected assembly errors. Assembly score <i>S</i> is computed as the mean of individual transcript scores.

In addition to transcript scores, CATS-rf provides a comprehensive set of assembly metrics, including transcript length and composition statistics, read mapping rates, positional coverage and accuracy profiles, and pair mapping consistency metrics.

CATS-rf consistently displays stronger performance than currently existing reference-free transcriptome assembly evaluation tools. Details and benchmarks can be found in the CATS [preprint](test)

# Installation 

## Installing from source

### Linux and Windows
CATS-rf consists of Bash and R scripts distributed via this repository. For the best compatibility and performance, we recommend running CATS-rf on:
- Any modern Linux distribution (e.g. Ubuntu, Debian, Fedora, etc.)
- WSL (i.e. Ubuntu on Windows)

### MacOS
If you’re using MacOS, Bash (version >= 4.0.) and GNU versions of core utilities are required. In this case, `PATH` variable should be adjusted so that CATS-rf uses the GNU versions of core utilities:

- Install **Bash ≥ 4.0** via [Homebrew](https://formulae.brew.sh/formula/bash):
```bash
brew install bash
```
- Install GNU core utilities:
```bash
brew install coreutils findutils gnu-sed gawk grep
```
- Add GNU tools to your `PATH` (adjust path depending on your architecture):
```bash
export PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/findutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/grep/libexec/gnubin:$PATH"
```
- Run the script using the installed Bash version (adjust path depending on your architecture):
```bash
/opt/homebrew/bin/bash CATS_rf
```

### Dependencies
All CATS-rf scripts must be made executable and included in the `PATH` environment variable.

The following dependencies are required and must also be included in `PATH`:

| **Dependency**      | **Tested Version** | **Homepage**                                   | **Conda Installation**                    |
|---------------------|--------------------|------------------------------------------------|-------------------------------------------|
| R                   | 4.4.2              | https://www.r-project.org/                     | `conda install conda-forge::r-base`       |
| data.table (R)      | 1.16.4             | https://cran.r-project.org/package=data.table  | `conda install conda-forge::r-data.table` |
| Bowtie2             | 2.5.4              | https://github.com/BenLangmead/bowtie2         | `conda install -c bioconda bowtie2`       |
| Samtools            | 1.21               | https://www.htslib.org/                        | `conda install -c bioconda samtools`      |
| kallisto            | 0.50.1             | https://github.com/pachterlab/kallisto         | `conda install -c bioconda kallisto`      |
| GNU Parallel        | 20220922           | https://www.gnu.org/software/parallel/         | `conda install conda-forge::parallel`     |
| bamToBed (bedtools) | 2.31.1             | https://github.com/arq5x/bedtools2             | `conda install -c bioconda bedtools`      |
| pysamstats          | 1.1.2              | https://github.com/alimanfoo/pysamstats        | `conda install -c bioconda pysamstats`    |

## Installing via conda

**CATS-rf** and its dependencies are available via [Bioconda](https://bioconda.github.io/):

```bash
conda install -c bioconda cats_rf
```

# Example usage 

CATS-rf requires a transcriptome assembly in FASTA format, along with short RNA-seq reads used during assembly in either FASTQ or FASTA format. Compressed (.gz) read files are supported.

CATS-rf supports both paired-end and single-end library configurations.

Example paired-end mode usage:

```bash
CATS_rf [OPTIONS] TRANSCRIPTOME READS1 READS2
```

Example single-end mode usage:

```bash
CATS_rf -C se -m 150 -s 100 [OTHER_OPTIONS] TRANSCRIPTOME READS1
```

Single-end mode requires three parameters to be specified: `C` for library configuration, `m` for mean fragment size, and `s` for standard deviation of fragment size. Note that single-end runs will output only general assembly statistics, read mapping metrics, and positional coverage/accuracy analysis.

# Detailed options

CATS-rf offers a comprehensive list of options which allow users to control the analysis parameters.

## Library options

`-C`: Paired- vs. signle-end library configuratiom: pe = paired-end, se = single-end, default: pe

`-S`: Library strandness, fr = forward-reverse, rf = reverse-forward, u = unstranded, a = automatic detection, default: u

CATS-rf can leverage strandness information when quantifying transcripts and calculating the local fidelity score component. When the automatic detection option is enabled, strandness is estimated using the first 100 000 read mappings.
While CATS-rf is primarily developed for Illumina reads, the analysis can be run on reads generated by other technologies. In such scenario, `S` should be accordingly adjusted.

`-Q`: Phred quality encoding of FASTQ files, 33 = phred33, 64 = phred64, default: 33

## Read mapping, transcript quantification, and read assignment options

`-R`: Random seed for read mapping, transcript quantification, and read assignment, default: 12345

Random seed is defined to ensure reproducible CATS-rf runs.

`-N`: Maximum number of distinct mappings per read, default: 10

Parameter `N` should be increased for complex transcriptome assemblies that contain a large number of isoforms, and decreased for simpler assemblies with fewer isoforms to maximize performance and accuracy. Note that mapping parameters are optimized to detect transcript sequence inaccuracies, while minimizing the number of false-positive mappings. Furthermore, secondary mappings of each read are filtered based on edit distance.

`-m:`: Estimated mean of fragment length needed for transcript quantification (single-end mode only)

`-s`: Estimated st. dev. of fragment length needed for transcript quantification (single-end mode only)

Fragment length distribution parameters `m` and `s` are required in single-end mode for transcript quantification by kallisto.

## Coverage analysis options

`-i`: Per-base coverage distribution breakpoints (specified with x,y,z...), default: "0,5,10,20,40,60,80,100"

Per-base coverage is split into intervals defined by `i` (e.g. [0-5>, [5-10>...). This category variable is used for plotting by the CATS_rf_compare script.

`-p`: Per-transcript proportion of covered bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of covered bases is split into intervals defined by `p` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-r`: Mean transcript coverage distribution breakpoints (specified with x,y,z...), default: "0,5,10,20,40,60,80,100"

Mean transcript coverage is split into intervals defined by `r` (e.g. [0-5>, [5-10>...). This category variable is used for plotting by the CATS_rf_compare script.

`-l`: Proportion of transcript length for positional relative coverage distribution analysis, default: 0.01

Transcripts are split into percentiles of size `l` for positional relative coverage distribution analysis. Coverage is expressed relative to the base with the highest coverage within the same trasncript.

`-n`: Proportion of transcript length for transcript end definition when calculating mean transcript end coverage, default: 0.02

Parameter `n` controls the relative size of transcript end regions when calculating mean transcript end coverage.

`-k`: Rolling window size for local coverage calculation (in bp) when defining low-coverage regions (LCR), default: 10 bp

`-z`: Local coverage threshold for LCR characterization, default: 3

LCRs are defined as rolling windows of size `k` with mean coverage lower than or equal to `z`.

`-u`: Per-transcript proportion of LCR bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of LCR bases is split into intervals defined by `u` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-w`: Coverage penalty weight, default: 1.5

`-e`: LCR extension penalty, default: 0.5

Parameters `w` and `e` adjust coverage penalties assigned to LCRs. Lower values of `w` and higher values of `e` increase the relative impact of LCR length on coverage penalty.

## Accuracy analysis options

`-I`: Per-base accuracy distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Accuracy is defined as the proportion of aligned read bases matching the transcript base. Per-base accuracy is split into intervals defined by `I` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-A`: Minimum accuracy for a base to be considered accurate, default: 0.95

`-P`: Per-transcript proportion of accurate bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Per-transcript proportion of accurate bases (bases with accuracy higher or equal to `A`) is split into intervals defined by `P` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-L`: Proportion of transcript length for positional accuracy distribution analysis, default: 0.01

Transcripts are split into percentiles of size `L` for positional accuracy distribution analysis.

`-K`: Rolling window size for local accuracy calculation (in bp) when defining low-accuracy regions (LAR), default: 10 bp

`-Z`: Local accuracy threshold for LAR characterization, default: 0.98

LARs are defined as rolling windows of size `K` with mean accuracy lower than or equal to `Z`.

`-U`: Per-transcript proportion of LAR bases distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1"

Per-transcript proportion of LAR bases is split into intervals defined by `U` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-E`: LAR extension penalty, default: 0.1

Extension penalty `E` adjusts accuracy penalties assigned to LARs. Higher values of `E` increases the relative impact of LAR length on accuracy penalty.

## Paired-end analysis options

These options should only be supplied in paired-end mode.

`-d`: Maximum distance from transcript ends for reads to be considered evidence of transcript end incompleteness or fragmentation (in bp), default: 40 bp

Reads mapping to transcript ends are considered evidence for transcript end incompleteness or fragmentation. Parameter `d` controls the relative size of transcript end regions when identifying such reads.

`-x`: Multiplicative factor for lower distance outlier threshold calculation, default: 8

`-X`: Multiplicative factor for higher distance outlier threshold calculation, default: 10

`-c`: Correction factor for distance outlier threshold calculation, default: 5

Parameters `x`, `X`, and `c` adjust distance penalty calculation. Read pairs are classified as mapping too far apart if their distance exceeds the lower distance threshold, defined as D<sub>1</sub> = Q<sub>3</sub>(d) + x × (IQR(d) + c). These reads are assigned a distance penalty P<sub>d</sub> = d / D<sub>2</sub>, where D<sub>2</sub> = Q<sub>3</sub>(d) + X × (IQR(d) + c), with the penalty capped at 1. Parameter `x` adjusts the classification of read pairs mapping too far apart, while  `X` controls distance penalty scailing. Parameter `c` ensures penalty stability in libraries with a high proportion of overlapping read pairs.

`-y`: Per-transcript proportion of improperly paired reads within a transcript distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Improperly paired reads include reads with pairs not mapped to the assembly, reads with pairs mapping to an inappropriate strand or in an unexpected orientation and read pairs mapping too far apart. Per-transcript proportion of improperly paired reads within a transcript is split into intervals defined by `y` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-f`: Minimum number of bridging events for transcripts to be considered fragmented, default: 3

A transcript is considered fragmented if more than `f` reads representing bridging events map to transcript end regions.

`-F`: Per-transcript proportion of reads whose pairs map to other transcripts distribution breakpoints (specified with x,y,z...), default: "0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,1"

Per-transcript proportion of reads whose pairs map to other transcripts is split into intervals defined by `F` (e.g. [0-0.2>, [0.2-0.4>...). This category variable is used for plotting by the CATS_rf_compare script.

`-a`: Alpha compression factor for the beta-like transformation applied to bridge index during integrity score component calculation, default: 7

`-b`: Beta compression factor for the beta-like transformation applied to bridge index during integrity score component calculation, default: 0.5

Bridge index measures the proportion of reads whose pair maps to a different transcript and considers the mapping distance of each read from the end of its respective transcript. This definition gives more weight to bridging events near transcript ends. Integrity score component is calculated using a beta-like transformation of bridge index. Compression factors `a` and `b` control the shape of the transformation: higher `a` increases sensitivity to fragmentation, while higher `b` reduces the likelihood of false-positive fragmentation penalties in transcripts with minimal bridging evidence.

## General options

`-t`: Number of CPU threads, default: 10

Several steps of CATS-rf pipeline support parallel processing. This includes read mapping, transcript quantification, read assignment, SAM file processing, positional coverage and accuracy calculation and analysis, as well as positional paired-end analysis. Recommended number of threads: 10-20.

`-G`: Percentage of available RAM used by GNU sort, default: 50

CATS-rf utilizes GNU sort in several steps of the pipeline. Higher values of `G` will ensure faster sorting, but may exhaust available RAM. In such scenarios, CATS-rf will resort to sorting with mimimum RAM usage.

`-M`: Memory block size for GNU Parallel, default: 512M

Parameter `M` controls the input block size used by GNU Parallel when splitting the mapping table for read assignment. If sufficient RAM is available, increasing `M` is recommended to minimize artifacts introduced by file splitting.

`-T`: Number of splits performed on positional and read pair mapping tables, default: 3

Positional and read pair mapping tables are split before analysis to reduce RAM usage. Increase `T` when working with limited memory to further reduce RAM demands.

`-D`: Output directory name, default: TRANSCRIPTOME_CATS_rf_dir

`-o`: Output file prefix, default: TRANSCRIPTOME_CATS_rf

`-O`: Overwrite the results directory, default: off

`-h`: Show usage information

Usage is printed when using the `h` option or calling CATS-rf without arguments.

# Output explanation

Summary files represent the main output of CATS-rf. In paired-end mode, four summary files are produced:

`assembly_score_summary.tsv` contains descriptive statistics of transcript score components and the overall assembly score. The contents of this file are printed when CATS-rf finishes running in paired-end mode.

`general_statistics_table.tsv` contains descriptive statistics of transcript length (mean, median, interquartile range, range, N50, L50, N90, L90), GC content, and read mapping rate.

`coverage_and_accuracy_analysis_summary.tsv` contains summarized results of coverage and accuracy analysis. The contents of this file are printed when CATS-rf finishes running in single-end mode.

`paired_end_analysis_summmary.tsv` contains summarized results of paired-end analysis, including local fidelity analysis and integrity analysis.

CATS-rf also produces several .tsv files containing detailed per-transcript metrics.

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
| `lcr_base_prop_cat`           | Proportion of bases in low-coverage regions category |
| `coverage_score_component`    | Coverage score component                             |

`per_base_coverage_distribution.tsv` contains distribution of per-base coverage in the assembly.

`relative_coverage_median_by_transcript_position.tsv` contains median values of mean relative coverage per transcript position percentile.

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

`per_base_accuracy_distribution.tsv` contains distribution of per-base accuracy in the assembly.

`accuracy_median_by_transcript_position.tsv` contains median values of mean accuracy per transcript position percentile.

`lar_list.tsv` contains low-accuracy region coordinates.

## Local fidelity analysis

`local_fidelity_stats.tsv` contains local fidelity analysis results for each transcript:

| **Column**                                 | **Description**                                                                                |
|--------------------------------------------|------------------------------------------------------------------------------------------------|
| `transcript`                               | Transcript name                                                                                |
| `unmapped_pair_read_N`                     | Number of reads with pairs not mapped to the assembly                                          |
| `unmapped_pair_read_prop`                  | Proportion of reads with pairs not mapped to the assembly                                      |
| `unmapped_pair_tr_end_read_N`              | Number of reads with pairs not mapped to the assembly on transcript ends                       |
| `unmapped_pair_tr_end_read_prop`           | Proportion of reads with pairs not mapped to the assembly on transcript ends                   |
| `improp_pair_strand_orientation_read_N`    | Number of reads whose pairs map to an inappropriate strand or in an unexpected orientation     |
| `improp_pair_strand_orientation_read_prop` | Proportion of reads whose pairs map to an inappropriate strand or in an unexpected orientation |
| `improp_pair_distance_read_N`              | Number of reads whose pairs map too far apart                                                  |
| `norm_distance_penalty`                    | Normalized transcript distance penalty                                                         |
| `improp_pair_within_tr_read_N`             | Number of improperly paired reads within a transcript                                          |
| `improp_pair_within_tr_read_prop`          | Proportion of improperly paired reads within a transcript                                      |
| `improp_pair_within_tr_read_prop_category` | Proportion of improperly paired reads within a transcript category                             |
| `local_fidelity_score_component`           | Local fidelity score component                                                                 |

`read_pairs_mapping_to_inconsistent_strand_or_orientation.tsv` contains coordinates of read pairs mapping to an inappropriate strand or in an unexpected orientation.

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

# Citation

CATS is an academic software distributed under the MIT license. if you use CATS, please cite the CATS preprint:

(add reference)

# Troubleshooting

Please report all potential bugs in the Issues tracker.

# Changelog

Version 1.0.0. Initial commit, June 3, 2025.
