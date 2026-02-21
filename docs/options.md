# Detailed options

CATS-rf offers a comprehensive list of options which allow users to control the analysis parameters.

## Library type options

`-C`: Paired- vs. single-end library configuration: pe = paired-end, se = single-end, default: pe

`-S`: Library strandness, fr = forward-reverse, rf = reverse-forward, u = unstranded, a = automatic detection, default: u

CATS-rf can leverage strandness information when quantifying transcripts and calculating the local fidelity score component. When the automatic detection option is enabled, strandness is estimated using the first 100 000 read mappings.

While CATS-rf was primarily tested on Illumina data, the analysis can be run on assemblies generated from other short-read platforms. In such scenario, `S` should be adjusted accordingly. If the strandness of the data is unknown, it is recommended to use either unstranded mode or automatic detection. Note that in unstranded mode, read pairs are expected to map to opposite strands. This is consistent with the behavior of virtually every short-read sequencing technology.

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