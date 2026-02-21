# Assembly comparison with `CATS_rf_compare`

CATS-rf also supports direct comparison of multiple analysed assemblies. The `CATS_rf_compare` script generates summary tables and visualizations that compare the most significant CATS-rf results of each assembly. As such, CATS-rf should be run on each individual assembly and the resulting CATS-rf output directories should then act as input to `CATS_rf_compare`.

## `CATS_rf_compare` dependencies

`CATS_rf_compare` requires the following dependencies:

| **Dependency**           | **Tested Version** | **Homepage**                                     | **Conda Installation**                    | **R installation**               |
|--------------------------|--------------------|--------------------------------------------------|-------------------------------------------|----------------------------------|
| R                        | 4.3.0.-4.4.3       | https://www.r-project.org                        | `conda install conda-forge::r-base`       | /                                |
| pandoc                   | 2.19.2             | https://pandoc.org/                              | `conda install conda-forge::pandoc`       | /                                |
| rmarkdown (R)            | 2.29               | https://cran.r-project.org/package=rmarkdown     | `conda install conda-forge::r-rmarkdown`  | `install.packages("rnarkdown)`   |
| data.table (R)           | 1.16.4             | https://cran.r-project.org/package=data.table    | `conda install conda-forge::r-data.table` | `install.packages("data.table")` |
| ggplot2 (R)              | 3.5.1              | https://cran.r-project.org/web/packages/ggplot2  | `conda install conda-forge::r-ggplot2`    | `install.packages("ggplot2")`    |
| ggdist (R)               | 3.3.2              | https://cran.r-project.org/web/packages/ggdist   | `conda install conda-forge::r-ggdist`     | `install.packages("ggdist")`     | 

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

Note on transcriptome assembly order and names: Assemblies will appear in the order they were provided on the command line when running the tool. For visualization purposes, assembly names are limited to a maximum of 20 characters; names exceeding this limit will be truncated. If multiple assemblies share the same name, a numeric suffix (e.g., .1, .2, etc.) will be appended to distinguish these assemblies.

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
