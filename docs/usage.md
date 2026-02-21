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