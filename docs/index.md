# CATS-rf

<img src="cats_rf_logo.png" alt="Logo" width="750" height="160"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](../LICENSE)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-brightgreen)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/cats-rf.svg)](https://anaconda.org/bioconda/cats-rf)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101/2025.07.22.666112v1-orange)](https://www.biorxiv.org/content/10.1101/2025.07.22.666112v1)

## Documentation

- [Introduction](#introduction)
- [Installation](installation.md)
- [Test data](test-data.md)
- [Example usage](usage.md)
- [Detailed options](options.md)
- [Output explanation](output.md)
- [Assembly comparison (`CATS_rf_compare`)](compare.md)
- [Citation](citation.md)
- [Troubleshooting](troubleshooting.md)
- [Changelog](changelog.md)

# Introduction

CATS-rf is the reference-free module of the CATS (Comprehensive Assessment of Transcript Sequences) framework. It evaluates the quality of transcriptomes assembled de novo from short reads, relying solely on RNA-seq reads used in the assembly construction. The pipeline maps reads back to the assembled transcripts and examines mapping evidence suggesting misassembly. Quality evaluation is performed at the transcript level, integrating four score components each targeting specific assembly errors:

| **Score Component**                                            | **Evidence**                                 | **Targeted Assembly Errors**                                      |
|----------------------------------------------------------------|----------------------------------------------|-------------------------------------------------------------------|
| Coverage&nbsp;component&nbsp;(<i>S<sub>c</sub></i>)            | Low-coverage regions                         | Insertions, redundancy                                            |
| Accuracy&nbsp;component&nbsp;(<i>S<sub>a</sub></i>)            | Low-accuracy regions                         | Sequence inaccuracy                                               |
| Local&nbsp;fidelity&nbsp;component&nbsp;(<i>S<sub>l</sub></i>) | Inconsistent pair mapping within transcripts | Structural errors (e.g. deletions, translocations, inversions...) |
| Integrity&nbsp;component&nbsp;(<i>S<sub>i</sub></i>)           | Pairs mapping to different transcripts       | Transcript fragmentation                                          |

Transcript quality score S<sub><i>t</i></sub> is calculated as the product of the described score components, equally weighting detected assembly errors. Assembly score <i>S</i> is computed as the mean of individual transcript scores. All components and scores are normalized to a range between 0 and 1, where higher values indicate better quality.

In addition to transcript scores, CATS-rf provides a comprehensive set of assembly metrics, including transcript length and composition statistics, read mapping rates, positional coverage and accuracy profiles, and pair mapping consistency metrics.

CATS-rf consistently displays stronger performance than currently existing reference-free transcriptome assembly evaluation tools. For detailed benchmarks and methodology, please refer to the CATS [preprint](https://www.biorxiv.org/content/10.1101/2025.07.22.666112v1).
