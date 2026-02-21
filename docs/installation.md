# Installation 

## Running CATS-rf via Docker

Starting from version 1.0.2, CATS-rf can be run using Docker:

Build the Docker image with:

```bash
docker pull bodulic/cats-rf
```

Run the container with:

```bash
docker run --rm -v "$PWD":/data -w /data bodulic/cats-rf CATS_rf

docker run --rm -v "$PWD":/data -w /data bodulic/cats-rf CATS_rf_compare
```

## Running CATS-rf via Singularity

Build the Singularity image with:

```bash
singularity build cats_rf.sif docker://bodulic/cats-rf:latest
```

Run the container with:

```bash
singularity run cats_rf.sif CATS_rf

singularity run cats_rf.sif CATS_rf_compare
```

## Installing CATS-rf via conda

CATS-rf and its dependencies can be directly installed via [Bioconda](https://anaconda.org/bioconda/cats-rf):

```bash
conda install -c bioconda cats-rf
```

In case of dependency conflicts, please see the Troubleshooting section.

## Installing CATS-rf from source

CATS-rf consists of Bash and R scripts located in the `scripts` directory of this repository. After cloning the repository, all CATS-rf scripts must be included in the `PATH` environment variable.
The following dependencies are required:

| **Dependency**      | **Tested Version** | **Homepage**                                   | **Conda Installation**                    |
|---------------------|--------------------|------------------------------------------------|-------------------------------------------|
| R                   | 4.3.0.-4.4.3       | https://www.r-project.org                      | `conda install conda-forge::r-base`       |
| data.table (R)      | 1.16.4             | https://cran.r-project.org/package=data.table  | `conda install conda-forge::r-data.table` |
| Bowtie2             | 2.5.4              | https://github.com/BenLangmead/bowtie2         | `conda install -c bioconda bowtie2`       |
| Samtools            | 1.21               | https://www.htslib.org                         | `conda install -c bioconda samtools`      |
| kallisto            | 0.50.1             | https://github.com/pachterlab/kallisto         | `conda install -c bioconda kallisto`      |
| GNU Parallel        | 20220922           | https://www.gnu.org/software/parallel          | `conda install conda-forge::parallel`     |
| bedtools (bamToBed) | 2.31.1             | https://github.com/arq5x/bedtools2             | `conda install -c bioconda bedtools`      |
| python              | 3.6-3.12           | https://github.com/arq5x/bedtools2             | `conda install -c bioconda bedtools`      |
| pysamstats          | 1.1.2              | https://github.com/alimanfoo/pysamstats        | `conda install -c bioconda pysamstats`    |

R, Bowtie2, Samtools, kallisto, GNU Parallel, bedtools (bamToBed), and pysamstats executables must be included in `PATH`. R package data.table can be installed via conda or directly in R with `install.packages("data.table")`

### MacOS
If you are using MacOS, Bash (version >= 4.0) and GNU versions of core utilities are required. In this case, `PATH` should be adjusted so that CATS-rf uses GNU versions of core utilities:

- Install Bash ≥ 4.0 via [Homebrew](https://formulae.brew.sh/formula/bash):

```bash
brew install bash
```

- Install GNU utilities:

```bash
brew install coreutils gnu-sed gawk
```

- Add Bash and GNU utilities to your `PATH` (adjust path depending on your architecture):

For Apple Silicon:
```bash
export PATH="/opt/homebrew/bin:$PATH"
export PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
```

For Intel-based configurations:
```bash
export PATH="/usr/local/bin:$PATH"
export PATH="/usr/local/opt/coreutils/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"
```

- Run CATS-rf using the installed Bash version:
  
```bash
bash CATS_rf
```

The stated changes can be made permanent by modifying the appropriate .rc file. 