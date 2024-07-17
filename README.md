HuGenVarDetective
===============
HuGenVarDetective (GVD) contains code and workflow for analysing human genes using the Genome Analysis Toolkit (GATK) version 4. The goal of this project was divided into two parts: [i] to leverage GATK4's powerful tools for variant discovery (SNP/MNP and INDELs) to conduct comprehensive genomic analyses and [ii] to annotate the discovered variants within our gene of interest using snpeff and snipsift.

Content
=======
* [Aim 1](#aim-1)
   * [About GATK4](#about-gatk4)
   * [Installation and Usage](#installation-and-usage)
   * [Getting Started (Project)](#getting-started-project)
   * [Pre-requisites before Analysis](#pre-requisites-before-analysis)
      * [Directory Structure and Configuration Data](#directory-structure-and-configuration-data)
      * [Input Data](#input-data)
      * [Supporting Files](#supporting-files)
      * [The Script](#the-script)
* [Aim 2](#aim-2)
   * [About SnpEff &amp; SnpSift](#about-snpeff--snpsift)

Aim 1
=====
## About GATK4
The Genome Analysis Toolkit (GATK) is a software package developed by the Broad Institute that focuses on variant discovery and genotyping. It aims to bring together well-established tools from the [GATK](http://www.broadinstitute.org/gatk) and [Picard](http://broadinstitute.github.io/picard/) codebases under a streamlined framework, and to enable selected tools to be run in a massively parallel way on local clusters or in the cloud using [Apache Spark](http://spark.apache.org/). The toolkit is widely used in the field of bioinformatics, especially for analysing high-throughput sequencing data.

## Installation and Usage
Documentation on the requirements, installation/building and running/usage of the software package can be found at [Broad Institute's GATK page](https://github.com/broadinstitute/gatk?tab=readme-ov-file#requirements)

## Getting Started (Project)
The bash script for this toolkit follows the author's best practice workflows. This project was conducted using both nanopore and illumina sequenced data.

## Pre-requisites before Analysis
A copy of the original directory containing reads (illumina/nanopore) is used. For this analysis, quality checks on sequenced reads and removal of adapters, poor quality and primers are recommended using QC tools (such as FastQC, fastp, trimmomatic, porechop, NanoPlot, MinionQC etc)

### Directory Structure and Configuration Data
The overall directory structure for this GATK analysis is illustrated below:

```
../
    experiment_name/
        aligned_data/
        data/
        reads/
        results/
        supporting_files/
```

Within the copied experiment folder `experiment_name`, the above directories are created to facilitate the pipeline analysis process.

### Input Data
For Nanopore sequenced data, the file structure of the input directory set up for this GATK analysis looks like:

```
../
    experiment_name/
        reads/
            A1.fastq.gz
            A2.fastq.gz
            A3.fastq.gz
            A4.fastq.gz
```

For Illumina sequenced data, the file structure of the input directory set up for this GATK analysis looks like:

```
../
    experiment_name/
        reads/
            A1_R1.fastq.gz
            A1_R2.fastq.gz
            A2_R1.fastq.gz
            A2_R2.fastq.gz
            A3_R1.fastq.gz
            A3_R2.fastq.gz
            A4_R1.fastq.gz
            A4_R2.fastq.gz
```

Note: The reads used for the analysis must be quality checked and trimmed to get accurate and high quality reads into the pipeline for analysis. The trimmed files may be labelled as `.trimmed.fastq.gz` for nanopore and `R1.trimmed.fastq.gz` for illumina data for easy identification of data.

### Supporting Files
Supporting files are one of the requirements for this project. They include the <reference.fasta> file, a sequence dictionary from the reference.fasta file, and <reference.fasta> indices produced for alignment and variant calling. These are placed in the the `supporting_files` directory and an example illustration is displayed below.

```
mmp3.nc_000011.10.reference.dict
mmp3.nc_000011.10.reference.fasta
mmp3.nc_000011.10.reference.fasta.amb
mmp3.nc_000011.10.reference.fasta.ann
mmp3.nc_000011.10.reference.fasta.bwt
mmp3.nc_000011.10.reference.fasta.fai
mmp3.nc_000011.10.reference.fasta.pac
mmp3.nc_000011.10.reference.fasta.sa
```

These files were generated from the MMP3 reference file `mmp3.nc_000011.10.reference.fasta`.

### The Script
Check the bash script and make specific changes to it (eg. the pathnames for the variables that are specific to the user)

Aim 2
=====
## About SnpEff & SnpSift

