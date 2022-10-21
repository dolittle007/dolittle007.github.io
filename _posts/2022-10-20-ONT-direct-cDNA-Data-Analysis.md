---
layout: post
title: "Nanopore direct cDNA data analysis"
date: 2022-10-20
category: tutorial
tags: [Nanopore, ONT, Long-reads, analysis]
---

An Introduction to Nanopore direct cDNA data analysis.

<!--more-->

### Software preparation

```bash
# Install Guppy CPU version
wget -c https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.3.8_linux64.tar.gz
tar zxvf ont-guppy-cpu_6.3.8_linux64.tar.gz

# Install Guppy GPU version
wget -c https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_6.3.8_linux64.tar.gz
tar zxvf ont-guppy_6.3.8_linux64.tar.gz

# add ont-guppy-cpu/bin to $PATH in .bashrc file
PATH=/path/to/ont-guppy-cpu/bin:$PATH

# install minimap2 and samtools
conda install -c bioconda minimap2 # paftools.js will be install automatically.
conda install -c bioconda samtools

# install pychopper and nhmmscan
conda install -c epi2melabs -c conda-forge -c bioconda "epi2melabs::pychopper"
```

### Step1: Basecalling

#### CPU-based basecalling
```bash
guppy_basecaller --input_path ./fast5 --save_path ./guppy_output --flowcell FLO-MIN106 --kit SQK-RNA002 --calib_detect --num_callers 16 --cpu_threads_per_caller 8 --compress_fastq --trim_strategy none
```
#### GPU-based basecalling
```bash
guppy_basecaller --input_path ./fast5 --save_path ./guppy_output --flowcell FLO-MIN106 --kit SQK-RNA002 --calib_detect --num_callers 16 ----gpu_runners_per_device 80 -x "cuda:all" --compress_fastq --trim_strategy none
```
### Step2: Identify full-length Nanopore cDNA reads
In this step, Pychopper v2 is used to identify, orient and trim full-length Nanopore cDNA reads. Pychopper can also rescue fused reads ([__chimeric reads__](https://yulijia.net/en/bioinformatics/2015/12/21/Linear-Chimeric-Supplementary-Primary-and-Secondary-Alignments.html)).
