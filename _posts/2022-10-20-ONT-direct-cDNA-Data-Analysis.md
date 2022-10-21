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

# install pychopper
conda install -c epi2melabs -c conda-forge -c bioconda "epi2melabs::pychopper" # nhmmscan will be install automatically.
```

### Annotation preparation
```bash
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gff3.gz
gunzip gencode.v39.annotation.gff3.gz

paftools.js gff2bed gencode.v39.annotation.gff3 > hg38.bigbed
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

:warning: `Do NOT Turn on trimming (setting *--trim_strategy none* ) during basecalling as it will remove the primers needed for classifying the reads!`

### Step2: Identify full-length Nanopore cDNA reads
In this step, [__Pychopper__](https://github.com/epi2me-labs/pychopper) is used to identify, orient and trim full-length Nanopore cDNA reads. Pychopper can also rescue fused reads ([__chimeric reads__](https://yulijia.net/en/bioinformatics/2015/12/21/Linear-Chimeric-Supplementary-Primary-and-Secondary-Alignments.html)).

![center](/figures/2022-10-20-ONT-direct-cDNA-Data-Analysis/SSP_VNP.png) 

| Term | Description |
| ------ | ----------- |
| __SSP__ | strand-switching primer |
| __VNP__ | anchored oligo(dT) VN primer |

#### Combine called FASTQ files
```bash
cat ./guppy_output/pass/*.gz > raw.fastq.gz
```

#### First round full-length cDNA reads identification with standard parameters using the default pHMM backend and autotuned cutoff parameters estimated from subsampled data:
```bash
pychopper -r report.pdf -k PCS109 -u unclassified.fq -w rescued.fq raw.fastq.gz full_length_output.fq
```
#### Second round full-length cDNA reads identification applied to the unclassified direct cDNA reads with DCS-specific read rescue enabled (parameter -x).

DCS109 can suffer from a specific reverse transcription artefact, which will lead to 2D-like reads with two VNP primers at the ends.

The -x DCS109 mode should be used on the unclassified reads only (while keeping the classified ones of course). This will classify the -VNP,VNP configurations as full length (also rescuing fused reads with this configuration, though they are not prevalent).

```bash
pychopper -r report_2.pdf -k PCS109 -x PCS109 -u unclassified_2.fq -w rescued_2.fq unclassified.fq full_length_output_2.fq
```
#### Full-length and rescued reads are merged and used for subsequent steps.
```bash
cat full_length_output.fq full_length_output_2.fq rescued.fq rescued_2.fq > full_length_cdna.fastq
```

### Step3: Aign to Genome
We currently recommend using [__minimap2__](https://github.com/lh3/minimap2) to align to the reference genome.

```bash
minimap2 -Y -t 8 -R "@RG\tID:Sample\tSM:hs\tLB:ga\tPL:ONT" --MD -ax splice -uf -k14 --junc-bed hg38.bigbed hg38.fasta full_length_cdna.fastq > aligned.sam
samtools sort -@ 8 -O BAM align.sam -o aligned.sort.bam
samtools index aligned.sort.bam
```

### References
* [Pychopper manual](https://github.com/epi2me-labs/pychopper)
* [Nanopore sequencing of RNA and cDNA molecules in Escherichia coli](https://rnajournal.cshlp.org/content/28/3/400.full)
* [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html)
* [Pychopper update (v2) in ONT community](https://community.nanoporetech.com/posts/pychopper-update-v2)
