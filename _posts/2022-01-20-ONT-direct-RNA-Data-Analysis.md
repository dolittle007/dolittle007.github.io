---
layout: post
title: "Nanopore direct RNA data analysis"
date: 2022-01-20
category: tutorial
tags: [Nanopore, ONT, Long-reads, analysis]
---

An Introduction to Nanopore direct RNA data analysis.

<!--more-->

### Software preparation

```bash
# Install Guppy CPU version
wget -c https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.5.7_linux64.tar.gz
tar zxvf ont-guppy-cpu_6.5.7_linux64.tar.gz

# Install Guppy GPU version
wget -c https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_6.5.7_linux64.tar.gz
tar zxvf ont-guppy_6.5.7_linux64.tar.gz

# add ont-guppy-cpu/bin to $PATH in .bashrc file
PATH=/path/to/ont-guppy-cpu/bin:$PATH

# install minimap2 and samtools
conda install -c bioconda minimap2 # paftools.js will be install automatically.
conda install -c bioconda samtools
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
guppy_basecaller --input_path ./fast5 --save_path ./guppy_output --flowcell FLO-MIN106 --kit SQK-RNA002 --calib_detect --num_callers 16 --cpu_threads_per_caller 8 --compress_fastq --reverse_sequence --u_substitution
```

#### GPU-based basecalling
```bash
guppy_basecaller --input_path ./fast5 --save_path ./guppy_output --flowcell FLO-MIN106 --kit SQK-RNA002 --calib_detect --num_callers 16 ----gpu_runners_per_device 80 -x "cuda:all" --compress_fastq --reverse_sequence --u_substitution
```
WARNING: Use RNA-specific parameters, --calib_detect, --reverse_sequence, --u_substitution.

#### Options 
```bash
--input_path  # The location of FAST5 files
--save_path # The location of output FASTQ files. It have three subfolders (pass, fail, and calibration_strands).
--calib_detect  # Enable RNA calibration strand (RCS) detection and filtering.
--reverse_sequence # Reverse the called sequence.
--u_substitution # Substitute 'U' for 'T' in the called sequence.
--compress_fastq # Compress fastq output files with gzip
--flowcell # flowcell name
--kit # kit name
```

List supported flowcells and kits:
```bash
guppy_basecaller --print_workflows
```
Alternatively, you can specific config file

```bash
guppy_basecaller --input_path ./fast5 --save_path ./guppy_output -c rna_r9.4.1_70bps_hac --calib_detect --num_callers 16 --cpu_threads_per_caller 8 --compress_fastq --reverse_sequence and --u_substitution
```

#### What is RNA Calibration Strand (RCS)?
The RNA CS (RCS) is the RNA Calibration Strand is the Enolase II from [__YHR174W__](http://useast.ensembl.org/Saccharomyces_cerevisiae/Gene/Summary?g=YHR174W;r=VIII:451327-452640;t=YHR174W_mRNA) at a concentration of 50 ng/μL. The reference fasta file for YHR174W ENO2 is available at __ont-guppy-cpu/data/YHR174W.fasta__.
RCS is included in included in the Direct RNA Sequencing kit, SQK-RNA002, and PCR-cDNA Barcoding Kit, SQK-PCB109

### Step3: Aign to Genome
We currently recommend using [__minimap2__](https://github.com/lh3/minimap2) to align to the reference genome.

```bash
minimap2 -Y -t 8 -R "@RG\tID:Sample\tSM:hs\tLB:ga\tPL:ONT" --MD -ax splice -uf -k14 --junc-bed hg38.bigbed hg38.fasta sample.fastq | samtools sort -@ 8 -O BAM -o aligned.bam -
samtools index aligned.bam
```


