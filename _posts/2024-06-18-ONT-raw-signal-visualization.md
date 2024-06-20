---
layout: post
title: "Nanopore raw data visualization using squigualiser"
date: 2024-06-18
category: tutorial
comments: true
tags: [Nanopore, ONT, Long-reads, visualization]
---

An Introduction to Nanopore raw data visualization using [squigualiser](https://hiruna72.github.io/squigualiser/)

<!--more-->

### Software preparation

#### Tool for converting raw data to BLOW5 format
* If the raw data is in POD5 format, [Blue-crab](https://github.com/Psy-Fer/blue-crab) is required.
* If the raw data is in FAST5 format, [slow5tools](https://hasindu2008.github.io/slow5tools/) is required.

```bash
# Create environment
mamba create --name ont python=3.9
mamba activate ont


# Install blue-crab
mamba install zstd
python3 -m pip install --upgrade pip
pip install blue-crab

# Install slow5tools
mamba install hdf5
mamba install slow5tools
```

#### Aligning raw signals to basecalled reads using [F5C](https://github.com/hasindu2008/f5c)

```bash
mamba install f5c=1.4

```

#### Signal-to-read visualization using [squigualiser](https://github.com/hiruna72/squigualiser/)

```bash
pip install squigualiser
```

### Step1: Converting data format

#### For multiple FAST5 files

```bash
slow5tools f2s ./fast5_dir -d blow5_dir # convert multiple FAST5 files to multiple BLOW5 files
slow5tools merge blow5_dir -o data.blow5 # merge BLOW5 into one
slow5tools get data.blow5 -l read_ids.txt --to blow5 -o target.blow5 # extract records from a blow5 file based on a list of read ids
slow5tools index target.blow5 # index BLOW5 file
```

#### For single POD5 file

```bash
blue-crab p2s data.pod5 -o data.blow5
```

### Step2: raw signals to basecalled reads alignment

```bash
f5c resquiggle -c --rna --pore r9 -o target.paf target.fastq target.blow5
```

### Step3: Signal-to-read visualization

```bash
squigualiser plot -f target.fastq -s target.blow5 -a target.paf -o out_dir --save_svg

```

#### Options

```bash
-r read_id # specify the read with read_id to plot
--rna # specify for RNA reads
```
