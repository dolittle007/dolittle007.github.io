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
```
