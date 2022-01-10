---
layout: post
title: "Operation on BigWig Files"
date: 2022-01-07
category: tutorial
tags: [BigWig]
---

Some operations on BigWig files.

<!--more-->

### Extract/subset bigwig file for a given genomic region

This is a solution in python version (3.0+) using a package called pyBigWig to extract a given genomic region from a whole genome bigwig file.

Prepare your input bigwig file:

```python
import pyBigWig

# First open bigwig file
bwIn = pyBigWig.open('input.bw')

# check bigwig file header
print(bwIn.header())
```

Prepare output, since your output doesn’t have a header, you need to add the header using the chosen chromosome size, here I’m using a region from chr18 as an example.

```python
bwOutput = pyBigWig.open('output.bw','w')
bwOutput.addHeader([('chr18',80373285)]) # chromosome size

for x in bwIn.intervals('chr18',62926563,63516911):
    bwOutput.addEntries(['chr18'],[x[0]],ends=[x[1]],values=[x[2]])

bwOutput.close()
```

### Merge bigwig files using average value

Software preparation

```bash
# install fetchChromSizes
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
# install wigToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
# install wiggletools
conda install -c bioconda wiggletools
```

```bash
fetchChromSizes hg38 > hg38.genome
wiggletools mean sample1.bw sample2.bw ... sampleN.bw | wigToBigWig stdin hg38.genome mean.bw
```

### References
* [extract bigwig regions](https://bioinfocore.com/blogs/extract-subset-bigwig-file-for-a-given-genomic-region/)
* [bigwig average](https://github.com/deeptools/deepTools/issues/723)
