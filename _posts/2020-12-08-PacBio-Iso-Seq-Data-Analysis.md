---
layout: post
title: "PacBio Iso-Seq data analysis"
date: 2020-12-08
category: tutorial
tags: [PacBio, Long-reads, analysis]
---

An Introduction to PacBio Iso-Seq data analysis.

<!--more-->

### First Post

PacBio Sequel Systems (Sequel I and Sequel II) have two modes, Continuous Long Reas (CLR) and Circular Consensus Sequencing (CCS).
![center](/figures/2020-12-08-PacBio-Iso-Seq-Data-Analysis/CLR_CCS.png) 


### Software preparation


```bash
conda install -c bioconda pbccs # install CCS version 5.0 or above
conda install -c bioconda lima
conda install -c bioconda isoseq3 # install isoseq v3.4.0 or above
```
#### Generate CCS
If you don't already have CCS reads, run

```bash
ccs [movie].subreads.bam [movie].ccs.bam 
```
Note that _ccs_ run polish by dedault unless you don't want it ( --skip-polish ), so __isoseq3 polish__ is no longer needed.

A typical command to run ccs will be like this, as followed.
```bash
ccs [movie].subreadset.xml [movie].consensusreadset.xml --log-level INFO --report-json [movie].report.json --hifi-summary-json [movie].hifi_summary.json --log-file [movie].ccs.log --report-file [movie].report.txt --metrics-json [movie].zmw_metrics.json.gz -j 64
```
OR
```bash
ccs [movie].subreadset.bam [movie].ccs.bam --log-level INFO --report-json [movie].report.json --hifi-summary-json [movie].hifi_summary.json --log-file [movie].ccs.log --report-file [movie].report.txt --metrics-json [movie].zmw_metrics.json.gz -j 64
```
One important changes for _ccs_ (>=v5.0.0) is that it has the --all mode. In this mode, _ccs_ outputs one representative sequence per productive ZMW, irrespective of quality and passes. 
Note that _ccs_ is now running on the Sequel IIe instrument, transferring HiFi reads directly off the instrument.
The on-instrument ccs version and also SMRT Link ≥v10 run in the --all mode by default. 
But don't worry, if you want to only use HiFi reads, _ccs_ automatically generates additional files for your convenience that only contain HiFi reads:
* hifi_reads.fastq.gz
* hifi_reads.fasta.gz
* hifi_reads.bam
