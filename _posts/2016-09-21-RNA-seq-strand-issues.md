---
layout: post
title: "RNA-seq data analysis (strand issues)"
date: 2016-09-21
category: opinion
tags: [RNA-seq, analysis, strand, dUTP]
---

I had been working on strand-specific paired-end reads from HiSeq lately and I had trouble mapping reads back to assembled transcripts using STAR as well as using RSEM to estimate transcript abundance. 

<!--more-->
### Bowtie and Tophat flags for strand-specific reads
Tophat uses --fr-firststrand for a library created by the dUTP method. This is stated clearly in the manual, so it is easy to understand. In contrast, Bowtie/Bowtie2 uses --fr, --rf, --ff to specify the orientation of paired-end reads.

--fr means the upstream read (/1) is from a forward strand and the downstream read (/2) is from a reverse strand.

--rf means the upstream read (/1) is from a reverse strand and the downstream read (/2) is from a forward strand.

--ff means both reads are from a forward strand.
![center](/figures/2016-09-21-RNA-seq-strand-issue/pe-orient.png) 

***
[**RSEM**](http://likit.github.io/running-bowtiebowtie2-rsem-and-tophat-on-dutp-strand-specific-reads.html "RNA-Seq by Expectation-Maximization")
