---
layout: post
title: "RNA-seq data analysis"
date: 2016-09-13
category: tutorial
tags: [RNA-seq, analysis, STAR, plot, R]
---

Below shows a general workflow for carrying out a RNA-Seq experiment. In this guide, I will focus on the pre-processing of NGS raw reads, mapping, quantification and identification of differentially expressed genes and transcripts.

<!--more-->

### RNA-Seq Analysis Workflow



{% highlight bash %} 
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir GENOME_data/star \
--genomeFastaFiles GENOME_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
{% endhighlight %} 
