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
![center](/figures/2016-09-13-RNA-seq-analysis/rna_seq_workflow.png) 

### Create Mapping Indices
Before we can perform NGS read mapping, we will create the genome indices using the genome FASTA file as input. You can re-use these indices in all your future short read mapping. However, if you wish to map to a different genome build/assembly, you have to re-run this step using different genome sequences and save the indices in a different directory.

#### Usage
{% highlight bash %}
STAR --runMode genomeGenerate --genomeDir path_to_genomedir --genomeFastaFiles reference_fasta_file(s)
{% endhighlight %}
#### Execute
{% highlight bash %}
mkdir GENOME_data/star
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir GENOME_data/star \
--genomeFastaFiles GENOME_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
{% endhighlight %} 
#### Options
{% highlight bash %}
--runThreadNdefines  # the number of threads to be used for genome generation.
--runMode genomeGenerate # directs STAR to run genome indices generation job.
--genomeDir # path to the directory where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to writing permissions. The file system needs to have at least 100GB of disk space available for a typical mammalian genome.
--genomeFastaFiles # one or more FASTA files with the genome reference sequences.
{% endhighlight %}

