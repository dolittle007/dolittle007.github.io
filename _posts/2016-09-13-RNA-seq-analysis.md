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
Here, we will create indices for STAR and RSEM

#### STAR

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

#### RESM

#### Usage
{% highlight bash %}
rsem-prepare-reference [options] reference_fasta_file(s) reference_name
{% endhighlight %}
#### Execute
{% highlight bash %}
mkdir GENOME_data/rsem
rsem-prepare-reference --gtf /work3/LSLNGS2015/GENOME_data/Homo_sapiens.GRCh38.82.gtf \
    /work3/LSLNGS2015/GENOME_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    /work3/LSLNGS2015/GENOME_data/rsem/rsem
{% endhighlight%}
#### Options
{% highlight bash %}
--gtf option specifies path to the gene annotations (in GTF format), and RSEM assumes the FASTA file contains sequence of a genome. If this option is off, RSEM will assume the FASTA file contains the reference transcripts. The name of each sequence in the Multi-FASTA files is its transcript_id.
{% endhighlight %}
### Mapping with STAR (2-pass mode)

#### Execute
{% highlight bash %}
STAR --genomeDir GENOME_data/star --sjdbGTFfile GENOME_data/Homo_sapiens.GRCh38.82.gtf \
    --readFilesIn RNASEQ_data/GM12878.rep1.R1.fastq.gz /RNASEQ_data/GM12878.rep1.R2.fastq.gz \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 \
    --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic \
    --runThreadN 20 --outFileNamePrefix "RNASEQ_data/star_GM12878_rep1/"
{% endhighlight %}
#### Options
{% highlight bash %}
--genomeDir   # path to the directory where genome files are stored.
--sjdbGTFfile # path to the GTF file with annotations.
--readFilesIn # paths to files that contain input read1 (and read2 if PE sequencing).
--readFilesCommand # command line to execute for each of the input file. For example: zcat to uncompress .gz files.
--outSAMtype # type of output, i.e. SAM or BAM.
--outFilterMultimapNmax # max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped。 Default is 10.
--outSAMunmapped # output of unmapped reads in the SAM format, None or Within SAM file.
--quantMode # types of quantification requested, i.e. GeneCounts(output ReadsPerGene.out.tab) or TranscriptomeSAM(output Aligned.toTranscriptome.out.bam)
--twopassMode # 2-pass mapping mode. In the first pass, the novel junctions are detected and inserted into the genome indices. In the second pass, all reads will be re-mapped using annotated (from the GTF file) and novel (detected in the first pass) junctions. While this doubles the run time, it significantly increases sensitivity to novel splice junctions.
--runThreadN # number of threads to run STAR.
--outFileNamePrefix # output files name prefix.
{% endhighlight %}

### Quantification with RSEM
In this tutorial, we use RSEM to quantify the expression of genes ans transcript. In the previous step, we instruct STAR to output genomic alignments in transcriptomic coordinates (i.e. Aligned.toTranscriptome.out.bam). We input this file to RSEM to produce gene and transcript expression levels.

#### Usage
{% highlight bash %}
rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --sam/--bam [--paired-end] input reference_name sample_name
{% endhighlight %}
#### Execute
{% highlight bash %}
rsem-calculate-expression --bam --no-bam-output -p 20 --paired-end --forward-prob 1 \ 
    RNASEQ_data/star_GM12878_rep1/Aligned.toTranscriptome.out.bam GENOME_data/rsem/rsem RNASEQ_data/rsem_GM12878_rep1/rsem >& \ 
    RNASEQ_data/rsem_GM12878_rep1/rsem.log
{% endhighlight %}
#### Options
{% highlight bash %}
--bam  # Input file is in BAM format.
--no-bam-output # Do not output any BAM file.
-p # Number of threads to use.
--paired-end # Input reads are paired-end reads.
--forward-prob # Probability of generating a read from the forward strand of a transcript. 1: strand-specific protocol where all (upstream) reads are derived from the forward strand; 0: strand-specific protocol where all (upstream) read are derived from the reverse strand; 0.5: non-strand-specific protocol.
{% endhighlight%}