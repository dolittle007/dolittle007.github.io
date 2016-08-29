---
layout: post
title: "RNA-seq normalization"
date: 2016-08-29
category: opinion
tags: [within-sample, between-sample, normalization, RNA-seq]
---

RNA-seq Normalization

<!--more-->

From a data analysis standpoint, PCA is used for studying one table of observations and variables with the main idea of transforming the observed variables into a set of new variables, the principal components, which are uncorrelated and explain the variation in the data. For this reason, PCA allows to reduce a “complex” data set to a lower dimension in order to reveal the structures or the dominant types of variations in both the observations and the variables.


### Preliminaries

Throughout this post “read” refers to both single-end or paired-end reads. The concept of counting is the same with either type of read, as each read represents a fragment that was sequenced.

When saying “feature”, I’m referring to an expression feature, by which I mean a genomic region containing a sequence that can normally appear in an RNA-Seq experiment (e.g. gene, isoform, exon).

Finally, I use the random variable \\(X_i\\) to denote the counts you observe from a feature of interest i. Unfortunately, with alternative splicing you do not directly observe \\(X_i\\), so often \\(E[X_i]\\) is used, which is estimated using the EM algorithm by a method like [eXpress](http://bio.math.berkeley.edu/express/), [RSEM](http://deweylab.biostat.wisc.edu/rsem/), [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/), [Cufflinks](http://cufflinks.cbcb.umd.edu/), or one of many other tools.

### Counts
“Counts” usually refers to the number of reads that align to a particular feature. I’ll refer to counts by the random variable X_i. These numbers are heavily dependent on two things: (1) the amount of fragments you sequenced (this is related to relative abundances) and (2) the length of the feature, or more appropriately, the effective length. Effective length refers to the number of possible start sites a feature could have generated a fragment of that particular length. In practice, the effective length is usually computed as:

\\( {\widetilde{l}_i} = l_i - \mu_{FLD} + 1 \\),

where $\mu_{FLD}$ is the mean of the fragment length distribution which was learned from the aligned read. If the abundance estimation method you’re using incorporates sequence bias modeling (such as eXpress or Cufflinks), the bias is often incorporated into the effective length by making the feature shorter or longer depending on the effect of the bias.

Since counts are NOT scaled by the length of the feature, all units in this category are not comparable within a sample without adjusting for the feature length. This means you can’t sum the counts over a set of features to get the expression of that set (e.g. you can’t sum isoform counts to get gene counts).

Counts are often used by differential expression methods since they are naturally represented by a counting model, such as a negative binomial (NB2).

####Effective counts

When eXpress came out, they began reporting “effective counts.” This is basically the same thing as standard counts, with the difference being that they are adjusted for the amount of bias in the experiment. To compute effective counts:

\\( \text{effCounts}_i = X_i \cdot \dfrac{l_i}{\widetilde{l}_i} \\).

The intuition here is that if the effective length is much shorter than the actual length, then in an experiment with no bias you would expect to see more counts. Thus, the effective counts are scaling the observed counts up.

####Counts per million

Counts per million (CPM) mapped reads are counts scaled by the number of fragments you sequenced (N) times one million. This unit is related to the FPKM without length normalization and a factor of 10^3:

\\( \text{CPM}_i = \dfrac{X_i}{\dfrac{N}{10^6}} = \dfrac{X_i}{N}\cdot 10^6 \\) 

I’m not sure where this unit first appeared, but I’ve seen it used with edgeR and talked about briefly in the limma voom paper.
