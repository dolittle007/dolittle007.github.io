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

Finally, I use the random variable X_i to denote the counts you observe from a feature of interest i. Unfortunately, with alternative splicing you do not directly observe X_i, so often \mathbb E[X_i] is used, which is estimated using the EM algorithm by a method like [eXpress](http://bio.math.berkeley.edu/express/), [RSEM](http://deweylab.biostat.wisc.edu/rsem/), [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/), [Cufflinks](http://cufflinks.cbcb.umd.edu/), or one of many other tools.
