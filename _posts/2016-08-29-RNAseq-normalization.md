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


### PCA in R

In R, there are several functions from different packages that allow us to perform PCA. In this post I'll show you 5 different ways to do a PCA using the following functions (with their corresponding packages in parentheses):
