---
layout: post
title: "Computer Practical Exercise on Penalized Regression Methods for Association Studies"
date: 2015-05-15
category: tutorial
tags: [linear regression model, analysis, logistic regression model, plot, R]
---
 

### Computer Practical Exercise on Penalized Regression Methods for Association Studies

### Overview 

### Purpose

In this exercise you will be carrying out case/control association analysis on a gene region assumed to have previously shown association with disease. 
The purpose is to determine which of the associated loci in the region are likely to be the variants causing disease (or be most associated with those variants causing the disease) using penalized regression methods to perform model selection.
 
### Methodology

We will use penalized regression analysis as implemented in the R package glmnet and grpreg, and the command line software HyperLasso. The glmnet package can perform lasso penalized regression, ridge regression, and the elastic net. The lasso (L1) penalty encourages sparsity while the ridge (L2) penalty encourages highly correlated variables to have similar coefficients. The elastic net is a linear combination of these two penalties. 

The grpreg package allows the user to group variables that are to be encouraged in and out of the model together. However, we will not use this functionality in this tutorial. Instead, we will use the minimax-concave penalty (MCP) provided by the package, which is a penalty that provides similar penalization as the lasso but has flat tails to give constant penalization after a user defined threshold. The HyperLasso uses a Bayesian inspired NEG penalty, which imposes heavy shrinkage on small coefficients, and relatively less shrinkage on large coefficients.



### Program documentation

####glmnet:

The R package glmnet has R documentation including a pdf manual, and can be downloaded at:
[http://cran.r-project.org/web/packages/glmnet/index.html](http://cran.r-project.org/web/packages/glmnet/index.html)

#### grpreg
The R package grpreg has R documentation including a pdf manual, and can be downloaded at:
[http://cran.r-project.org/web/packages/grpreg/index.html](http://cran.r-project.org/web/packages/grpreg/index.html)


#### R documentation:

The R website is at
[http://www.r-project.org/](http://www.r-project.org/) 

From within R, one can obtain help on any command **xxxx** by typing '**help(xxxx)**'


####HyperLasso:
Full documentation can be found at:

[http://www.ebi.ac.uk/projects/BARGEN/] (http://www.ebi.ac.uk/projects/BARGEN/)



 #### PUMA

The PUMA software can be found at:
 [http://mezeylab.cb.bscb.cornell.edu/Software.aspx](http://mezeylab.cb.bscb.cornell.edu/Software.aspx)


### Data overview


We will be analyzing simulated data consisting of 228 SNP markers from the CTLA4 gene region. We have 1000 cases and 1000 controls. Since this is simulated data, we know where the true causal loci are located.




### Appropriate data

Appropriate data for this exercise is genotype data from 
a region of interest, such as fine mapping, re-sequencing, or imputation data, typed in a number of unrelated individuals. These analyses can be done for either a dichotomous or a quantitative trait. Here we will consider case/control data.

***
## Instructions


### Data files

You should ensure you have the following files saved to an appropriate directory (folder) on your machine:

[Genotypes.txt](/data/Genotypes.txt)

[Phenotypes.txt](/data/Phenotypes.txt)

[SingleMarkerPvalues.txt](/data/SingleMarkerPvalues.txt)

[SNPNames.txt](/data/SNPNames.txt)

[pumadata.tfam](/data/pumadata.tfam)

[pumadata.tped](/data/pumadata.tped)



Within Linux, open a terminal window, and move into the directory where the data files are e.g. by typing 

{% highlight bash %} 
cd xxxxx
{% endhighlight %} 

(where **xxxxx** is replaced by the name of the appropriate folder).
