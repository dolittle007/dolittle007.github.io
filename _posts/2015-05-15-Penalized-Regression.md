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


### Data format

We will need to reformat the data into the format each particular software requires. Currently, the genotype file contains one row per individual and one column per marker. The genotypes are coded as 0/1/2 to refer to the number of reference alleles. The phenotype file contains one row per individual, with controls coded as 0 and cases coded as 1.


Take a look at the data files **Genotypes.txt** and  **Phenotypes.txt** (e.g. using the **more** command), and check that you understand how the data is coded.






###Step-by-step instructions

####1. Analysis with glmnet


To perform the analysis, we will need to open an R terminal.  Open up a new terminal window, move to the directory where your files are located, and start R (by typing **R**). To exit R at anytime type ** quit()**.


Now (within R) read in the genotype data by typing:
{% highlight r %} 
geno<-read.table("Genotypes.txt")
{% endhighlight %} 

This command reads the genotypes into a dataframe named "geno". To see the size of the data frame, type:

{% highlight r %} 
dim(geno)
{% endhighlight %} 

The data frame has 2000 rows, one for each individual, and 228 columns, one for each marker. To see the top five lines and first ten columns of this dataframe, type:

{% highlight r %} 
geno[1:5,1:10]
{% endhighlight %} 

Next, we can read in the phenotype file into the dataframe "pheno":

{% highlight r %} 
pheno<-read.table("Phenotypes.txt")
{% endhighlight %} 
To look at the first few lines of the file, type:

{% highlight r %} 
 head(pheno) 
{% endhighlight %} 

This file should have 2000 rows. Before we use penalized regression methods to analyze our data, let\'s look at why these methods may be beneficial. We have previously run single marker analysis (the Armitage Trend test) on our data and saved the p-values in the file "SingleMarkerPvalues.txt". We can open the file in R and plot the -log10 p-values. Since the data were simulated, we know at which marker the causal loci are located.

To read the single marker results into R, type:

{% highlight r %} 
single<-read.table("SingleMarkerPvalues.txt")
{% endhighlight %} 

To plot these, type:

{% highlight r %} 
plot(1:228,-log10(single[,1]))
{% endhighlight %} 

This command will work provided you have an X-windows connection to the Linux server (or if you are running R on your own laptop). Otherwise, to see this plot you will have to save it as a file and transfer it back to your own laptop to view. To save the plot as a pdf file in R (rather than plotting it on the screen), type:

{% highlight r %} 
pdf("TrendTest.pdf")
plot(1:228,-log10(single[,1]))
dev.off()
{% endhighlight %} 


You will have to use a similar set of commands (<tt>pdf("file.pdf")</tt> and <tt>dev.off() </tt>) before and after each plot command,
if you need to save the plot as a file in order to transfer it over to your laptop. To make things easier, we have included these commands before and after each plot command below. If you have an X-windows connection to the Linux server, you can just miss out the commands

{% highlight r %} 
pdf("file.pdf")
dev.off()
{% endhighlight %} 


Is the number and location of the causal loci obvious? We can store the known location of each causal locus in the vector "locus" and draw vertical lines at each causal locus:

{% highlight r %} 
locus<-c(14,46,98,164,176)

pdf("TrendTest2.pdf")
plot(1:228,-log10(single[,1]))
abline(v=locus,col=1:5)
dev.off()
{% endhighlight %} 

The entire region is very significant and it is difficult not only to distinguish which loci are the causal loci but to see that there are 5 different causal loci. Hopefully, our penalized regression methods will select a group of markers that can explain most of the signal.
