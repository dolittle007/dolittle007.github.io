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



To run glmnet, we first need to open the R library:

{% highlight r %} 
library(glmnet)
{% endhighlight %}

First, we need to convert our genotype data from a dataframe into a matrix by typing:

{% highlight r %} 
geno1=as.matrix(geno)
#and convert our phenotype data into a vector:
pheno1=pheno[,1]
{% endhighlight %}

We can run the lasso, elastic net, and ridge regression in glmnet. Recall that the elastic net penalty is: lambda*{(1-alpha)/2 ||beta||22 + alpha ||beta||1. To run the lasso, we set **alpha=1**.
To analyze a dichotomous outcome such as case/control status we use **family="binomial"**. We will attempt to run the lasso at 100 different values of lambda (the penalty strength) by setting **nlambda=100**. The glmnet default is to standardize the genotypes. We can store our results in the variable "fit_lasso":

{% highlight r %} 
fit_lasso <- glmnet(geno1,pheno1,family="binomial",alpha=1,nlambda=100)

# (This may take a little while to run).
#Once it is finished, to see the headers for the information stored in the dataframe "fit_lasso", type:
names(fit_lasso)
{% endhighlight %} 

This lists all the different pieces of information that are stored in the dataframe "fit_lasso".
We are interested in the coefficients stored in "fit_lasso$beta". To see these, and the values of lambda to which they correspond, type:

{% highlight r %} 
fit_lasso$beta
fit_lasso$lambda
{% endhighlight %} 

Although we asked the program to try 100 different values of lambda, in fact only 99 values were evaluated. This is because glmnet sometimes terminates before nlambda values of lambda have been used, because of numerical instabilities near a saturated fit.


The information in  "fit_lasso$beta" is hard to see as it corresponds to 228 coefficients (one for each marker) at each of 99 values of lambda. You can check this by typing:

{% highlight r %} 
dim(fit_lasso$beta)
{% endhighlight %} 

To plot the coefficents for several (six) different values of lambda indexed by "select_l", adding lines at each causal variant, we type:

{% highlight r %} 
select_l<-c(5,10,20,30,50,90)

pdf("Lasso.pdf")
par(mfrow=c(3,2))
for(i in 1:6){ 
plot(abs(fit_lasso$beta[,select_l[i]]),main=paste("Lambda=",fit_lasso$lambda[select_l[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)") 
abline(v=locus,col=1:5) } 
dev.off()
{% endhighlight %} 

In the lasso, the coefficents of many variables are driven to zero. As you can see, as lambda becomes smaller (our penalty becomes weaker), not only do more variables enter the model, the coefficients become larger (further from zero). It looks as if, for these data, a value of lambda of around 0.02-0.05 does quite well at picking out a single marker to represent each true causal variant.


Next, we can run ridge regression by setting <tt>alpha=0</tt> and plot our results:
{% highlight r %}
fit_ridge<-glmnet(geno1,pheno1,family="binomial",alpha=0,nlambda=100)
select_r<-c(5,10,20,30,50,100)

pdf("Ridge.pdf")
par(mfrow=c(3,2))
for(i in 1:6){ 
  plot(abs(fit_ridge$beta[,select_r[i]]),main=paste("Lambda=",fit_ridge$lambda[select_r[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")
  abline(v=locus,col=1:5)
}
dev.off()
{% endhighlight %} 


Although ridge regression does not perform model selection, for large values of lambda, some of the coefficients have so much shrinkage that we cannot distinguish them from zero. (Note that the y axes in these plots have very different ranges). Notice that as lambda is relaxed, most variables have non-zero coefficients, although they still may be small! Additionally, ridge regression assigns similar coefficients to highly correlated variables rather than selecting one of the group as the lasso does.


We can combine sparsity of lasso regression and the grouping effect of ridge regression by using the elastic net. If we set the mixing parameter **alpha=0.4** and plot:

{% highlight r %}
fit_en<-glmnet(geno1,pheno1,family="binomial",alpha=0.4,nlambda=100)
select_e<-c(2,10,20,35,50,100)

 
pdf("EN.pdf")
par(mfrow=c(3,2))
for(i in 1:6){ 
  plot(abs(fit_en$beta[,select_e[i]]),main=paste("Lambda=",fit_en$lambda[select_e[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")
  abline(v=locus,col=1:5) 
} 
dev.off()
{% endhighlight %} 

We can see that test elastic net allows for clustering, but the clusters are not as dense as the ridge. The second and third plots (lambda=0.1258 and 0.0496)
 give us the clearest idea of the number of independent signals and their approximate locations.
Additionally, one can produce a coefficient profile plot of the coefficient paths for a fitted glmnet object using plot, eg. 


{% highlight r %}
pdf("path.pdf")
par(mfrow=c(1,1))
plot(fit_lasso)
dev.off()
{% endhighlight %} 

However, these plots are very messy when you have a large number of variables since the path is plotted for each variable.

####2. Analysis with grpreg


The package grpreg fit paths for group lasso, group bridge, or group MCP at a grid of values of the penalty parameter lambda for linear or logistic regression models. Recall that MCP is similar to the lasso except that its flat tails apply less shrinkage to larger coefficients. First, we need to load the grpreg package into R:

We can run grpreg on the same genotype matrix and phenotype vector we used for glmnet: 

{% highlight r%}
fit_mcp<-grpreg(geno1, pheno1, family="binomial", penalty="gMCP", nlambda=100, lambda.min=.1)

#(This may take some time).

#Again, we can plot our results for selected lambda values:

select_m<-c(2,10,20,30,50,100)

pdf("MCP.pdf")
par(mfrow=c(3,2))
for(i in 1:6){ 
  plot(abs(fit_mcp$beta[-c(1),select_m[i]]),main=paste("Lambda=",fit_mcp$lambda[select_m[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)") 
  abline(v=locus,col=1:5) } 
dev.off()
{% endhighlight %}
