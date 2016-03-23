---
layout: post
title: "Penalized Regression Methods for Association Studies"
date: 2015-05-15
category: tutorial
tags: [linear regression model, analysis, logistic regression model, plot, R]
---

The default graphical display of most plotting functions in R is very limited (and usually not very pretty). But that doesn’t mean that we should conform with those crude figures.

<!--more-->
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

```bash 
cd xxxxx
``` 

(where **xxxxx** is replaced by the name of the appropriate folder).


### Data format

We will need to reformat the data into the format each particular software requires. Currently, the genotype file contains one row per individual and one column per marker. The genotypes are coded as 0/1/2 to refer to the number of reference alleles. The phenotype file contains one row per individual, with controls coded as 0 and cases coded as 1.


Take a look at the data files **Genotypes.txt** and  **Phenotypes.txt** (e.g. using the **more** command), and check that you understand how the data is coded.






###Step-by-step instructions

####1. Analysis with glmnet


To perform the analysis, we will need to open an R terminal.  Open up a new terminal window, move to the directory where your files are located, and start R (by typing **R**). To exit R at anytime type ** quit()**.


Now (within R) read in the genotype data by typing:
```r
geno<-read.table("Genotypes.txt")
``` 

This command reads the genotypes into a dataframe named "geno". To see the size of the data frame, type:

```r 
dim(geno)
```

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
```{r}
fit_ridge<-glmnet(geno1,pheno1,family="binomial",alpha=0,nlambda=100)
select_r<-c(5,10,20,30,50,100)

pdf("Ridge.pdf")
par(mfrow=c(3,2))
for(i in 1:6){ 
  plot(abs(fit_ridge$beta[,select_r[i]]),main=paste("Lambda=",fit_ridge$lambda[select_r[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")
  abline(v=locus,col=1:5)
}
dev.off()
```


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

{% highlight r %}
library(grpreg)
{% endhighlight %}

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

Note: The grpreg package has an option to select the lambda satisfying the BIC, AIC or GVC criterion using e.g. 

{% highlight r %}
select(fit_mcp, criterion="BIC")
{% endhighlight %}

Using BIC, the "best" value of lambda seems to be 0.1862722. This fits with the fact that, visually, it seems as if the second and third plots (lambda=0.219 and lambda=0.174) do quite well at picking out a single marker to represent each true causal variant.


####3. Analysis with HyperLasso
The HyperLasso software requires a genotype file with one row per individual and one column per SNP, along with a header row of snp names. While we are still in R, we can create a genotype file named "HlassoGenotypes.txt" in the appropriate format from our current genotype file stored in the variable "geno". We have made a file of SNP names called "SNP_Names.txt" to add as a header line to the "geno" variable: 

{% highlight r %}
snp_names<-read.table("SNPNames.txt")
names(geno)<-snp_names[,1]
write.table(geno,"HlassoGenotypes.txt",col.names=T,row.names=F,quote=F)
{% endhighlight %}

If we now look the the first few lines and columns of our "geno" dataframe (e.g. by typing **geno[1:10,1:10]**), we can see that the column names are now our SNP names.

Make sure to leave your R terminal open so that we can use our previous results to compare the methods.  The HyperLasso is a command line software, so we will need to open a new terminal (i.e. open a new connection to the Linux server) and change into the appropriate directory.

Previously, we analyzed the data over a range of lambda values, but how do we know which lambda is appropriate for model selection? One way to select lambda is to choose that value of lambda that gives the appropriate false positive rate for your data set. This can be done by permuting case/control status (so generating data under the null hypothesis) and determining the number of false positives for each value of lambda. We can repeat this process many times and use the mean, median, or maximum lambda to give us the false positive rate we desire. We have already done this for you and selected <tt>lambda=90.0</tt>. In Linux, type:

{% highlight bash %}
runHLasso -genotypes HlassoGenotypes.txt -target Phenotypes.txt -o HlassoResults.txt -shape 1.0 -lambda 90.0 -std
{% endhighlight %}

The Hyperlasso has two penalty parameters, the shape and the scale parameter. However, setting the scale parameter is equivalent to setting a value for lambda, as the two are related. We have set the shape parameter to 1.0 as recommended by Vignal et. al. (2011). They use the scale parameter to control for the false positive rate, whereas we have chosen to do the same using lambda. 

The output file has a special format so we need use the R command, <tt>dget</tt>, to read in the file. Then we will need to do a little fiddling to extract the information we require. 


First look at the file in Linux:

{% highlight bash %}
more HlassoResults.txt 
{% endhighlight %}

Within R, open the file in R using **dget** and look at it:

{% highlight r %}
hlasso<-dget("HlassoResults.txt")
hlasso<-as.data.frame(hlasso)
hlasso<-as.matrix(hlasso)
hlasso
{% endhighlight %}
The first column records the number of cycles required to find the mode, the log-posterior and the log-likelihood of the mode (the final entry of this column is always 0). 
Subsequent columns describe the covariates selected in the model:

{% highlight r %}
	row 1 - position (marker) of the variable in each of the input files, starting at zero.
	row 2 - the variable name supplied in the input files.
	row 3 - the value of the regression coefficient.
	row 4 - For genetic variables, the genetic model:

			  0 - additive (genotypes coded 0,1,2),
			  1 - recessive (genotypes coded 0,0,1),
			  2 - dominant (genotypes coded 0,1,1),
			  3 - heterozygous (genotypes coded 0,1,0).
{% endhighlight %}

The last column contains the information for the intercept. We want to exclude the first and the last column by creating a vector "exclude" that contains the number 1 and the number of the last column. To extract the information we need, type:

{% highlight r %}
exclude<-c(1,dim(hlasso)[2])
coeff<-as.numeric(hlasso[3,-exclude])

#Remember that the indexing for the first row of marker indexes begins with zero so we must add one to the marker variable "hloci":

hloci<-as.numeric(hlasso[1,-exclude])+1  #INDEXING STARTS AT ZERO 

#Finally, we can plot our results:

pdf("HLasso.pdf")
par(mfrow=c(1,1))
plot(hloci,abs(coeff),xlim=c(1,228),xlab="Abs(Coeff)",ylab="Marker")
abline(v=locus,col=1:5)
dev.off()
{% endhighlight %}

The HyperLasso appears to have similar sparsity to the lasso and the MCP. To see any subtle difference between all the methods, we can plot them together, trying to use the "best" value of lambda for each method:

{% highlight r %}
pdf("AllMethods.pdf")
par(mfrow=c(3,2))

plot(1:228,-log10(single[,1]),main="ATT",xlab="Marker",ylab="-log10(pvalues)")
abline(v=locus,col=1:5)

plot(abs(fit_lasso$beta[,10]),main="Lasso",xlab="Marker",ylab="Abs(Coeff)")
abline(v=locus,col=1:5)

plot(abs(fit_en$beta[,10]),main="Elastic Net",xlab="Marker",ylab="Abs(Coeff)")
abline(v=locus,col=1:5)

plot(abs(fit_ridge$beta[,10]),main="Ridge",xlab="Marker",ylab="Abs(Coeff)")
abline(v=locus,col=1:5)

plot(abs(fit_mcp$beta[-c(1),36]),main="MCP",xlab="Marker",ylab="Abs(Coeff)")
abline(v=locus,col=1:5)

plot(hloci,coeff,xlim=c(1,228),main="HyperLasso",xlab="Abs(Coeff)",ylab="Marker")
abline(v=locus,col=1:5)
dev.off() 
{% endhighlight %}

We can see quite strong similarity between the ATT and the ridge, and between the lasso, MCP, and HyperLasso. The elastic net falls in the middle of these 2 groups.

#### 4. Analysis with PUMA 

WARNING: This last part of the exercise is very slow to run! I've included instructions 
in case you want to try it out at home, but I don't recommend running it during the course.
The most computationally intensive part is the calculation using the NEG penalty (which involves
a 2-dimensional search over penalty parameter values). If you want to implement 
a faster analysis, just use the lasso penalty (i.e. omit the word "NEG" at the end of the
puma command line  below).

The PUMA software requires input files in PLINK's transposed  (.tfam and .tped) format. See 

[http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr)

for details. We have prepared these files (**pumadata.tfam** and  **pumadata.tped**) for you, as this file format is quite different from the format we have been using for the other programs.

To run (penalized) logistic regression in PUMA, assuming lasso and NEG penalties, from the linux command line type:

{% highlight bash %}
puma --tped pumadata.tped --tfam pumadata.tfam --name pumaresults --regression LOGISTIC --penalty LASSO NEG
{% endhighlight %}

Once this has finished running, the lasso results should be in the file **results\_pumaresults\_LASSO.R** and the NEG results should be 
in the file **results\_pumaresults\_NEG.R**


To read in and visualise the  lasso results (for example), use the following R script:

{% highlight r %}
# Read R file
result = dget("results_pumaresults_LASSO.R")

# Use AIC information criterion to pick penalty parameter value
criterion = "AIC"

# Evalaute the number of features selected by each model
modelSize = unlist(sapply(1:length(result), function(i){ncol(result[[i]]$features)}))

# Specify the sample size and apply the 1.5 sqrt(n) rule
n=2000
index = which(modelSize < 1.5 * sqrt(n))

# Identify model with the smallest AIC, and which satisfies the 1.5 sqrt(n) rule
j = which(names(result[[1]]$criteria) == criterion)
criteria = unlist(sapply(1:length(result), function(i){result[[i]]$criteria[j]}))
i = which.min(criteria[index])

# Display P-values and coefficients for best model 
result[[i]]

# Convert into dataframe and merge in SNP numbers 
coeffs<-t(result[[i]]$features)
newcoeffs<-data.frame(row.names(coeffs),coeffs)
names(newcoeffs)<-c("snpname", "beta", "pval") 
snp_names<-read.table("SNPNames.txt")
snps<-data.frame(1:228,snp_names)
names(snps)<-c("snpno","snpname")
final<-merge(newcoeffs,snps, by="snpname", all.x=T, sort=F)

# Display final results
final

# Plot final results 
locus<-c(14,46,98,164,176)
plot(final$snpno, final$beta)
abline(v=locus,col=1:5)
{% endhighlight %}

You can use a similar sequence of commands to read in and visualise the results from PUMA using the
NEG penalty.

### Other Methods

A variety of other methods and software packages exist for performing penalized regression with various different penalty functions.
 Model selection can be done by either cross validation, permutation, goodness of fit or parsimony. 
Other related methods include stepwise regression, forward/backward regression, and subset selection.


***
##References

* Ayers KL, and Cordell HJ (2010) SNP Selection in Genome-Wide and Candidate Gene Studies via Penalized Logistic Regression. Genet Epidemiology 34:879-891.
* Breheny P, Huang J (2009) Penalized methods for bi-level variable selection. Statistics and its interface 2:369-380.
* Friedman J, Hastie T and Tibshirani R (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent. J Statist Software 33:1-22.
* Hoffman GE, Logsdon BA, Mezey JG (2013) PUMA: A unified framework for penalized multiple regression analysis of GWAS data. PLOS Computational Biology in press.
* Hoggart CJ, Whittaker JC, De Iorio M, and Balding DJ (2008) Simultaneous analysis of all SNPs in genome-wide and re-sequencing studies. PLoS Genet 4(7):e10000130.
* Vignal CM, Bansal AT, and Balding DJ (2011) Using Penalised Logistic Regression to Fine Map HLA Variants for Rheumatoid Arthritis. Annals of Human Genet 75:655-664.

[Source](https://www.staff.ncl.ac.uk/heather.cordell/PenalizedRegressionTutorial3.html)
