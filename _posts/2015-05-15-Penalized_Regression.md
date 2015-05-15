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

<a href="/data/pumadata.tped" target="mainwindow">pumadata.tped</a>



Within Linux, open a terminal window, and move into the directory where the data files are e.g. by typing 

{% highlight linux %} 
cd xxxxx
{% endhighlight %} 

(where **xxxxx** is replaced by the name of the appropriate folder).


### Data format

We will need to reformat the data into the format each particular software requires. Currently, the genotype file contains one row per individual and one column per marker. The genotypes are coded as 0/1/2 to refer to the number of reference alleles. The phenotype file contains one row per individual, with controls coded as 0 and cases coded as 1.


Take a look at the data files **Genotypes.txt** and  **Phenotypes.txt** (e.g. using the **more** command), and check that you understand how the data is coded.






	<h2>Step-by-step instructions</h2>

        <h3>1. Analysis with glmnet </h3>


<p>


To perform the analysis, we will need to open an R terminal.  Open up a new terminal window, move to the directory where your files are located, and start R (by typing <tt>R</tt>). To exit R at anytime type <tt> quit() </tt>.


<br><br>

Now (within R) read in the genotype data by typing:

<br> <br>
<tt>geno<-read.table("Genotypes.txt") </tt>
<br> <br>

This command reads the genotypes into a dataframe named "geno". To see the size of the data frame, type:

<br> <br>
<tt>dim(geno) </tt>
<br> <br>

The data frame has 2000 rows, one for each individual, and 228 columns, one for each marker. To see the top five lines and first ten columns of this dataframe, type:

<br> <br>
<tt>geno[1:5,1:10] </tt>
<br> <br>

Next, we can read in the phenotype file into the dataframe "pheno":

<br> <br>
<tt>pheno<-read.table("Phenotypes.txt") </tt>
<br> <br>

To look at the first few lines of the file, type:

<br> <br>
<tt> head(pheno) </tt>
<br> <br>

This file should have 2000 rows. Before we use penalized regression methods to analyze our data, let's look at why these methods may be beneficial. We have previously run single marker analysis (the Armitage Trend test) on our data and saved the p-values in the file "SingleMarkerPvalues.txt". We can open the file in R and plot the -log10 p-values. Since the data were simulated, we know at which marker the causal loci are located.
<br> <br>
To read the single marker results into R, type:

<br> <br>
<tt>single<-read.table("SingleMarkerPvalues.txt") </tt>
<br> <br>

To plot these, type:

<br> <br>
<tt>plot(1:228,-log10(single[,1]))</tt>
<br> <br>

This command will work provided you have an X-windows connection to the Linux server (or if you are running R on your own laptop). Otherwise, to see this plot you will have to save it as a file and transfer it back to your own laptop to view. To save the plot as a pdf file in R (rather than plotting it on the screen), type:

<br> <br>
<tt>pdf("TrendTest.pdf")</tt> <br>
<tt>plot(1:228,-log10(single[,1]))</tt> <br>
<tt>dev.off() </tt>
<br> <br>


You will have to use a similar set of commands (<tt>pdf("file.pdf")</tt> and <tt>dev.off() </tt>) before and after each plot command,
if you need to save the plot as a file in order to transfer it over to your laptop. To make things easier, we have included these commands before and after each plot command below. If you have an X-windows connection to the Linux server, you can just miss out the commands
<tt>pdf("file.pdf")</tt> and <tt>dev.off() </tt>.


<br> <br>
Is the number and location of the causal loci obvious? We can store the known location of each causal locus in the vector "locus" and draw vertical lines at each causal locus:

<br> <br>
<tt>locus<-c(14,46,98,164,176) </tt>  <br>
<br> 
<tt>pdf("TrendTest2.pdf")</tt> <br>
<tt>plot(1:228,-log10(single[,1]))</tt> <br> 
<tt>abline(v=locus,col=1:5)</tt> <br>
<tt>dev.off() </tt>
<br> <br>

The entire region is very significant and it is difficult not only to distinguish which loci are the causal loci but to see that there are 5 different causal loci. Hopefully, our penalized regression methods will select a group of markers that can explain most of the signal.

<br> <br>

To run glmnet, we first need to open the R library:


<br> <br>
<!-- <tt>library(glmnet,lib.loc="/home/nkla/Library") </tt> -->
<tt>library(glmnet) </tt>
<br> <br>

First, we need to convert our genotype data from a dataframe into a matrix by typing:

<br> <br>
<tt>geno1=as.matrix(geno)</tt>
<br> <br>
and convert our phenotype data into a vector:
<br> <br>
<tt>pheno1=pheno[,1]</tt>
<br> <br>

We can run the lasso, elastic net, and ridge regression in glmnet. Recall that the elastic net penalty is: lambda*{(1-alpha)/2 ||beta||<SUB>2</SUB><sup>2</sup> + alpha ||beta||<SUB>1</SUB>}. To run the lasso, we set <tt>alpha=1</tt>. To analyze a dichotomous outcome such as case/control status we use <tt> family="binomial"</tt>. We will attempt to run the lasso at 100 different values of lambda (the penalty strength) by setting <tt>nlambda=100</tt>. The glmnet default is to standardize the genotypes. We can store our results in the variable "fit_lasso":

<br> <br>
<tt>  fit_lasso <- glmnet(geno1,pheno1,family="binomial",alpha=1,nlambda=100)</tt>
<br> <br>

(This may take a little while to run).

<br> <br>
Once it is finished, to see the headers for the information stored in the dataframe "fit_lasso", type:

<br> <br>
<tt>  names(fit_lasso)</tt>
<br> <br>

This lists all the different pieces of information that are stored in the dataframe "fit_lasso".
We are interested in the coefficients stored in "fit_lasso$beta". To see these, and the values of lambda to which they correspond, type:

<br> <br>
<tt>  fit_lasso$beta</tt>  <br>
<tt>  fit_lasso$lambda</tt>
<br> <br>

Although we asked the program to try 100 different values of lambda, in fact only 99 values were evaluated. This is because glmnet sometimes terminates before nlambda values of lambda have been used, because of numerical instabilities near a saturated fit.

<br> <br>
The information in  "fit_lasso$beta" is hard to see as it corresponds to 228 coefficients (one for each marker) at each of 99 values of lambda. You can check this by typing:

<br> <br>
<tt>  dim(fit_lasso$beta)</tt> 
<br> <br>

To plot the coefficents for several (six) different values of lambda indexed by "select_l", adding lines at each causal variant, we type:

<br> <br>
<tt> select_l<-c(5,10,20,30,50,90) </tt><br>
<br> 
<tt>pdf("Lasso.pdf")</tt> <br>
<tt> par(mfrow=c(3,2))  </tt><br>
<tt> for(i in 1:6){ </tt><br>
<tt> plot(abs(fit_lasso$beta[,select_l[i]]),main=paste("Lambda=",fit_lasso$lambda[select_l[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")  </tt><br>
<tt> abline(v=locus,col=1:5) } </tt><br>
<tt>dev.off() </tt>
<br> <br>
In the lasso, the coefficents of many variables are driven to zero. As you can see, as lambda becomes smaller (our penalty becomes weaker), not only do more variables enter the model, the coefficients become larger (further from zero). It looks as if, for these data, a value of lambda of around 0.02-0.05 does quite well at picking out a single marker to represent each true causal variant.

<br> <br>
Next, we can run ridge regression by setting <tt>alpha=0</tt> and plot our results:
<br> <br>
<tt>fit_ridge<-glmnet(geno1,pheno1,family="binomial",alpha=0,nlambda=100)</tt><br>
<tt> select_r<-c(5,10,20,30,50,100) </tt><br>


<br> 
<tt>pdf("Ridge.pdf")</tt> <br>
<tt> par(mfrow=c(3,2))  </tt><br>
<tt> for(i in 1:6){ </tt><br>
<tt> plot(abs(fit_ridge$beta[,select_r[i]]),main=paste("Lambda=",fit_ridge$lambda[select_r[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")  </tt><br>
<tt> abline(v=locus,col=1:5) } </tt><br>
<tt>dev.off() </tt>
<br> <br>
Although ridge regression does not perform model selection, for large values of lambda, some of the coefficients have so much shrinkage that we cannot distinguish them from zero. (Note that the y axes in these plots have very different ranges). Notice that as lambda is relaxed, most variables have non-zero coefficients, although they still may be small! Additionally, ridge regression assigns similar coefficients to highly correlated variables rather than selecting one of the group as the lasso does.

<br> <br>
We can combine sparsity of lasso regression and the grouping effect of ridge regression by using the elastic net. If we set the mixing parameter <tt>alpha=0.4</tt> and plot:

<br> <br>
<tt> fit_en<-glmnet(geno1,pheno1,family="binomial",alpha=0.4,nlambda=100)</tt><br>
<tt> select_e<-c(2,10,20,35,50,100) </tt><br>

<br> 
<tt>pdf("EN.pdf")</tt> <br>
<tt> par(mfrow=c(3,2))  </tt><br>
<tt> for(i in 1:6){ </tt><br>
<tt> plot(abs(fit_en$beta[,select_e[i]]),main=paste("Lambda=",fit_en$lambda[select_e[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)")  </tt><br>
<tt> abline(v=locus,col=1:5) } </tt><br>
<tt>dev.off() </tt>
<br> <br>
We can see that test elastic net allows for clustering, but the clusters are not as dense as the ridge. The second and third plots (lambda=0.1258 and 0.0496)
 give us the clearest idea of the number of independent signals and their approximate locations.
Additionally, one can produce a coefficient profile plot of the coefficient paths for a fitted glmnet object using plot, eg. 


<br> <br>
<tt>pdf("path.pdf")</tt> <br>
<tt> par(mfrow=c(1,1))  </tt><br>
<tt>plot(fit_lasso)</tt> <br>
<tt>dev.off() </tt>
<br> <br>


However, these plots are very messy when you have a large number of variables since the path is plotted for each variable.
<br> <br> <br>

<p>

        <h3>2. Analysis with grpreg </h3>
<p>

The package grpreg fit paths for group lasso, group bridge, or group MCP at a grid of values of the penalty parameter lambda for linear or logistic regression models. Recall that MCP is similar to the lasso except that its flat tails apply less shrinkage to larger coefficients. First, we need to load the grpreg package into R:


<br> <br>
<!-- <tt>library(grpreg,lib.loc="/home/nkla/Library") </tt>  -->
<tt>library(grpreg) </tt>
<br> <br>
We can run grpreg on the same genotype matrix and phenotype vector we used for glmnet: 
<br> <br>

<tt>fit_mcp<-grpreg(geno1, pheno1, family="binomial", penalty="gMCP", nlambda=100, lambda.min=.1)</tt>
<br> <br>

(This may take some time).
<br> <br>
Again, we can plot our results for selected lambda values:
<br> <br>
<tt>select_m<-c(2,10,20,30,50,100)</tt><br>


<br> 
<tt>pdf("MCP.pdf")</tt> <br>
<tt>par(mfrow=c(3,2))  </tt><br>
<tt>for(i in 1:6){ </tt><br>
<tt>plot(abs(fit_mcp$beta[-c(1),select_m[i]]),main=paste("Lambda=",fit_mcp$lambda[select_m[i]],sep=""),xlab="Marker",ylab="Abs(Coeff)") </tt><br>
<tt>abline(v=locus,col=1:5) } </tt><br>
<tt>dev.off() </tt>
<br> <br>
Note: The grpreg package has an option to select the lambda satisfying the BIC, AIC or GVC criterion using e.g. 

<br> <br>
<tt>select(fit_mcp, criterion="BIC") <br> </tt>
<br> <br>

Using BIC, the "best" value of lambda seems to be 0.1862722. This fits with the fact that, visually, it seems as if the second and third plots (lambda=0.219 and lambda=0.174) do quite well at picking out a single marker to represent each true causal variant.
<br><br> <br>


<p>
        <h3>3. Analysis with HyperLasso </h3>
<p>

The HyperLasso software requires a genotype file with one row per individual and one column per SNP, along with a header row of snp names. While we are still in R, we can create a genotype file named "HlassoGenotypes.txt" in the appropriate format from our current genotype file stored in the variable "geno". We have made a file of SNP names called "SNP_Names.txt" to add as a header line to the "geno" variable: 

<br> <br>
<tt>snp_names<-read.table("SNPNames.txt") </tt><br>
<tt>names(geno)<-snp_names[,1] </tt><br>
<tt>write.table(geno,"HlassoGenotypes.txt",col.names=T,row.names=F,quote=F) </tt>

<br> <br>
If we now look the the first few lines and columns of our "geno" dataframe (e.g. by typing <tt>geno[1:10,1:10]</tt>), we can see that the column names are now our SNP names.

Make sure to leave your R terminal open so that we can use our previous results to compare the methods.  The HyperLasso is a command line software, so we will need to open a new terminal (i.e. open a new connection to the Linux server) and change into the appropriate directory.
<br>
<br>
Previously, we analyzed the data over a range of lambda values, but how do we know which lambda is appropriate for model selection? One way to select lambda is to choose that value of lambda that gives the appropriate false positive rate for your data set. This can be done by permuting case/control status (so generating data under the null hypothesis) and determining the number of false positives for each value of lambda. We can repeat this process many times and use the mean, median, or maximum lambda to give us the false positive rate we desire. We have already done this for you and selected <tt>lambda=90.0</tt>. In Linux, type:

<br> <br>
<!-- <tt>/home/nkla/Library/HyperLasso/runHLasso -genotypes HlassoGenotypes.txt -target Phenotypes.txt -o HlassoResults.txt -shape 1.0 -lambda 90.0 -std </tt>  -->
<tt>runHLasso -genotypes HlassoGenotypes.txt -target Phenotypes.txt -o HlassoResults.txt -shape 1.0 -lambda 90.0 -std </tt>  
<br> <br>
The Hyperlasso has two penalty parameters, the shape and the scale parameter. However, setting the scale parameter is equivalent to setting a value for lambda, as the two are related. We have set the shape parameter to 1.0 as recommended by Vignal et. al. (2011). They use the scale parameter to control for the false positive rate, whereas we have chosen to do the same using lambda. 

The output file has a special format so we need use the R command, <tt>dget</tt>, to read in the file. Then we will need to do a little fiddling to extract the information we require. 

<br> <br>
First look at the file in Linux:

<br> <br>
<tt> more HlassoResults.txt </tt>
<br> <br>

Within R, open the file in R using <tt>dget</tt> and look at it:

<br> <br>
<tt> hlasso<-dget("HlassoResults.txt") </tt><br>
<tt> hlasso<-as.data.frame(hlasso) </tt><br>
<tt> hlasso<-as.matrix(hlasso) </tt><br>
<tt> hlasso </tt>
<br> <br>
The first column records the number of cycles required to find the mode, the log-posterior and the log-likelihood of the mode (the final entry of this column is always 0). 
Subsequent columns describe the covariates selected in the model:

<br> <br>
</p>


<pre>
	row 1 - position (marker) of the variable in each of the input files, starting at zero.
	row 2 - the variable name supplied in the input files.
	row 3 - the value of the regression coefficient.
	row 4 - For genetic variables, the genetic model:

			  0 - additive (genotypes coded 0,1,2),
			  1 - recessive (genotypes coded 0,0,1),
			  2 - dominant (genotypes coded 0,1,1),
			  3 - heterozygous (genotypes coded 0,1,0).
</pre>

<p>

<br> <br>
The last column contains the information for the intercept. We want to exclude the first and the last column by creating a vector "exclude" that contains the number 1 and the number of the last column. To extract the information we need, type:
<br> <br>
<tt> exclude<-c(1,dim(hlasso)[2])  </tt><br>
<tt> coeff<-as.numeric(hlasso[3,-exclude]) </tt><br>
<br>
Remember that the indexing for the first row of marker indexes begins with zero so we must add one to the marker variable "hloci":
<br>
<br>
<tt> hloci<-as.numeric(hlasso[1,-exclude])+1  #INDEXING STARTS AT ZERO </tt><br>
<br>
Finally, we can plot our results:

<br><br>
<tt>pdf("HLasso.pdf")</tt> <br>
<tt> par(mfrow=c(1,1)) </tt><br>
<tt> plot(hloci,abs(coeff),xlim=c(1,228),xlab="Abs(Coeff)",ylab="Marker") </tt><br>
<tt> abline(v=locus,col=1:5)   </tt><br>
 <tt>dev.off() </tt>
<br> <br>
The HyperLasso appears to have similar sparsity to the lasso and the MCP. To see any subtle difference between all the methods, we can plot them together, trying to use the "best" value of lambda for each method:
<br> <br>

<tt>pdf("AllMethods.pdf")</tt> <br>
<tt>par(mfrow=c(3,2))</tt><br> 

<tt>plot(1:228,-log10(single[,1]),main="ATT",xlab="Marker",ylab="-log10(pvalues)")</tt><br>
<tt>abline(v=locus,col=1:5) </tt><br>

<tt>plot(abs(fit_lasso$beta[,10]),main="Lasso",xlab="Marker",ylab="Abs(Coeff)")</tt><br> 
<tt>abline(v=locus,col=1:5) </tt><br>

<tt>plot(abs(fit_en$beta[,10]),main="Elastic Net",xlab="Marker",ylab="Abs(Coeff)") </tt><br>
<tt>abline(v=locus,col=1:5) </tt><br>

<tt>plot(abs(fit_ridge$beta[,10]),main="Ridge",xlab="Marker",ylab="Abs(Coeff)")</tt><br> 
<tt>abline(v=locus,col=1:5) </tt><br>

<tt>plot(abs(fit_mcp$beta[-c(1),36]),main="MCP",xlab="Marker",ylab="Abs(Coeff)")</tt><br> 
<tt>abline(v=locus,col=1:5) </tt><br>

<tt>plot(hloci,coeff,xlim=c(1,228),main="HyperLasso",xlab="Abs(Coeff)",ylab="Marker")</tt><br>
<tt>abline(v=locus,col=1:5) </tt><br>
 <tt>dev.off() </tt>
<br><br>

We can see quite strong similarity between the ATT and the ridge, and between the lasso, MCP, and HyperLasso. The elastic net falls in the middle of these 2 groups.

<br> <br><br>

<p>
        <h3>4. Analysis with PUMA </h3>
<p>
<br>
WARNING: This last part of the exercise is very slow to run! I've included instructions 
in case you want to try it out at home, but I don't recommend running it during the course.
The most computationally intensive part is the calculation using the NEG penalty (which involves
a 2-dimensional search over penalty parameter values). If you want to implement 
a faster analysis, just use the lasso penalty (i.e. omit the word "NEG" at the end of the
puma command line  below).



<br>
<br> <br>
The PUMA software requires input files in PLINK's transposed  (.tfam and .tped) format. See 

<br><br>
<a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr" target="mainwindow"> http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr </a><br> <br>

for details. We have prepared these files (<tt>pumadata.tfam </tt> and  <tt> pumadata.tped</tt>) for you, as this file format is quite different from the format we have been using for the other programs.

<br> <br>

To run (penalized) logistic regression in PUMA, assuming lasso and NEG penalties, from the linux command line type:


<br><br>
 <tt> puma --tped pumadata.tped --tfam pumadata.tfam --name pumaresults --regression LOGISTIC --penalty LASSO NEG </tt>
<br><br>

Once this has finished running, the lasso results should be in the file <tt> results_pumaresults_LASSO.R </tt> and the NEG results should be 
in the file <tt> results_pumaresults_NEG.R </tt>

<br><br>
To read in and visualise the  lasso results (for example), use the following R script:


<br><br>
<tt># Read R file </tt><br>
<tt>result = dget("results_pumaresults_LASSO.R") </tt><br>
<br>
<tt># Use AIC information criterion to pick penalty parameter value </tt><br>
<tt>criterion = "AIC" </tt><br>
<br>
<tt># Evalaute the number of features selected by each model  </tt><br>
<tt>modelSize = unlist(sapply(1:length(result), function(i){ncol(result[[i]]$features)})) </tt><br>
<br>
<tt># Specify the sample size and apply the 1.5 sqrt(n) rule </tt><br>
<tt>n=2000 </tt><br>
<tt>index = which(modelSize < 1.5 * sqrt(n)) </tt><br>
<br>
<tt># Identify model with the smallest AIC, and which satisfies the 1.5 sqrt(n) rule </tt><br>
<tt>j = which(names(result[[1]]$criteria) == criterion) </tt><br>
<tt>criteria = unlist(sapply(1:length(result), function(i){result[[i]]$criteria[j]})) </tt><br>
<tt> i = which.min(criteria[index]) </tt><br>
<br>
<tt># Display P-values and coefficients for best model </tt><br>
<tt>result[[i]] </tt><br>
<br>
<tt># Convert into dataframe and merge in SNP numbers </tt><br>
<tt>coeffs<-t(result[[i]]$features) </tt><br>
<tt>newcoeffs<-data.frame(row.names(coeffs),coeffs)  </tt><br>
<tt>names(newcoeffs)<-c("snpname", "beta", "pval") </tt><br>
<tt>snp_names<-read.table("SNPNames.txt")</tt><br>
<tt>snps<-data.frame(1:228,snp_names) </tt><br>
<tt>names(snps)<-c("snpno","snpname")</tt><br>
<tt>final<-merge(newcoeffs,snps, by="snpname", all.x=T, sort=F)</tt><br>
<br>
<tt># Display final results </tt><br>
<tt>final </tt><br>
<br>
<tt># Plot final results </tt><br>
<tt>locus<-c(14,46,98,164,176)</tt><br>
<tt>plot(final$snpno, final$beta)</tt><br>
<tt>abline(v=locus,col=1:5)</tt><br>
<br><br>
You can use a similar sequence of commands to read in and visualise the results from PUMA using the
NEG penalty.


</p>
	<h2>Other Methods</h2>

<p> 
A variety of other methods and software packages exist for performing penalized regression with various different penalty functions.
 Model selection can be done by either cross validation, permutation, goodness of fit or parsimony. 
Other related methods include stepwise regression, forward/backward regression, and subset selection.
<br>
<br>
</p>


<hr>
<h1>References</h1>
<p>
Ayers KL, and Cordell HJ (2010) SNP Selection in Genome-Wide and Candidate Gene Studies via Penalized Logistic Regression. Genet Epidemiology 34:879-891.
<br><br>
Breheny P, Huang J (2009) Penalized methods for bi-level variable selection. Statistics and its interface 2:369-380.
<br><br>

Friedman J, Hastie T and Tibshirani R (2008) Regularization Paths for Generalized Linear
Models via Coordinate Descent. J Statist Software 33:1-22.
<br><br>
Hoffman GE, Logsdon BA, Mezey JG (2013) PUMA: A unified framework for penalized multiple regression analysis of GWAS data. PLOS Computational Biology in press.
<br><br>
Hoggart CJ, Whittaker JC, De Iorio M, and Balding DJ (2008) Simultaneous analysis of all SNPs in genome-wide and re-sequencing studies. PLoS Genet 4(7):e10000130.
<br><br>
Vignal CM, Bansal AT, and Balding DJ (2011) Using Penalised Logistic Regression to Fine Map HLA Variants for Rheumatoid Arthritis. Annals of Human Genet 75:655�664.

<br><br>
</p>


<hr>
<address>
Exercises prepared by: Kristin Ayers and Heather Cordell <br>
	     Checked/ammended by: Heather Cordell<br>
Programs used: R, glmnet, grpreg, hyperlasso, PUMA <br>
Last updated:
</address>
</body>
</html>
