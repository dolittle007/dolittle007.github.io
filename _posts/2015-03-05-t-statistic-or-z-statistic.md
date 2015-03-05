---
layout: post
title: "When is a z-statistic and when a t-statistic used?"
date: 2015-03-05
category: opinion
tags: [linear regression model, analysis, logistic regression model, plot, R]
---

### When is a z-statistic and when a t-statistic used?

The choice between a z-value or a t-value depends on how the standard error of the coefficients has been calculated. Because the Wald statistic is asymptotically distributed as a standard normal distribution, we can use the z-score to calculate the p-value. When we, in addition to the coefficients, also have to estimate the residual variance, a t-value is used instead of the z-value. In ordinary least squares (OLS, normal linear regression), the variance-covariance matrix of the coefficients is \\( Var[\hat{\beta}|X]=\sigma^2(X^{\prime}X)^{-1} \\) where \\( \sigma^2 \\) is the variance of the residuals (which is unknown and has to be estimated from the data) and \\( X \\) is the design matrix. The standard errors of the coefficients are the square roots of the diagonal elements of the variance-covariance matrix. Because we don't know σ2, we have to replace it by its estimate \\( \hat{\sigma}^2=s^2 \\), so, \\( \widehat{se}(\hat{\beta_{j}}) \\) = \\( \sqrt{s^2(X'X)_{jj}^{-1}} \\). Now that's the point: Because we have to estimate the variance of the residuals to calculate the standard error of the coefficients, we need to use a t-value and the t-distribution.

In logistic (and poisson) regression, the variance of the residuals is related to the mean. If \\( Y∼Bin(n,p) \\), the mean is \\( E(Y)=np \\) and the variance is \\( Var(Y)=np(1−p) \\) so the variance and the mean are related. In logistic and poisson regression but not in regression with gaussian errors, we know the expected variance and don't have to estimate it separately. The dispersion parameter \\( \phi \\) indicates if we have more or less than the expected variance. If \\( \phi=1 \\) this means we observe the expected amount of variance, whereas \\( \phi < 1 \\) means that we have less than the expected variance called underdispersion and \\( \phi > 1 \\) means that we have extra variance beyond the expected (called overdispersion). The dispersion parameter in logistic and poisson regression is fixed at 1 which means that we can use the z-score. The dispersion parameter . In other regression types such as normal linear regression, we have to estimate the residual variance and thus, a t-value is used for calculating the p-values. In R, look at these two examples:

### Logistic regression

{% highlight r %}
mydata <- read.csv("http://databeauty.com/data/binary.csv")

mydata$rank <- factor(mydata$rank)

fit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")

summary(fit)

Coefficients:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept) -3.989979   1.139951  -3.500 0.000465 ***
gre          0.002264   0.001094   2.070 0.038465 *  
gpa          0.804038   0.331819   2.423 0.015388 *  
rank2       -0.675443   0.316490  -2.134 0.032829 *  
rank3       -1.340204   0.345306  -3.881 0.000104 ***
rank4       -1.551464   0.417832  -3.713 0.000205 ***
   ---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

(Dispersion parameter for binomial family taken to be 1)

{% endhighlight %}

Note that the dispersion parameter is fixed at 1 and thus, we get z-values.

### Normal linear regression (OLS)

{% highlight r %}
summary(lm(Fertility~., data=swiss))

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)      66.91518   10.70604   6.250 1.91e-07 ***
Agriculture      -0.17211    0.07030  -2.448  0.01873 *  
Examination      -0.25801    0.25388  -1.016  0.31546    
Education        -0.87094    0.18303  -4.758 2.43e-05 ***
Catholic          0.10412    0.03526   2.953  0.00519 ** 
Infant.Mortality  1.07705    0.38172   2.822  0.00734 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7.165 on 41 degrees of freedom

{% endhighlight %}

Here, we have to estimate the residual variance (denoted as "Residual standard error") and hence, we use t-statistic instead of z-statistic.
Of course, in large samples, the t-distribution approximates the normal distribution and the difference doesn't matter.
