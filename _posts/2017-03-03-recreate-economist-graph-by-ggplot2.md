---
layout: post
title: "Recreate Economist graph by ggplot2"
date: 2017-03-03
category: tutorial
tags: [ggplot, Economist, plot, R]
---

Principal Component Analysis ([PCA](http://en.wikipedia.org/wiki/Principal_component)) is a multivariate technique that allows us to summarize the systematic patterns of variations in the data. 

<!--more-->

From a data analysis standpoint, PCA is used for studying one table of observations and variables with the main idea of transforming the observed variables into a set of new variables, the principal components, which are uncorrelated and explain the variation in the data. For this reason, PCA allows to reduce a “complex” data set to a lower dimension in order to reveal the structures or the dominant types of variations in both the observations and the variables.


### Challenge: Recreate This Economist Graph
Graph source: [http://www.economist.com/node/21541178](http://www.economist.com/node/21541178 "economist")

Building off of the graphics you created in the previous exercises, put the finishing touches to make it as close as possible to the original economist graph.

![center](/figures/2017-03-03-recreate-economist-graph-by-ggplot2/Economist1.png) 

### Challenge Solution

Lets start by creating the basic scatter plot, then we can make a list of things that need to be added or changed. The basic plot loogs like this:
{% highlight r %}

dat <- read.csv("/data/EconomistData.csv")

  pc1 <- ggplot(dat, aes(x = CPI, y = HDI, color = Region))
  pc1 + geom_point()
  
{% endhighlight %}
![center](/figures/2017-03-03-recreate-economist-graph-by-ggplot2/pc1.png) 

To complete this graph we need to:

* add a trend line
* change the point shape to open circle
* change the order and labels of Region
* label select points
* fix up the tick marks and labels
* move color legend to the top
* title, label axes, remove legend title
* theme the graph with no vertical guides
* add model R2 (hard)
* add sources note (hard)
* final touches to make it perfect (use image editor for this)

### Adding the trend line

Adding the trend line is not too difficult, though we need to guess at the model being displyed on the graph. A little bit of trial and error leds to

{% highlight r %}
(pc2 <- pc1 +
     geom_smooth(aes(group = 1),
                 method = "lm",
                 formula = y ~ log(x),
                 se = FALSE,
                 color = "red")) +
     geom_point()
{% endhighlight %}
![center](/figures/2017-03-03-recreate-economist-graph-by-ggplot2/pc2.png) 

Notice that we put the geom_line layer first so that it will be plotted underneath the points, as was done on the original graph

### Use open points
This one is a little tricky. We know that we can change the shape with the shape argument, what what value do we set shape to? The example shown in ?shape can help us:

{% highlight r %}
 # A look at all 25 symbols
  df2 <- data.frame(x = 1:5 , y = 1:25, z = 1:25)
  s <- ggplot(df2, aes(x = x, y = y))
  s + geom_point(aes(shape = z), size = 4) + scale_shape_identity()
  # While all symbols have a foreground colour, symbols 19-25 also take a
  # background colour (fill)
  s + geom_point(aes(shape = z), size = 4, colour = "Red") +
    scale_shape_identity()
  s + geom_point(aes(shape = z), size = 4, colour = "Red", fill = "Black") +
    scale_shape_identity()
 {% endhighlight %}
![center](/figures/2015-02-25-PCA-in-R/shapes.png) 

This shows us that shape 1 is an open circle, so

{% highlight r %}
pc2 +
    geom_point(shape = 1, size = 4)
{% endhighlight %}
![center](/figures/2015-02-25-PCA-in-R/pc2.png) 
That is better, but unfortunately the size of the line around the points is much narrower than on the original. This is a frustrating aspect of ggplot2, and we will have to hack around it. One way to do that is to multiple point layers of slightly different sizes.

{% highlight r %}
(pc3 <- pc2 +
     geom_point(size = 4.5, shape = 1) +
     geom_point(size = 4, shape = 1) +
     geom_point(size = 3.5, shape = 1))
{% endhighlight %}
![center](/figures/2015-02-25-PCA-in-R/pc3.png) 

#### Reference
[R graphics tutorials](http://www.grroups.com/blog/r-graphics-tutorial-series-part-6-ggplot2) from Ankit Agarwal
[R graphics ipython notebook](tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.ipynb#Putting-It-All-Together)
