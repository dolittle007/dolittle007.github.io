---
layout: post
title: "Probabilistic Graphical Models for Fraud Detection"
date: 2020-09-15
category: tutorial
tags: [Probabilistic Graphical Models, R]
---

Probabilistic Graphical Models for Fraud Detection.

<!--more-->

Bayesian networks are useful tools for probabilistically computing the interdependencies and outcomes of real-world systems given limited information. Here we describe their use in fraud detection.

In this new technical series I will introduce the basic concepts of Bayesian networks a.k.a. belief networks a.k.a. Probabilistic Graphical Models (PGMs) and explain how they can be used to detect anomalous or fraudulent behaviour in insurance. In particular I will work through specific examples of how we might spot anomalies in medical non-disclosure in the life insurance underwriting process.

PGMs are a very interesting and rigorous way to define the causal relations between different events in a system or process. PGMs allow us to state our prior assumptions and account for new evidence in real-time, thus uncovering patterns and customer behaviours that are not always easily clear.

Before we start, some terms:

- A PGM is a probabilistic model wherein we use a network or graph to express the structure of conditional dependence between random variables i.e if event A happens, then event B will happen with X% probability and so on.
- In particular, a Bayesian network or belief network is a PGM where the the graph is a Directed Acyclic Graph (DAG)
- Variables can be directly observed, e.g. the car battery is Y% charged or unobserved a.ka. latent e.g the customer has fraudulent intent
- Variables in the graph can be categorical: event A happens True or False, color in is Red, Blue or Green, or they can be ordinal, integers or floats. For simplicity, here we will restrict our models to categorical variables only.

## Conditional Probability and Conditional Independence
Before we can discuss Bayesian Networks, we first need to introduce the concepts of conditional probability and conditional independence.

## Conditional Probability

This is simply the probability of observing an event given that we already know the realisation of some other event.

For example, suppose we want to know the probability of seeing a total of $11$ from the roll of two 6-sided dice. There are two ways to score $11$1, so absent any other information, the probability of this event is:

$$ P(T = 11) = \frac{2}{36} = 0.0555 $$

Now suppose we already know that one of the dice has rolled a $5$. Now in order to score $11$ total, there's only one option: to roll a $6$. Thus the probability of the total equalling $11$ is:

$$ P(T = 11 \; | \; D_{1} = 5) = \frac{1}{6} = 0.1667 $$

... where $D_{1}$ is the score of first die.

Knowledge of the outcome of the first die roll $D_{1}$ means we now have a different expectation for the final total $T$.

This concept can be expanded further: instead of having complete information about $D_{1}$, we might have incomplete information and only know that $D_{1}$ is either a $4$ or $5$, with a $50\%$ probability of each. With a PGM we can propagate this knowledge forward to adjust the probabilities for $T$.2

## Conditional Dependence

This is a closely related topic, but more subtle. Let's say we have three variables, $A$, $B$ and $C$, and the outcome of $C$ is (directly) dependent upon $A$ and also (directly) dependent upon $B$. Now we can state that $A$ and $B$ are conditionally dependent based on our knowledge of the outcome $C$.

To illustrate, consider the 3 random variables related to the 2 dice rolls, $T$, $D_{1}$ and $D_{2}$. The total $T$ is dependent on both $D_{1}$ and $D_{2}$, and both dice are physically independent of one another.

If our friend rolls the dice and tells us only the value of $T = 4$ we intuitively know we can make a guess about the dice values, since $T = D_{1} + D_{2}$.

We could make a dependency table based on $T=4$ about the viable combinations (marked X):

|      | D2         |
|      |  1 | 2 | 3 |...
|------|--------------
|D1  1 |  o | o | X |
|    2 |  o | X | o |
|    3 |  X | o | o |
|  ... |    |   |   |
  
Furthermore, if we then learn the value of one die, e.g. $D_{1} = 1$ then we can use our table to deduce that $D_{2} = 3$. The values of the dice are conditionally dependent on the value of their total.

Thinking about systems in this way can be initially counter-intuitive, since here the dice are not physically changed and are truly independent in a physical sense.

More formally, the joint distribution of $A$, $B$, and $C$ can be factorised into two expressions, one containing $A$ and $C$, and the other containing $B$ and $C$, i.e.

$$ P(A, B, C) = f(A, C) \; g(B, C) $$

## Conditional Independence
We see above how 2 independent variables can become conditionally dependent from knowing the value of a 3rd, mutually dependent variable.

Lets say we also have a deck of cards $C$ and draw a single card $C_1$ based on nothing but blind chance. Our friend is interested in the values $T$ and $C_1$ so we report them separately as separate systems, nothing more. I.e. $T$ and $C_1$ are independent, and our primary variables $D_1 + D_{2}$ and $C_1$ are conditionally independent.

If we created a final, super, overarching total $S = D_{1} + D_{2} + C_{1}$ then that would change the graph, and $D_{1}$, $D_{2}$, $C_{1}$ would become conditionally dependent on $S$.

---
## Bayesian Networks and the Sprinkler Network
Now suppose we are trying to create a model with lots of variables. Rather than try to consider all the variables at once, we try to model the system by building up chains of dependency. This allows us to break up the modelling process into smaller chunks, making the whole system much more manageable.

Our Bayesian network or PGM is a directed acyclic graph (DAG). Each variable is a node on the graph, and an edge or vertex represents the conditional relationship between the variables, with the edge pointing towards the dependent variable.

## The Sprinkler Network
One of the simplest Bayesian networks is the Sprinkler network. In this network we have three binary-valued variables:

$R(ain)$, representing whether or not it is raining;
$S(prinkler)$, representing if the lawn sprinkler is turned on; and
$G(rass)$, representing if the grass is wet.
We can't control the rain, but if it is raining we are less likely to need to use the sprinkler, so it is less likely to be turned on. Both rain and an active sprinkler affect the likelihood of the grass being wet.

We can represent all of the above as a DAG:

sprinkler_network

Now that we have the structure in place, we need to think about probabilities. This is generally done by specifying conditional probability tables (CPTs). Each variable has its own CPT, specifying different probabilities for outcomes given the different values of the conditional variables.

In the Sprinkler network:

$R$ is totally independent (no arrows point into it) and so we just need to specify a binomial probability for it.3.
$S$ depends on the value of $R$, and so we need to specify two probabilities, one for each value of $R$.
Finally, $G$ depends on both $S$ and $R$ and so we need 4 probabilities, one for each combination of $S$ and $R$.
We can now fully specify the sprinkler network by specifying the CPTs for each variable. Note the nomenclature $P(A|\neg B)$ means the conditional probability that $A$ happens given that $B$ does not happen.

First the probabilities for $R$:

$$ P(R) = 0.2 $$

Then for $S$:

$$ P(S|R) = 0.01 \;\; P(S|\neg R) = 0.4 $$

And finally for $G$:

$$ P(G|S,R) = 0.99 \;\; P(G|S, \neg R) = 0.9 $$

$$ P(G|\neg S, R) = 0.8 \;\; P(G|\neg S, \neg R) = 0 $$

Using the PGM
Ideally, we can use this network to answer (simple) questions like:

What is the likelihood of the grass being wet?
If the grass is wet, how likely is it to be raining?
If the sprinkler is not on, how likely is it to not be raining?
In order to compute the answers, we need to combine all the CPTs. When the overall size of the network is small, with few nodes (variables) and few edges (outcomes), we can use the 'brute force' method: just multiply out all the conditional tables and calculate the joint distribution of all the variables.

With the full, joint distribution we know everything that can be known about the system: questions are answered by conditioning on the evidence observed and then marginalising out all the other variables.

This is not conceptually complex but can be a laborious accounting exercise, making errors all too common. Thankfully, this is exactly the kind of task computers excel at, so we will use some R packages to handle all the details, allowing us to focus on the interesting stuff.

The package I am going to use from now on is the gRain package, available on CRAN.

## Computing the PGM
We specify the network by first specifying the CPTs using the function cptable(). You may notice dependencies are described using R's formula notation, aiding familiarity. The code used to create this Sprinkler network is shown below, please get in touch if you would like access to the full code repository.
```r
r <- cptable(~rain  
            ,levels = c("yes", "no")
            ,values = c(0.2, 0.8));

s <- cptable(~sprinkler | rain  
            ,levels = c("yes", "no")
            ,values = c(0.01, 0.99, 0.4, 0.6));

g <- cptable(~grass | sprinkler + rain  
            ,levels = c("yes", "no")
            ,values = c(0.99, 0.01, 0.9, 0.1, 0.8, 0.2, 0, 1));

sprinkler.grain <- grain(compileCPT(list(r, s, g)));
```
Once we have the network, it is very straightforward to ask questions:
```r
> ### Q1: What is the likelihood of the grass being wet?
> querygrain(sprinkler.grain, nodes = 'grass')
$grass
 grass
     yes      no
 0.43618 0.56382
 ```
The probability of the grass being wet is $0.436$.
```r
 > ### Q2: If the grass is wet, how likely is it to be raining?
 querygrain(sprinkler.grain, nodes = 'rain', evidence = list(grass = 'yes'))
 $rain
 rain
      yes       no
 0.413086 0.586914
 ```
If the grass is wet, the probability that it is raining is $0.413$.

I will leave answering the third question as an exercise for the reader.

---
## Expanding the Model

Now that we have a simple and working model for the Sprinkler system, we can move on to building a little more complexity. Rather than adding nodes to the system, how do we add additional outcomes to the variables we have?

It makes sense to keep the sprinkler as a binary variable. On/Off is intuitive. Instead, how about we add three levels to the Grass variable, denoted three levels of 'wetness' - 'Dry', 'Damp' and 'Wet'. How does this affect our system?

In terms of specifying conditional probability, we do not need to change a whole lot with our network. $R$ and $S$ remain the same, we still have 4 combinations with which to condition upon, but we now need to specify three values for $G$, one for each of the three levels of $G$.

The code to do this is shown below:
```r
r <- cptable(~rain  
            ,levels = c("rain", "norain")
            ,values = c(0.2, 0.8));

s <- cptable(~sprinkler | rain  
            ,levels = c("on", "off")
            ,values = c(0.01, 0.99, 0.4, 0.6));

g <- cptable(~grass | sprinkler + rain  
            ,levels = c("dry", "damp", "wet")
            ,values = c(0.01, 0.20, 0.79  ## ON, Rain
                       ,0.05, 0.30, 0.65  ## OFF, Rain
                       ,0.10, 0.70, 0.20  ## ON, NoRain
                       ,0.90, 0.09, 0.01  ## OFF, NoRain

expandedsprinkler1.grain <- grain(compileCPT(list(r, s, g)));
```
We can now ask slightly more interesting questions:
```r
 > ### What is the likelihood that it is raining given that the grass is damp?
querygrain(expandedsprinkler1.grain  
           ,nodes='rain'                                   
           ,evidence = list(grass = 'damp'));           
$rain
 rain
     rain   norain
 0.182875 0.817125
 ```
What about the probability of rain given that the grass is not dry (thus either damp or wet)? We cannot calculate the two probabilities and add them, that quantity is not a probability as they are not mutually exclusive events.4

Instead, we need to calculate:

$$ P(R = \text{rain} \; | \; G = \text{damp}) \; P(G = \text{damp}) + P(R = \text{rain} \; | \; G = \text{wet}) \; P(G = \text{wet}) $$

An alternative method - one that is probably easier in this case - is to inspect the joint distribution and add the probabilities.
```r
### What is the likelihood that is raining given that the grass is either
### damp or wet?
 > ftable(querygrain(expandedsprinkler1.grain
                  ,nodes = c("rain", "sprinkler", "grass")
                  ,type = 'joint'),
        row.vars = 'grass');
       rain         rain          norain
       sprinkler      on     off      on     off
 grass
 dry             0.00002 0.00990 0.03200 0.43200
 damp            0.00040 0.05940 0.22400 0.04320
 wet             0.00158 0.12870 0.06400 0.00480
 ```
So, to calculate those probabilities we add up the probabilities where $G$ is either "damp" or "wet" and $R$ is "rain":
```r
       rain         rain          norain
       sprinkler      on     off      on     off
 grass
 dry                   .       .       .       .
 damp            0.00040 0.05940       .       .
 wet             0.00158 0.12870       .       .
 ```
Thus our probability is

$$P(R | G \neg dry) = 0.00040 + 0.00158 + 0.05940 + 0.12870 = 0.19008$$

Finally, what if we also have three levels of $R$? Instead of 'rain' and 'norain', we have 'norain', 'light', and 'heavy'. We need to expand the CPTs, but otherwise the network is unchanged.
```r
### Create the network with three levels for R:
r <- cptable(~rain  
            ,levels = c("norain", "light", "heavy")
            ,values = c(0.7, 0.2, 0.1));

s <- cptable(~sprinkler | rain  
            ,levels = c("on", "off")
            ,values = c(0.4, 0.6, 0.15, 0.85, 0.01, 0.99));

g <- cptable(~grass | sprinkler + rain  
            ,levels = c("dry", "damp", "wet")
            ,values = c(0.10, 0.70, 0.20  ## ON, NoRain
                       ,0.88, 0.10, 0.02  ## OFF, NoRain
                       ,0.05, 0.60, 0.35  ## ON, Light
                       ,0.10, 0.60, 0.30  ## OFF, Light
                       ,0.01, 0.25, 0.74  ## ON, Heavy
                       ,0.05, 0.25, 0.70  ## OFF, Heavy
                        ));

expandedsprinkler2.grain <- grain(compileCPT(list(r, s, g))); 
```

Given that we have wet grass, what is the probability table for the 3 different levels of rain?

```r
 ### What are the probability values for rain given that the
 ### grass is wet?
 > querygrain(expandedsprinkler2.grain
           ,nodes = "rain"
           ,evidence = list(grass = 'wet'))
$rain
 rain
   norain    light    heavy
 0.328672 0.313872 0.357456
```

So, according to our CPTs and the network, wet grass means all three levels of rain are roughly equally likely, not something I would have expected! In this case it would be very helpful to gather evidence on the sprinkler, and include that in the computation.

This is an excellent illustration of why PGMs and Bayesian networks are so useful. Once you start bringing in evidence and condition probability, it is unlikely you will have a good intuition on what the effects are. Human brains are rarely wired well for this kind of problem.

## Building Bayesian Networks from Data
I really like the power of Bayesian networks. Despite being built from simple blocks, the interaction of those probabilities often result in non-intuitive outcomes. They do have limitations, and it is important to consider those.

We will discuss other downsides in the next article, but one big one is probably becoming obvious: constructing the CPTs can be tedious. Bayesian networks often suffer from poor scaling. Creating a CPT for a variable with 5 possible values that depends on three other variables, each with 4 different levels is a nightmare to create. There are $4^{3} = 64$ combinations of conditions and each condition needs a probability for 5 levels.

We can mitigate this a little if we already have a lot of data. That way, we can create the DAG without any probabilities, and then create the network from the data. The grain() function can take data as a parameter and use this to calculate the CPTs.

This saves a lot of time and effort, but has its own issues. Some combinations may not appear in the data, requiring smoothing and introduce errors. More seriously, in many cases Bayesian networks have unobserved, missing, or latent variables, making it much harder to properly estimate the CPTs.

We will return to these issues later in the series.
---
## Summary
In this first article I introduced Bayesian networks, illustrated the concepts using the simple Sprinkler system, and performed some conditional probability calculations using different levels of evidence.

Now that we are comfortable with the basics, we can move on to the primary model I want to discuss in this series: attempting to model medical non-disclosure in underwriting life insurance. This is where customers can fail to disclose pre-existing medical conditions and lifestyle choices when applying for life insurance products. Such non-disclosure is not necessarily fraud (quite often it can be a simple, honest mistake), but insurers are very well-advised to spot it as soon as possible to minimise portfolio risk.

