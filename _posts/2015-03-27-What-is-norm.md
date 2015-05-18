---
layout: post
title: "What is norm?"
date: 2015-03-27
category: opinion
---

All about norm in statistics.

<!--more-->

### What is a norm?

Mathematically a norm is a total size or length of all vectors in a vector space  or matrices. For simplicity, we can say that the higher the norm is, the bigger the (value in) matrix or vector is. Norm may come in many forms and many names, including these popular name: Euclidean distance, Mean-squared Error, etc.

Most of the time you will see the norm appears in a equation like this:

$\\|x\\|$ where x can be a vector or a matrix.

For example, a Euclidean norm of a vector $a = \begin{bmatrix}  3  \\\\  -2  \\\\  1  \end{bmatrix}$ is $\\|a\\|_2=\sqrt{3^2+(-2)^2+1^2}=3.742$ which is the size of vector a

The above example shows how to compute a Euclidean norm, or formally called an l2-norm. There are many other types of norm that beyond our explanation here, actually for every single real number, there is a norm correspond to it (Notice the emphasised word real number, that means it not limited to only integer.)

Formally the ln-norm of x is defined as:

$\\|x\\|\_n = \sqrt[n]{\sum\_{i}|x_i|^n}$ where $n \epsilon \mathbb{R}$

That’s it!  A n-th-root of a summation of all elements to the n-th power is what we call a norm.

The interesting point is even though every ln-norm is all look  very similar to each other, their mathematical properties are very different and thus their application are dramatically different too. Hereby we are going to look into some of these norms in details.



### l0-norm 

The first norm we are going to discuss is a l0-norm. By definition, l0-norm of x is

$\\|x\\|\_0 = \sqrt[0]{\sum\_{i}{x_i}^0}$

This norm is a bit tricky because there is a present of zeroth-power and zeroth-root in it. Obviously any x > 0 will become one, but the problems of the definition of zeroth-power and especially zeroth-root is messing things around here. So in reality, most mathematicians and engineers use this definition of l0-norm instead:

$\|x\|\_0 = \\#(i | x_i \neq 0) $

that is a total number of non-zero elements in a vector.

Because it is a number of non-zero element, there is so many applications that use l0-norm. Lately it is even more in focus because of the rise of the Compressive Sensing scheme, which is try to find the sparsest solution of the under-determined linear system. The sparsest solution means the solution which has fewest non-zero entries, i.e. the lowest l0-norm. This problem is usually regarding as a optimisation problem of l0-norm or l0-optimisation.

### l0-optimisation

Many application, including Compressive Sensing, try to minimise the l0-norm of a vector corresponding to some constraints, hence called “l0-minimisation”. A standard minimisation problem is formulated as:

min \left \| x \right \|_0 subject to Ax = b

However, doing so is not an easy task. Because the lack of l_0-norm’s mathematical representation, l_0-minimisation is regarded by computer scientist as an NP-hard problem, simply says that it’s too complex and almost impossible to solve.

In many case, l_0-minimisation problem is relaxed to be higher-order norm problem such as l_1-minimisation and l_2-minimisation.

### l1-norm

Following the definition of norm, l_1-norm of x is defined as

\left \| x \right \|_1 = \sum_{i} \left | x_i \right |

This norm is quite common among the norm family. It has many name and many forms among various fields, namely Manhattan norm is it’s nickname. If the l_1-norm is computed for a difference between two vectors or matrices, that is

SAD(x_1,x_2) = \left \| x_1-x_2 \right \|_1 = \sum \left | x_{1_i}-x_{2_i} \right |

it is called Sum of Absolute Difference (SAD) among computer vision scientists.

In more general case of signal difference measurement, it may be scaled to a unit vector by:

MAE(x_1,x_2) = \frac{1}{n} \left \| x_1-x_2 \right \|_1 = \frac {1} {n} \sum \left | x_{1_i} - x_{2_i} \right | where n is a size of x.

which is known as Mean-Absolute Error (MAE).

### l2-norm

The most popular of all norm is the l_2-norm. It is used in almost every field of engineering and science as a whole. Following the basic definition, l_2-norm is defined as

\left \| x \right \|_2 = \sqrt{\sum_{i}x_i^2}

l_2-norm is well known as a Euclidean norm, which is used as a standard quantity for measuring a vector difference. As in l_1-norm, if the Euclidean norm is computed for a vector difference, it is known as a Euclidean distance:

\left \| x_1-x_2 \right \|_2 = \sqrt{\sum_i (x_{1_i}-x_{2_i})^2} 

or in its squared form, known as a Sum of Squared Difference (SSD) among Computer Vision scientists:

SSD(x_1,x_2) = \left \| x_1-x_2 \right \|_2^2 = \sum_i (x_{1_i}-x_{2_i})^2

It’s most well known application in the signal processing field is the Mean-Squared Error (MSE) measurement, which is used to compute a similarity, a quality, or a  correlation between two signals. MSE is

MSE(x_1,x_2) = \frac{1}{n} \left \| x_1-x_2 \right \|_2^2 = \frac{1}{n} \sum_i (x_{1_i}-x_{2_i})^2

As previously discussed in l_0-optimisation section, because of many issues from both a computational view and a mathematical view, many l_0-optimisation problems relax themselves to become l_1- and l_2-optimisation instead. Because of this, we will now discuss about the optimisation of l_2.

### l2-optimisation

As in l0-optimisation case, the problem of minimising l2-norm is formulated by

min \left \| x \right \|_2 subject to Ax = b

Assume that the constraint matrix A has full rank, this problem is now a underdertermined system which has infinite solutions. The goal in this case is to draw out the best solution, i.e. has lowest l_2-norm, from these infinitely many solutions. This could be a very tedious work if it was to be computed directly. Luckily it is a mathematical trick that can help us a lot in this work.

By using a trick of Lagrange multipliers, we can then define a Lagrangian

\mathfrak{L}(\boldsymbol{x}) = \left \| \boldsymbol{x} \right \|_2^2+\lambda^{T}(\boldsymbol{Ax}-\boldsymbol{b})

where λ is the introduced Lagrange multipliers. Take derivative of this equation equal to zero to find a optimal solution and get

\hat{\boldsymbol{x}}_{opt} = -\frac{1}{2} \boldsymbol{A}^{T} \lambda

plug this solution into the constraint to get

\boldsymbol{A}\hat{\boldsymbol{x}}_{opt} = -\frac{1}{2}\boldsymbol{AA}^{T}\lambda=\boldsymbol{b}

\lambda=-2(\boldsymbol{AA}^{T})^{-1}\boldsymbol{b}

and finally

\hat{\boldsymbol{x}}_{opt}=\boldsymbol{A}^{T} (\boldsymbol{AA}^{T})^{-1} \boldsymbol{b}=\boldsymbol{A}^{+} \boldsymbol{b}

By using this equation, we can now instantly compute an optimal solution of the L2-optimisation problem. This equation is well known as the Moore-Penrose Pseudoinverse and the problem itself is usually known as Least Square problem, Least Square regression, or Least Square optimisation.

However, even though the solution of Least Square method is easy to compute, it’s not necessary be the best solution. Because of the smooth nature of l_2-norm itself,  it is hard to find a single, best solution for the problem.
