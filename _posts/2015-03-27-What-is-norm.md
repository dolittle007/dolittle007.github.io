---
layout: post
title: "What is norm?"
date: 2015-03-27
category: opinion
---

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

$\\|x\\|\_0 = \sqrt[0]{\sum\_{i}x_i^0$

This norm is a bit tricky because there is a present of zeroth-power and zeroth-root in it. Obviously any x > 0 will become one, but the problems of the definition of zeroth-power and especially zeroth-root is messing things around here. So in reality, most mathematicians and engineers use this definition of l0-norm instead:

$\|x\|\_0 = \#(i | x_i \neq 0)$

that is a total number of non-zero elements in a vector.

Because it is a number of non-zero element, there is so many applications that use l0-norm. Lately it is even more in focus because of the rise of the Compressive Sensing scheme, which is try to find the sparsest solution of the under-determined linear system. The sparsest solution means the solution which has fewest non-zero entries, i.e. the lowest l0-norm. This problem is usually regarding as a optimisation problem of l0-norm or l0-optimisation.
