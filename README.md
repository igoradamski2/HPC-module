# High Performance Computing Module
This repository contains code for assessment of the M3C High Performance Computing module, taught in 2017 by Prasun Ray at Imperial College London. 

Below are listed some of the more useful and insightful parts of the work, with some graphs to make things interesting.

## Optimization (Newton's method, LBFGS and bracket descent)

In the [cw2/hw2.f90](cw2/hw2.f90) file you can find simple implementations of the Newton optimisation algorithm (second order method) and an implementation of a simple bracket descent algorithm. Both these methods require full knowledge of the cost function (for our toy example given in [cw2/cost.f90](cw2/cost.f90)). To use as a module in python compile the hw2.f90 file using [f2py](https://docs.scipy.org/doc/numpy/f2py), so that a hw2mod.cpython-36m-darwin.so file is created.

The bracket descent algorithm works by randomly placing a triangle on the initial guess point and recentering it on the vertex with the highest/lowest value of the cost function it wants to optimize. It requires a lot of steps to converge but the Fortran implementation makes it extremely quick. Here's graphs comparing it to LBFGS and Newton:

LBFGS vs Bracket Descent  |  Newton vs Bracket Descent
:-------------------------:|:-------------------------:
![](cw2/figs/hw245.png =250x)  |  ![](cw2/figs/hw231.png =250x)

![Optimizer comparison](cw2/figs/hw245.png)



