# NewtonRaphsonSafe Julia Package

[![Build Status](https://travis-ci.org/cossio/NewtonRaphsonSafe.jl.svg?branch=master)](https://travis-ci.org/cossio/NewtonRaphsonSafe.jl)
[![Coverage Status](https://coveralls.io/repos/github/cossio/NewtonRaphsonSafe.jl/badge.svg?branch=master)](https://coveralls.io/github/cossio/NewtonRaphsonSafe.jl?branch=master)


Uni-dimensional root finding, based on rtsafe routine from the book:

Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007.
Numerical Recipes: The Art of Scientific Computing. 3 edition. Cambridge: Cambridge University Press.

Notably, this implementation allows you to pass a *single* function that computes `f(x)` and `f'(x)` simultaneously.
If `f(x)` and `f'(x)` do some computations in common, this is more efficient than computing them separately.

We also admit an initial guess for the root that is independent of the bracket.