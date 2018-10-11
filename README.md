# NewtonRaphsonSafe Julia Package

[![pipeline status](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/badges/master/pipeline.svg)](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/commits/master)
[![coverage report](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/badges/master/coverage.svg)](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/commits/master)

Unidimensional root finding, based on rtsafe routine from the book:

Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007.
Numerical Recipes: The Art of Scientific Computing. 3 edition. Cambridge: Cambridge University Press.

Notably, this implementation allows you to pass a *single* function that computes `f(x)` and `f'(x)` simultaneously.
If `f(x)` and `f'(x)` do some computations in common, this is more efficient than computing them separately.
