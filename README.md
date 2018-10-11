# NewtonRaphsonSafe Julia Package

[![pipeline status](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/badges/master/pipeline.svg)](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/commits/master)
[![coverage report](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/badges/master/coverage.svg)](https://gitlab.com/PhageDisplayInference/NewtonRaphsonSafe.jl/commits/master)

Unidimensional root finding, based on rtsafe routine from the book:

Press, William H., Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. Numerical Recipes: The Art of Scientific Computing. 3 edition. Cambridge: Cambridge University Press.

Notably, this implementation allows to pass a *single* function to compute f(x) and f'(x), thus being more efficient if there are computations in common.