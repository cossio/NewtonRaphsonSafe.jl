export rtnewton


"""
Solve f(x) == 0.

If r is the actual root, and x the returned value,

    abs(x - r) ≤ max(xatol, xrtol * abs(x))

or

    abs(f(x)) ≤ fatol

The function is specified by the argument fdf, which computes
the objective function and its derivative in a single call, as:

    f(x), f'(x) = fdf(x)

x0 is an initial guess for the root.
"""
function rtnewton(fdf,                    # f(x), f'(x) = fdf(x)
                  x0;                     # initial root guess
                  maxiter::Integer = 100, # maximum number of iterations
                  xatol::Real = 1e-10,    # absolute tolerance in abscisa
                  xrtol::Real = 1e-10,    # relative tolerance in abscisa
                  fatol::Real = 0         # absolute tolerance in f(x)
                  )

    #= Based on rtsafe, from Press, William H., et al. 2007. Numerical Recipes. 3rd edition. 
    At each iteration this routine does the best between a Newton step or a bissection step.
    Thus it is guaranteed to converge to a root. =#

    if !(xatol ≥ 0 && xrtol ≥ 0 && fatol ≥ 0)
        throw(ArgumentError("tolerances must be non-negative; got xatol=$xatol, xrtol=$xrtol, fatol=$fatol"))
    elseif maxiter < 0
        throw(ArgumentError("maxiter must be a non-negative integer; got maxiter=$maxiter"))
    elseif !isfinite(x0)
        throw(ArgumentError("x0 must be finite; got x0=$x0"))
    end

    x = x0 # current root guess
    dxold = Inf     # step size before last
    dx = dxold      # last step size

    f, df = fdf(x)

    if isnan(f) || isnan(df)
        @error "f or df NaN" f df x
    end

    for iter = 1 : maxiter
            
        dxold = dx
        dx = f / df
        xold = x
        x -= dx

        #@info "Newton with" x0 xl xh dx

        xold == x && return x;
    
        if abs(dx) ≤ xatol || abs(dx) ≤ xrtol * abs(x)   # convergence criterion
            return x
        end
        
        f, df = fdf(x)

        if isnan(f) || isnan(df)
            @error "f or df NaN" f df x
        elseif abs(f) ≤ fatol 
            return x
        end
    end
    
    @warn "Maximum number of iterations exceeded in rtnewton" x dx dxold x0
    
    return x
end
