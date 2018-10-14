export rtsafe


"""
Solve f(x) == 0.

The returned x is such that, if r is the actual root, then:

    abs(x - r) ≤ max(xatol, xrtol * abs(x))

or

    abs(f(x)) ≤ fatol

The function is specified by the argument fdf, which computes
the objective function and its derivative in a single call, as:

    f(x), f'(x) = fdf(x)

xl, xh bracket the root and must satisfy F(xl) ≤ 0 ≤ F(xh).
x0 is an initial guess for the root and must satisfy min(xl,xh) ≤ x0 ≤ max(xh,xl).
"""
function rtsafe(fdf,                    # f(x), f'(x) = fdf(x)
                x0, xl, xh;             # initial root guess and bracket
                maxiter::Integer = 100, # maximum number of iterations
                xatol::Real = 1e-10,    # absolute tolerance in abscisa
                xrtol::Real = 1e-10,    # relative tolerance in abscisa
                fatol::Real = 0,        # absolute tolerance in f(x)
                )

    #= Based on rtsafe, from Press, William H., et al. 2007. Numerical Recipes. 3rd edition. 
    At each iteration this routine does the best between a Newton step or a bissection step.
    Thus it is guaranteed to converge to a root. =#

    x00, xl0, xh0 = x0, xl, xh # save for logging purposes

    if !(xatol ≥ 0 && xrtol ≥ 0 && fatol ≥ 0)
        throw(ArgumentError("tolerances must be non-negative; got xatol=$xatol, xrtol=$xrtol, fatol=$fatol"))
    end

    if maxiter < 0
        throw(ArgumentError("maxiter must be a non-negative integer; got maxiter=$maxiter"))
    end

    if !(xl ≤ x0 ≤ xh) && !(xh ≤ x0 ≤ xl) || !isfinite(x0)
        throw(ArgumentError("x0 must be finite and ∈ [min(xl,xh), max(xl,xh)]; got x0=$x0, xl=$xl, xh=$xh"))
    end

    if !isfinite(xl) || !isfinite(xh)
        throw(ArgumentError("xl = $xl or xh = $xh is not finite"))
    end


    dxold = abs(xh - xl)    # step size before last
    dx = dxold              # last step sized

    @debug begin
        flo, _ = fdf(xl)
        fhi, _ = fdf(xh)
        if !(flo ≤ 0 ≤ fhi)
            @error "invalid bracket" flo fhi xl xh
            error()
        end
        "bracket ok"
    end
    
    f, df = fdf(x0)

    if isnan(f) || isnan(df)
        error("f, df must not be NaN; got f=$f, df=$f at x=$x0")
    end

    for iter = 1 : maxiter
        if (((x0 - xh) * df - f) * ((x0 - xl) * df - f) > 0.0) || (abs(2f) > abs(dxold*df))
            # do bissection
            dxold = dx
			dx = (xh - xl)/2
            x0 = xl + dx
            #@info "did Bissection" x0 xl xh dx
            xl == x0 && return x0     # change in root is negligible
        else
            # Newton step acceptable; take it
            dxold = dx
			dx = f / df
			temp = x0
			x0 -= dx
            #@info "did Newton" x0 xl xh dx
			temp == x0 && return x0;
        end
        
        if abs(dx) ≤ xatol || abs(dx) ≤ xrtol * abs(x0)   # convergence criterion
            return x0
        end
        
        f, df = fdf(x0)

        abs(f) ≤ fatol && return x0
        
        # maintain the bracket
        if f < 0
            xl = x0
        else
            xh = x0
        end
    end
    
    @warn "Maximum number of iterations exceeded in rtsafe" x00 xl0 xh0 x0 xl xh
    
    return x0
end
