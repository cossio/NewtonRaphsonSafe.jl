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

(a0,b0) bracket the root and must satisfy F(a0) * F(b0) ≤ 0.
x0 is an initial guess for the root and must satisfy
    
    x0 ∈ [a0, b0]
"""
function rtsafe(fdf,                    # f(x), f'(x) = fdf(x)
                x0, a0, b0;             # initial root guess and bracket
                maxiter::Integer = 100, # maximum number of iterations
                xatol::Real = 1e-10,    # absolute tolerance in abscisa
                xrtol::Real = 1e-10,    # relative tolerance in abscisa
                fatol::Real = 0,        # absolute tolerance in f(x)
                )

    #= Based on rtsafe, from Press, William H., et al. 2007. Numerical Recipes. 3rd edition. 
    We use false position instead of bissection. At each iteration this routine does the best 
    between a Newton step or a false position step. Thus it is guaranteed to converge to a root. =#

    @assert xatol ≥ 0 && xrtol ≥ 0 && fatol ≥ 0
    @assert maxiter ≥ 0
    @assert a0 ≤ x0 ≤ b0 || b0 ≤ x0 ≤ a0
    @assert isfinite(x0) && isfinite(a0) && isfinite(b0)

    a, b = a0, b0               # current bracket
    x  = xold  = x0             # current and last guess
    dx = dxold = abs(b - a)     # last step size and before last


    call_count = 0  # counts number of fdf calls
    F(x) = begin 
        call_count += 1
        f, df = fdf(x)
        @assert !isnan(f) && !isnan(df)
        f, df
    end

    f, df = F(x)
    fa, dfa = F(a)
    fb, dfb = F(b)

    # bracket condition
    @assert fa ≤ 0 ≤ fb || fb ≤ 0 ≤ fa

    # make sure f(a) ≤ 0
    if fa > 0
        a, b = b, a
        fa, fb = fb, fa
        dfa, dfb = dfb, dfa
    end

    false_count = newton_count = 0

    for iter = 1 : maxiter
        xold  = x
        dxold = dx

        newtondx = f / df
        newtonx = x - newtondx

        if a ≤ newtonx ≤ b && abs(2f) ≤ abs(dxold * df)
            # newton step good to go
            x  = newtonx
            dx = newtondx
            newton_count += 1
        else
            # do false position
            x  = a + (b - a) * fa / (fa - fb)
            dx = (b - a) / 2
            false_count += 1
        end

        f, df = F(x)

        # maintain bracket
        if f < 0
            a = x
        else
            b = x
        end

        # convergence criterion
        if abs(f) ≤ fatol || abs(dx) ≤ xatol || abs(dx) ≤ xrtol * abs(x0) || x == xold
            @debug "rtsafe converged" x dx a b f df call_count false_count newton_count
            @debug "rtsafe convergence reason" xatol xrtol fatol (abs(f) ≤ fatol) (abs(dx) ≤ xatol) (abs(dx) ≤ xrtol * abs(x0))
            return x
        end

    end
    
    @error "Maximum number of iterations exceeded in rtsafe" x0 a0 b0 x a b f df call_count false_count newton_count
    return x
end
