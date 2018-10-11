export rtsafe


"""
Solve F(x) == 0.

The returned x is such that, if x0 is the actual root, then

abs(x0 - x) ≤ max(xatol, xrtol * max(abs(x)))

or

abs(f(x)) ≤ max(fatol, frtol * abs(f(x)))

The function is specified by the arguments F, dF, FdF:

f  = F(x),  f' = dF(x),  (f, f') = FdF(x)

x1, x2 bracket the root and must satisfy F(x1) ≤ 0 ≤ F(x2) or F(x1) ≤ 0 ≤ F(x2).
"""
function rtsafe(F, dF, FdF, 
                x1::Real, x2::Real;         # root bracket
                maxiter::Integer = 100,     # maximum number of iterations
                xatol::Real = 1e-10,        # absolute tolerance in abscisa
                xrtol::Real = 1e-10,        # relative tolerance in abscisa
                fatol::Real = 0,            # absolute tolerance in f(x)
                frtol::Real = 0             # relative tolerance in f(x)
                )

    #= Based on rtsafe, from Press, William H., et al. 2007. Numerical Recipes. 3rd edition. =#
    if !(xatol ≥ 0 && xrtol ≥ 0 && fatol ≥ 0 && frtol ≥ 0)
        throw(ArgumentError("tolerances must be non-negative; got xatol=$xatol, xrtol=$xrtol, fatol=$fatol, frtol=$frtol"))
    end
    if maxiter < 0
        throw(ArgumentError("maxiter must be a non-negative integer; got maxiter=$maxiter"))
    end
    if !isfinite(x1) || !isfinite(x2)
        throw(ArgumentError("x1, x2 must be finite; got x1=$x1, x2=$x2"))
    end

    fl = F(x1)
    fh = F(x2)

    if isnan(fl) || isnan(fh)
        error("F(x1), F(x2) must not be NaN; got F(x1)=$fl, F(x2)=$fh")
    end

    if fl > 0 && fh > 0 || fl < 0 && fh < 0
        throw(ArgumentError("F(x1) and F(x2) must satisfy F(x1) ≤ 0 ≤ F(x2) or F(x2) ≤ 0 ≤ F(x1); got F(x1)=$fl, F(x2)=$fh"))
    end

    iszero(fl) && return x1
    iszero(fh) && return x2

    if fl > 0
        # Orient the search so that f(x1) < 0
        xl, xh = x2, x1
        # fl, fh = fh, fl   # no need for this because we don't use fl,fh below
    end

    rts = (x1 + x2)/2       # initial root guess
    dxold = abs(x2 - x1)    # step size before last
    dx = dxold              # last step size

    f, df = FdF(rts)

    if isnan(f) || isnan(df)
        error("f, df must not be NaN; got f=$f, df=$f at x=$rts from FdF")
    end

    for iter = 1 : maxiter
        if (((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (abs(2f) > abs(dxold*df))
            dxold = dx
			dx = (xh - xl)/2
            rts = xl + dx
            xl == rts && return rts     # change in root is negligible
        else        # Newton step acceptable; take it
            dxold = dx
			dx = f / df
			temp = rts
			rts -= dx
			temp == rts && return rts;
        end
        
        if dx ≤ xatol || dx ≤ xrtol * rts   # convergence criterion
            return rts
        end
        
        f, df = FdF(rts)

        f ≤ fatol && return rts
        
        # maintain the bracket
        if f < 0
            xl = rts
        else
            xh = rts
        end
    end

    @warn "Maximum number of iterations exceeded in rtsafe"
    return rts
end

function rtsafe(f, df, a::Real, b::Real; kwargs...)
    fdf(x) = (f(x), df(x))
    rtsafe(f, df, fdf, a, b; kwargs...)
end

function rtsafe(fdf, a::Real, b::Real; kwargs...)
    f(x) = first(fdf(x))
    df(x) = last(fdf(x))
    rtsafe(f, df, fdf, a, b; kwargs...)
end
